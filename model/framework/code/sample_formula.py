"""Formula-preserving molecule generation, adapted from CoCoGraph's
sample_scripts/sample_molecules_FPSmodel_unseed.py (github.com/manurubo/CoCoGraph).

The source script batches many *different* input molecules across a
ProcessPoolExecutor for large-scale dataset generation. Here we need one
CPU-only, single-process call: one input SMILES in, N generated SMILES out,
all sharing that SMILES's exact molecular formula. This module lifts the
per-molecule denoising loop from `process_batch` (same source file) and
calls the underlying helper functions directly instead of through a process
pool, but performs the same computation, step for step.
"""

import json
import math
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

import torch
from rdkit import Chem, RDLogger
from torch_geometric.data import Data

from lib_functions.adjacency_utils import components_to_graph, nx_to_rdkit
from lib_functions.config import device
from lib_functions.data_preparation_utils import embed_edges_manuel
from lib_functions.formula_utils import build_gt_from_formula
from lib_functions.models import (
    GINEdgeQuadrupletPredictor_MorganFP,
    GINETimePredictor_MorganFP,
)
from lib_functions.sample_utils import calculate_data_molecule_fps, sample_step_graph

RDLogger.DisableLog("rdApp.*")

_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
_CHECKPOINTS_DIR = os.path.abspath(os.path.join(_CODE_DIR, "..", "..", "checkpoints"))
_DATA_DIR = os.path.join(_CODE_DIR, "Data")

_DIFFUSION_CHECKPOINT = os.path.join(
    _CHECKPOINTS_DIR, "FPS_diffusion", "model_epoch_1_slice_22.pth"
)
_TIME_CHECKPOINT = os.path.join(
    _CHECKPOINTS_DIR, "FPS_time", "model_epoch_2_slice_22.pth"
)

# Sigma controls how many double-edge-swap steps are planned for a molecule
# (num_swaps = ceil(sigma * edges)); 0.5 matches the source script's default.
_SIGMA = 0.5
_MAX_GENERATE_RETRIES = 3


def _load_valence_weights(path):
    with open(path, "r") as f:
        raw = json.load(f)
    valid_valences = {}
    valence_weights = {}
    for k, vals in raw.items():
        sym, chs = k.split("__")
        key = (sym, int(chs))
        if isinstance(vals, list):
            valid_valences[key] = set(int(v) for v in vals)
        elif isinstance(vals, dict):
            valid_valences[key] = set(int(v) for v in vals.keys())
            valence_weights[key] = {int(v): float(p) for v, p in vals.items()}
        else:
            valid_valences[key] = set()
    return valence_weights


def _load_models():
    model = GINEdgeQuadrupletPredictor_MorganFP()
    time_model = GINETimePredictor_MorganFP()

    checkpoint = torch.load(_DIFFUSION_CHECKPOINT, map_location=device)
    checkpoint_time = torch.load(_TIME_CHECKPOINT, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])
    time_model.load_state_dict(checkpoint_time["model_state_dict"])

    model = model.to(device)
    time_model = time_model.to(device)
    model.eval()
    time_model.eval()
    return model, time_model


_VALENCE_WEIGHTS = _load_valence_weights(os.path.join(_DATA_DIR, "valid_valences.json"))
with open(os.path.join(_DATA_DIR, "charge_symbol_weights.json"), "r") as _f:
    _CHARGE_SYMBOL_WEIGHTS = json.load(_f)
with open(os.path.join(_DATA_DIR, "radical_symbol_weights.json"), "r") as _f:
    _RADICAL_SYMBOL_WEIGHTS = json.load(_f)

_MODEL, _TIME_MODEL = _load_models()


def _generate_one(formula: str) -> str:
    """Generates a single new molecule sharing the given molecular formula."""

    g_ruido = build_gt_from_formula(
        formula,
        randomize_swaps=0,
        valence_weights=_VALENCE_WEIGHTS,
        charge_symbol_weights=_CHARGE_SYMBOL_WEIGHTS,
        max_sampling_retries=100,
        allow_radicals=(_RADICAL_SYMBOL_WEIGHTS is not None),
        radical_weights=_RADICAL_SYMBOL_WEIGHTS,
    )

    edges_noisy = g_ruido.number_of_edges()
    num_swaps = math.ceil(_SIGMA * edges_noisy)

    tensor, _, _ = embed_edges_manuel(g_ruido, list(g_ruido.nodes()))

    all_smiles_seen = set()
    prediction_time = 0.5
    best_time = 0.5
    best_tensor = tensor.clone()

    for step in range(num_swaps):
        (
            processed_graph,
            tensor_out,
            mol,
            gemb,
            nemb,
            distances,
            edge_index,
            edge_attr,
            dosd_positions,
            componentes_ant,
            fingerprint,
        ) = calculate_data_molecule_fps(g_ruido, tensor, num_swaps, step)

        d = Data(
            x=nemb,
            edge_index=edge_index,
            y=tensor_out,
            xA=gemb,
            edge_attr=edge_attr,
            noiselevel=torch.tensor(prediction_time, device=device),
            distances=torch.Tensor(distances),
            dosd_distances=dosd_positions,
            morgan_fp=fingerprint,
        ).to(device)

        with torch.no_grad():
            _, _, probs_quadrupletas_mod = _MODEL(d)
            time_pred = _TIME_MODEL(d)
        prediction_time = time_pred.detach().cpu().item()

        if step == 0 or prediction_time < best_time:
            best_time = prediction_time
            best_tensor = tensor.clone()

        tensor, smiles_next = sample_step_graph(
            g_ruido, tensor.clone(), probs_quadrupletas_mod.detach(),
            all_smiles_seen, 0, num_swaps, step,
        )
        if smiles_next is None:
            # No valid swap found this step; keep the current tensor and continue.
            continue

    g_gen = components_to_graph(g_ruido.nodes(data=True), best_tensor)
    mol_des = nx_to_rdkit(g_gen, False)
    return Chem.MolToSmiles(mol_des)


def _generate_one_with_retries(formula: str) -> str:
    """Generates one molecule for `formula`, retrying on failure before giving up.

    Module-level (not a closure) so it can be pickled and dispatched to worker
    processes by reference.
    """
    generated = ""
    for _attempt in range(_MAX_GENERATE_RETRIES):
        try:
            generated = _generate_one(formula)
            break
        except Exception:
            continue
    return generated


def _worker_init():
    # Each worker runs one generation at a time; without this, every worker
    # would try to use all cores for its own PyTorch ops, oversubscribing the
    # machine across `n_jobs` processes at once.
    torch.set_num_threads(1)


_POOL = None
_POOL_N_JOBS = None


def _get_pool(n_jobs: int) -> ProcessPoolExecutor:
    global _POOL, _POOL_N_JOBS
    if _POOL is None or _POOL_N_JOBS != n_jobs:
        if _POOL is not None:
            _POOL.shutdown(wait=True)
        spawn_context = multiprocessing.get_context("spawn")
        _POOL = ProcessPoolExecutor(
            max_workers=n_jobs, mp_context=spawn_context, initializer=_worker_init
        )
        _POOL_N_JOBS = n_jobs
    return _POOL


def generate_from_smiles(smiles: str, n: int = 100, n_jobs: int = -1) -> list:
    """Generates `n` new molecules sharing the input SMILES's exact molecular formula.

    `n_jobs` controls parallelism across the `n` independent generations: -1 uses
    a fixed pool of 4 workers (a persistent worker pool is created once and
    reused across calls) — each worker independently loads PyTorch plus both
    checkpoints, so even half the available cores (8 on a 16-core machine) was
    enough to exceed 2GB RSS per worker within seconds and get the whole
    process OOM-killed alongside typical desktop load; 1 runs sequentially
    in-process with no pool overhead. A positive integer requests that many
    workers exactly.

    Returns a list of `n` SMILES strings. Any generation attempt that fails after
    retries is filled with an empty string so a single bad draw never fails the
    whole call.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [""] * n
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    if n_jobs == 1:
        return [_generate_one_with_retries(formula) for _ in range(n)]

    resolved_n_jobs = 4 if n_jobs == -1 else n_jobs
    pool = _get_pool(resolved_n_jobs)
    return list(pool.map(_generate_one_with_retries, [formula] * n))
