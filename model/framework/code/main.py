# imports
import sys

from ersilia_pack_utils.core import read_smiles, write_out

from sample_formula import generate_from_smiles

N_OUTPUTS = 100


def my_model(smiles_list):
    return [generate_from_smiles(smi, n=N_OUTPUTS, n_jobs=-1) for smi in smiles_list]


if __name__ == "__main__":
    # Required for spawn-based multiprocessing in sample_formula.py: this guard
    # stops each worker process from re-running the script body when it
    # re-imports __main__ on startup.

    # parse arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # read SMILES from .csv file, assuming one column with header
    _, smiles_list = read_smiles(input_file)

    # run model
    outputs = my_model(smiles_list)

    # check input and output have the same length
    assert len(smiles_list) == len(outputs)

    header = [f"smi_{str(i).zfill(2)}" for i in range(N_OUTPUTS)]

    # write output in a .csv file
    write_out(outputs, header, output_file)
