# CoCoGraph Formula-Preserving Generation

Generates 100 new, chemically valid molecules that share the exact atomic composition of an input compound, using a constrained graph diffusion process trained on 2.25 million molecules from PubChem, ChEMBL, ZINC and NIST. Validity is guaranteed by design rather than learned, since every diffusion step swaps bonds while preserving atom valences. In a blind test, 121 organic chemistry experts correctly distinguished generated from real molecules only 62% of the time. Training and inference are restricted to molecules of 5-70 total atoms, hydrogens included.



## Information
### Identifiers
- **Ersilia Identifier:** `eos55vx`
- **Slug:** `cocograph-formula`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Compound generation`, `Chemical graph model`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Variable`
- **Interpretation:** 100 newly generated molecules sharing the input compound's exact molecular formula and valid chemical structure.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Generated molecule 0 sharing the input compound's exact molecular formula |
| smi_01 | string |  | Generated molecule 1 sharing the input compound's exact molecular formula |
| smi_02 | string |  | Generated molecule 2 sharing the input compound's exact molecular formula |
| smi_03 | string |  | Generated molecule 3 sharing the input compound's exact molecular formula |
| smi_04 | string |  | Generated molecule 4 sharing the input compound's exact molecular formula |
| smi_05 | string |  | Generated molecule 5 sharing the input compound's exact molecular formula |
| smi_06 | string |  | Generated molecule 6 sharing the input compound's exact molecular formula |
| smi_07 | string |  | Generated molecule 7 sharing the input compound's exact molecular formula |
| smi_08 | string |  | Generated molecule 8 sharing the input compound's exact molecular formula |
| smi_09 | string |  | Generated molecule 9 sharing the input compound's exact molecular formula |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`

### Resource Consumption


### References
- **Source Code**: [https://doi.org/10.5281/zenodo.18940151](https://doi.org/10.5281/zenodo.18940151)
- **Publication**: [https://doi.org/10.1038/s42256-026-01229-5](https://doi.org/10.1038/s42256-026-01229-5)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2026`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos55vx
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos55vx
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
