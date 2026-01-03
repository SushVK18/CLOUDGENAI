# BioNeMo Poc

This repository contains content from the CapgeminiXBayer PoC for BioNeMo. 

## Structure

The Project is seperated into the following main parts:

```
project/
├── notebooks-eda/
├── notebooks-eval/
│   ├── AF2/
│   ├── diffdock/
│   ├── genmol/
│   ├── molmim/
│   ├── Z00_Midpoint_Hurdle/
│   └── Z01_Final_Eval/
├── terraform/
└── README.md
```

## Folder Descriptions

- **notebooks-eda/**: Contains Notebooks used for an initial EDA as well as for data download (from S3 and from Uniprot)

- **notebooks-eval/**: Contains all notebooks and helper scripts to Start Nims, Inference and initial outcome analysis
    - **AF2/**: Alphafold2 inference scripts and notebooks - for structure prediction
    - **diffdock/**: DiffDock inference scripts and notebooks - for protein-ligand docking
    - **genmol/**: GenMoL inference scripts and notebooks - for Molecule generation and optimization
    - **molmim/**: MolMiM inference scripts and notebooks - for Molecule generation and optimization
    - **Z00_Midpoint_Hurdle/**: Notebooks for initial analysis and result visualization
    - **Z01_Final_Eval//**: Notebooks for final analysis and result visualization

- **terraform/**: Contains all terraform scripts fo a deployment via EC2

- **README.md**: Documentation for the project, including setup instructions and usage.


## Deployment 

Refer to `terraform/README.md` for Terraform deployment instructions.
