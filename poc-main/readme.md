# GenAI Models Poc

This repository contains content from the GenAI PoC . 

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

## CI/CD workflow

This repository includes a GitHub Actions pipeline for Terraform validation and deployment.

### Workflow behavior
- On every push or pull request to `main`, the workflow runs Terraform validation.
- It checks Terraform formatting with `terraform fmt -check -recursive`.
- It initializes and validates the dev environment under `poc-main/terraform/envs/dev`.
- The `workflow_dispatch` trigger allows a manual deployment run from the GitHub Actions UI.

### Required secrets
Set these repository secrets in GitHub before running the manual deployment:
- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `AWS_SESSION_TOKEN` (if your AWS credentials are temporary)
- `TF_VAR_certificate_arn`
- `TF_VAR_route53_zone_id`
- `TF_VAR_client_id`
- `TF_VAR_client_secret`
- `TF_VAR_ngc_api_key`

### How to trigger deployment
1. Open the repository in GitHub.
2. Go to `Actions`.
3. Select the `Terraform CI/CD` workflow.
4. Click `Run workflow` and choose the branch to deploy.

This pipeline is designed to validate infrastructure changes automatically and apply the Terraform configuration only when you trigger a deployment manually.
