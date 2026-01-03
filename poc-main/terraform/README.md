# BioNemo Deployment Blueprint
This repository contains the Terraform code necessary for provisioning cloud infrastructure on AWS. It supports BioNemo framework and NIMs deployment by creating a robust, scalable, and secure environment, leveraging AWS ECS with GPU-supported EC2 instance, an Application Load Balancer for efficient traffic distribution, and a VPC with both public and private subnets for security and resource optimization.


## Repository Structure
```
├── terraform                               <- Root module which combines submodules for a complete infrastructure setup
│   ├── envs                                <- Environment
│   │   └── dev                             <- Development environment
│   │       ├── locals.tf                   <- Configure special variables that are evaluated during execution of the plan phase
│   │       ├── main.tf                     <- Main Terraform file to tie modules and resources together
│   │       ├── outputs.tf                  <- Defines outputs from the Terraform modules
│   │       ├── secrets.tfvars.template     <- Template for sensitive variables that are evaluated during execution of the plan phase
│   │       └── variables.tf                <- Contains variable definitions that are populated during the execution
│   ├── modules                             <- Contains individual Terraform modules
│   │   ├── api-gateway-lambda              <- Module to provide (HTTP) API Gateway and Lambda Authorizer
│   │   ├── ec2-asg                         <- Module to set up the EC2 Auto Scaling Group
│   │   ├── ecr                             <- Module to set up the private ECR repository
│   │   ├── ecs-ec2-deploy                  <- Module to set up the ECS cluster and EC2 services
│   │   ├── iam-roles                       <- Module for provisioning all necessary roles in the deployment blueprint
│   │   ├── launch-template                 <- Module to configure EC2 launch template
│   │   ├── load-balancer                   <- Module for provisioning an Application Load Balancer
│   │   ├── s3-bucket                       <- Module for provisioning a S3 bucket
│   │   ├── secrets-manager                 <- Module for provisioning a secret manager
│   │   └── vpc                             <- Module to create a VPC with public and private subnets across multiple availability zones
└── .gitignore                              <- Excludes specific files and folders for Git
--------------------------------------------------------------------
```


## Prerequisites
1. AWS user/role with sufficient permission configured
2. Configure AWS ACM - SSL certificate.
    - https://confluence.cloud.bayer.com/spaces/DECSCI/pages/541557785/AWS+ACM+-+Requesting+SSL+Certificate
3. Generate NVIDIA NGC_API_KEY and subscribe to developer program.
    - https://org.ngc.nvidia.com/setup/api-key
    - https://developer.nvidia.com/developer-program 
4. Register/integrate the application with Bayer Entra ID
    - https://docs.int.bayer.com/cloud/general/cloud_iam/applications/ 


## Project related variables
Create a `secrets.tfvars` file based on the template `secrets.tfvars.template`. Edit the values on this file to match your project-related variable values.
- **certificate_arn**: ARN of of the AWS ACM SSL/TLS certificate.
- **route53_zone_id**: Route53 Hosted Zone ID.
- **client_id**: The client ID for the OIDC registered client.
- **client_secret**: The secret associated with the OIDC client ID.
- **ngc_api_key**: NVIDIA NGC API key.


## Azure AD in the ALB for authentication
In order to perform the authentication via Azure AD in the ALB, you must first register your Azure AD application with a BEAT ID. Then you need to configure the Azure OIDC/OAuth2 for Web-UIs and APIs. For guidance on registering the application and configuring the Azure AD see:
- https://docs.int.bayer.com/cloud/general/cloud_iam/applications/
- https://docs.int.bayer.com/cloud/general/cloud_iam/howto_auth_auth/
- https://docs.int.bayer.com/cloud/general/cloud_iam/howto_app_configuration_api/
- https://docs.int.bayer.com/cloud/general/cloud_iam/entraid_managed_apis/


## Terraform Execution
- **Preparation:** Change directory to the Terraform working directory. Example: `cd ./terraform/envs/dev`. 
- **Initialization:** Run `terraform init` to initialize the Terraform working directory.
- **Planning:** Run `terraform plan -var-file="secrets.tfvars"` to create an execution plan.
- **Applying:** Run `terraform apply -var-file="secrets.tfvars"` to apply the changes required to reach the desired state of the configuration.
- **Destruction (Caution):** Run `terraform destroy -var-file="secrets.tfvars"` to destroy the Terraform-managed infrastructure.