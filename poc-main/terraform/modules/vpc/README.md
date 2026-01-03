# AWS VPC Module

Terraform module which creates VPC resources on AWS.

## Prerequisites

- Github SSH access configured
- AWS user/role with sufficient permission configured

## Usage

```hcl

terraform {
  required_providers {
    aws = {
      source = "hashicorp/aws"
      version = ">= 5.0, <= 5.2.0"
    }
  }
}


module "vpc-network" {
    source          = "git::git@github.com:bayer-int/efds-terraform-modules.git//vpc?ref=v1.1.0"

    vpc_name        = "<NAME-OF-THE-VPC>"
    ### for a full list of variables, have a look at variables.tf file ###
}
```

## Requirements

| Name                                                                       | Version   |
|----------------------------------------------------------------------------|-----------|
| <a name="requirement_terraform"></a> [terraform](#requirement\_terraform)  | >= 1.0 |
| <a name="requirement_aws"></a> [aws provider](#requirement\_aws)                    | >= 5.0, <= 5.2.0   |

## Providers

| Name                                              | Version |
|---------------------------------------------------|---------|
| <a name="provider_aws"></a> [aws](#provider\_aws) | >= 5.0, <= 5.2.0 |

## Inputs

| Name            | Description                                                      | Type     | Default     | Required |
|-----------------|------------------------------------------------------------------|----------|-------------|:--------:|
| vpc\_name       | provide a name for you vpc which will be created                 | `str`    |             |   yes    |
| cidr\_block     | CIDR block of the VPC   | `string` | "10.0.0.0/16" |   no, if not provided default value will be used    |
| az\_ids     | Available zone IDs i.e. us-east-1a   | list(`string`) |  ["us-east-1a", "us-east-1b"] |   no, if not provided default value will be used |
| private\_subnets    | A list of private subnets cidr block/s inside the VPC   | list(`string`) |  ["10.0.10.0/24", "10.0.20.0/24"] |   no, if not provided default value will be used |
| public\_subnets    | A list of public subnets cidr block/s inside the VPC   | list(`string`) |  ["10.0.101.0/24", "10.0.102.0/24"] |   no, if not provided default value will be used |

Refer to ``variables.tf`` to find out details about input variables.

## Outputs

| Name | Description |
|------|-------------|
| vpc_id | ID of the VPC |
| private_subnets | List of IDs of private subnets |
| public_subnets | List of IDs of public subnets |
| public_subnets_cidr_blocks | list of public subnet CIDR block |
| private_subnets_cidr_blocks | list of private subnet CIDR block |
| vpc_cidr_block | VPC CIDR block |
private_route_table_ids | Private route table IDs |