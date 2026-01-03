variable "aws_region" {
  description = "The AWS region where resources will be deployed"
  type        = string
  default     = "us-east-1"
}

variable "environment" {
  description = "The name of the environment (e.g., dev, staging, prod)"
  type        = string
  default     = "dev"
}

variable "certificate_arn" {
  type        = string
  description = "The unique identifier (ARN) of the AWS ACM SSL/TLS certificate"
  sensitive   = true
}

variable "route53_zone_id" {
  description = "Route53 Hosted Zone ID to create DNS records for both certificate validation and ALB alias."
  type        = string
  sensitive   = true
}

variable "client_id" {
  type        = string
  description = "The client ID for the OIDC registered client."
  sensitive   = true
}

variable "client_secret" {
  type        = string
  description = "The secret associated with the OIDC client ID."
  sensitive   = true
}

variable "ngc_api_key" {
  type        = string
  description = "The NGC API key."
  sensitive   = true
}

variable "project_name" {
  description = "The name of the project for tagging purposes."
  type        = string
  default     = "smol-bionemo-poc"
}

variable "owner" {
  description = "The name of the owner for tagging purposes."
  type        = string
  default     = "CG-ext"
}