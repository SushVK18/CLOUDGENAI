variable "vpc_id" {
  description = "The VPC ID in which the load balancer and security groups will be created."
  type        = string
}

variable "public_subnets_ids" {
  type        = list(string)
  description = "List of IDs for the public subnets where the ECS service will be deployed"
}

variable "project_prefix" {
  type        = string
  description = "Prefix used for naming resources in the project"
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

variable "session_timeout" {
  type        = string
  description = "The maximum duration of the authentication session, in seconds. Default value is 3600 seconds"
  default     = 3600
}

variable "target_group_settings" {
  description = "Map of target group configurations. Each key identifies a target group and its settings."
  type = any
  default = {
    app-tf-molmim-tg = {
      protocol              = "HTTP"
      port                  = 8000
      health_check_protocol = "HTTP"
      health_check_path     = "/v1/health/ready"
      protocol_version      = "HTTP1"
    },
    app-tf-diffdock-tg = {
      protocol              = "HTTP"
      port                  = 8000
      health_check_protocol = "HTTP"
      health_check_path     = "/v1/health/ready"
      protocol_version      = "HTTP1"
    },
    app-tf-genmol-tg = {
      protocol              = "HTTP"
      port                  = 8000
      health_check_protocol = "HTTP"
      health_check_path     = "/v1/health/ready"
      protocol_version      = "HTTP1"
    },
    app-tf-framework-tg = {
      protocol              = "HTTP"
      port                  = 8888
      health_check_protocol = "HTTP"
      health_check_path     = "/lab"
      protocol_version      = "HTTP1"
    }
  }
}

variable "tags" {
  type        = map(string)
  description = "Additional tags to be applied to the resources in this module."
  default     = {}
}