variable "cluster_name" {
  description = "ECS cluster name."
  type        = string
  default     = "bionemo-cluster"
}

variable "launch_template_id" {
  description = "ID of the EC2 launch template used for provisioning container instances."
  type        = string
}

variable "instance_type" {
  description = "Instance type to use for the Auto Scaling Group via mixed instances policy."
  type        = string
  default     = "g5.2xlarge"
}

variable "ecs_asg_arn"{
  description = "Auto Scaling Group ARN."
  type        = string
}

variable "capacity_provider_name" {
  description = "Name for the ECS capacity provider."
  type        = string
}

variable "ecs_task_execution_role_arn" {
  description = "IAM role ARN for ECS task execution (e.g. ecsTaskExecutionRole)."
  type        = string
}

variable "task_role_arn" {
  description = "The IAM role ARN used by containers to call AWS services."
  type        = string
  default     = null
}

variable "bedrock_knowledge_base_id" {
  description = "Optional Bedrock knowledge base ID exposed to task containers."
  type        = string
  default     = null
}

variable "bedrock_generation_model_id" {
  description = "Optional Bedrock foundation model ID exposed to task containers."
  type        = string
  default     = null
}

variable "bedrock_region" {
  description = "Optional AWS region for Bedrock requests made by task containers."
  type        = string
  default     = null
}

variable "ngc_api_key" {
  description = "NGC API Key."
  type        = string
}

variable "service_subnets" {
  description = "List of subnet IDs for ECS service tasks (typically private subnets)."
  type        = list(string)
}

variable "service_security_group" {
  description = "Security group ID for ECS service tasks (e.g., service-app-sg)."
  type        = string
}

variable "target_group_arns" {
  description = "Map of target group ARNs keyed by task key. For example: { \"molmim-deploy\" = \"<ARN>\", \"diffdock-deploy\" = \"<ARN>\", \"genmol-deploy\" = \"<ARN>\", \"framework-deploy\" = \"<ARN>\" }"
  type        = map(string)
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
