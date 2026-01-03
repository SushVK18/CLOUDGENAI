variable "asg_name" {
  description = "The name of the Auto Scaling Group."
  type        = string
  default     = "bionemo-asg"
}

variable "launch_template_id" {
  description = "The ID of the launch template to use for the ASG."
  type        = string
}

variable "instance_types" {
  description = "List of allowed instance types. The primary instance type should be 'g5.2xlarge'."
  type        = list(string)
  default     = ["g5.2xlarge"]
}

variable "private_subnets" {
  description = "List of private subnet IDs (from the VPC module) for the ASG."
  type        = list(string)
}

variable "desired_capacity" {
  description = "The desired number of instances in the ASG."
  type        = number
  default     = 1
}

variable "min_size" {
  description = "The minimum number of instances in the ASG."
  type        = number
  default     = 1
}

variable "max_size" {
  description = "The maximum number of instances in the ASG."
  type        = number
  default     = 3
}

variable "on_demand_percentage_above_base_capacity" {
  description = "Percentage of on-demand instances above the base capacity in the mixed instances policy (100 means all on-demand)."
  type        = number
  default     = 100
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
