variable "launch_template_name" {
  description = "Name for the launch template."
  type        = string
}

variable "ami_id" {
  description = "AMI ID to use for the launch template."
  type        = string
}

variable "key_pair_name" {
  description = "The name to assign to the generated key pair."
  type        = string
}

variable "ebs_device_name" {
  description = "Name of the device for the root EBS volume."
  type        = string
  default     = "/dev/xvda"
}

variable "ebs_volume_size" {
  description = "EBS volume size in GB."
  type        = number
  default     = 512
}

variable "ebs_volume_type" {
  description = "EBS volume type."
  type        = string
  default     = "gp3"
}

variable "iam_instance_profile_arn" {
  description = "The IAM instance profile arn to attach."
  type        = string
}

variable "security_group_ids" {
  description = "List of security group IDs to assign to the instances."
  type        = list(string)
}

variable "template_tags" {
  description = "A map of tags applied to the EC2 instance created by the launch template."
  type        = map(string)
}

variable "resource_tags" {
  description = "A map of tags applied to the EBS volumes attached to the instance."
  type        = map(string)
}

variable "ecs_cluster_name" {
  description = "Name of the ECS cluster where GPU needs to be enabled"
  type        = string
}
