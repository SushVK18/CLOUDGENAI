variable "role_name" {
  description = "The name for the ECS IAM role."
  type        = string
}

variable "tags" {
  description = "A map of tags to assign to the IAM role."
  type        = map(string)
  default     = {}
}
