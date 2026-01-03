variable "repositories" {
  description = "List of ECR repository names to create."
  type        = list(string)
  default     = [""]
}

variable "tags" {
  description = "A map of tags to assign to the ECR repositories."
  type        = map(string)
  default     = {}
}
