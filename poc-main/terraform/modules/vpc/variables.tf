variable "vpc_name" {
  type        = string
  description = "name of the vpc to be created"
}
variable "cidr_block" {
  type        = string
  description = "list of cidr block for the vpc to be created"
  default     = "10.0.0.0/16"
}
variable "az_ids" {
  type        = list(string)
  description = "list of availability zones vpc subnets will be spreaded into"
  default     = ["us-east-1a", "us-east-1b"]
}
variable "private_subnets" {
  type        = list(string)
  description = "A list of private subnets cidr block/s inside the VPC"
  default     = ["10.0.10.0/24", "10.0.20.0/24"]
}
variable "public_subnets" {
  type        = list(string)
  description = "A list of public subnets cidr block/s inside the VPC"
  default     = ["10.0.101.0/24", "10.0.102.0/24"]
}

variable "tags" {
  type        = map(string)
  description = "Additional tags to be applied to the resources in this module."
  default     = {}
}