variable "secret_name" {
  description = "The name of the secret in Secrets Manager."
  type        = string
}

variable "secret_description" {
  description = "A description of the secret."
  type        = string
  default     = "Secret containing API keys and tokens"
}

variable "client_id" {
  description = "The API client ID."
  type        = string
  sensitive   = true
}

variable "client_secret" {
  description = "The API client secret."
  type        = string
  sensitive   = true
}

variable "ngc_api_key" {
  description = "The NGC API key."
  type        = string
  sensitive   = true
}

variable "tags" {
  description = "A map of tags to assign to the secret."
  type        = map(string)
  default     = {}
}
