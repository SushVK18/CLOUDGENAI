variable "name_prefix" {
  description = "Lowercase project prefix used for Bedrock, OpenSearch, and S3 resource names."
  type        = string
  default     = "bionemo-tf-rag"

  validation {
    condition     = can(regex("^[a-z][a-z0-9-]{1,19}$", var.name_prefix))
    error_message = "name_prefix must start with a lowercase letter and contain only lowercase letters, digits, or hyphens (2-20 characters)."
  }
}

variable "embedding_model_id" {
  description = "Bedrock foundation model ID used to create document embeddings."
  type        = string
  default     = "amazon.titan-embed-text-v2:0"
}

variable "generation_model_id" {
  description = "Bedrock foundation model ID that application tasks may invoke for RAG responses."
  type        = string
  default     = "amazon.nova-lite-v1:0"
}

variable "tags" {
  description = "Tags to assign to supported resources."
  type        = map(string)
  default     = {}
}
