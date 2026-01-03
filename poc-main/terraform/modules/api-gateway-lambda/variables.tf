variable "api_name" {
  description = "The name of the API Gateway."
  type        = string
}

variable "lambda_function_name" {
  description = "The name of the lambda function."
  type        = string
  default     = "azure-api-gateway-authorizer"
}

variable "lambda_exection_role_arn" {
  description = "The ARN of the lambda execution IAM role."
  type        = string
}

variable "lambda_zip_path" {
  description = "The path to the zipped lambda function."
  type        = string
}

variable "algorithms" {
  description = "Token signing algorithm."
  type        = string
  default     = "RS256"
}

variable "audience" {
  description = "oken audience, set this for extra security."
  type        = string
  default     = ""
}

variable "issuersurls" {
  description = "Base URL(s) of the Azure AD tenant."
  type        = string
  default     = "https://sts.windows.net/fcb2b37b-5da0-466b-9b83-0014b67a7c78/"
}

variable "tags" {
  description = "A map of tags to assign to the resource."
  type        = map(string)
  default     = {}
}
