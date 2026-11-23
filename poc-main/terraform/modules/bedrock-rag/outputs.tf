output "knowledge_base_id" {
  description = "ID of the Bedrock knowledge base."
  value       = aws_bedrockagent_knowledge_base.this.id
}

output "knowledge_base_arn" {
  description = "ARN of the Bedrock knowledge base."
  value       = aws_bedrockagent_knowledge_base.this.arn
}

output "data_source_id" {
  description = "ID of the S3-backed Bedrock data source."
  value       = aws_bedrockagent_data_source.documents.data_source_id
}

output "documents_bucket_name" {
  description = "Name of the private S3 bucket that stores source documents."
  value       = aws_s3_bucket.documents.bucket
}

output "application_role_arn" {
  description = "ECS task role with scoped Bedrock knowledge base retrieval permissions."
  value       = aws_iam_role.application.arn
  depends_on  = [aws_iam_role_policy.application]
}

output "generation_model_id" {
  description = "Bedrock foundation model ID configured for RAG response generation."
  value       = var.generation_model_id
}
