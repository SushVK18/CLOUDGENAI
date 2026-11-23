output "bedrock_knowledge_base_id" {
  description = "ID of the Bedrock RAG knowledge base."
  value       = module.bedrock_rag.knowledge_base_id
}

output "bedrock_data_source_id" {
  description = "ID of the S3 document data source configured for the knowledge base."
  value       = module.bedrock_rag.data_source_id
}

output "bedrock_documents_bucket" {
  description = "S3 bucket where source documents should be uploaded under documents/."
  value       = module.bedrock_rag.documents_bucket_name
}