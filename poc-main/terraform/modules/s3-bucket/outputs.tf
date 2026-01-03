output "s3_bucket_arn" {
  description = "Arn of the created s3 bucket"
  value       = module.s3_bucket.s3_bucket_arn
}