output "repository_arns" {
  description = "A map of ECR repository ARNs keyed by repository name."
  value       = { for name, repo in aws_ecr_repository.ecr_repos : name => repo.arn }
}
