output "ecs_cluster_id" {
  description = "The ID of the ECS cluster."
  value       = aws_ecs_cluster.this.id
}

output "ecs_cluster_name" {
  description = "The name of the ECS cluster."
  value       = aws_ecs_cluster.this.name
}

output "ecs_service_names" {
  description = "Map of ECS service names keyed by task definition."
  value       = { for k, svc in aws_ecs_service.service : k => svc.name }
}

output "ecs_capacity_provider_name" {
  description = "The name of the ECS capacity provider."
  value       = aws_ecs_capacity_provider.this.name
}
