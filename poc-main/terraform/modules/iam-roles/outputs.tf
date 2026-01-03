output "role_arn" {
  description = "The ARN of the ECS IAM role."
  value       = aws_iam_role.ecs_role.arn
}

output "role_name" {
  description = "The name of the ECS IAM role."
  value       = aws_iam_role.ecs_role.name
}

output "instance_profile_arn" {
  description = "The ARN of the instance profile associated with the ECS IAM role."
  value       = aws_iam_instance_profile.ecs_instance_profile.arn
}

output "task_execution_role_arn" {
  description = "The ARN of the ECS task execution IAM role."
  value       = aws_iam_role.ecs_task_role.arn
}

output "task_execution_role_name" {
  description = "The name of the ECS task execution IAM role."
  value       = aws_iam_role.ecs_task_role.name
}

output "task_execution_instance_profile_arn" {
  description = "The ARN of the instance profile associated with the ECS task execution IAM role."
  value       = aws_iam_instance_profile.ecs_task_instance_profile.arn
}

output "lambda_execution_role_arn" {
  description = "The ARN of the lambda execution IAM role."
  value       = aws_iam_role.lambda_role.arn
}

output "lambda_execution_role_name" {
  description = "The name of the lambda execution IAM role."
  value       = aws_iam_role.lambda_role.name
}

output "lambda_execution_instance_profile_arn" {
  description = "The ARN of the instance profile associated with the lambda execution IAM role."
  value       = aws_iam_instance_profile.lambda_role_instance_profile.arn
}