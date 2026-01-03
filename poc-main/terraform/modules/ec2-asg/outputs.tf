output "asg_name" {
  description = "The name of the Auto Scaling Group."
  value       = aws_autoscaling_group.asg.name
}

output "asg_id" {
  description = "The ID of the Auto Scaling Group."
  value       = aws_autoscaling_group.asg.id
}

output "asg_arn" {
  description = "The ARN of the Auto Scaling Group."
  value       = aws_autoscaling_group.asg.arn
}
