output "alb_listener" {
  description = "alb_listener"
  value       = aws_lb_listener.https.arn
}

output "alb_app_sg_id" {
  description = "The ID of the ALB application security group."
  value       = aws_security_group.alb_app_sg.id
}

output "service_app_sg_id" {
  description = "The ID of the service application security group."
  value       = aws_security_group.service_app_sg.id
}

output "molmim_target_group_arn" {
  description = "aws_lb_target_group app-tf-molmim-tg arn"
  value       = aws_lb_target_group.tg["app-tf-molmim-tg"].arn
}

output "diffdock_target_group_arn" {
  description = "aws_lb_target_group app-tf-diffdock-tg arn"
  value       = aws_lb_target_group.tg["app-tf-diffdock-tg"].arn
}

output "genmol_target_group_arn" {
  description = "aws_lb_target_group app-tf-genmol-tg arn"
  value       = aws_lb_target_group.tg["app-tf-genmol-tg"].arn
}

output "framework_target_group_arn" {
  description = "aws_lb_target_group app-tf-framework-tg arn"
  value       = aws_lb_target_group.tg["app-tf-framework-tg"].arn
}

output "molmim_endpoint_name" {
  description = "The connection endpoint for molmim"
  value       = aws_route53_record.molmim_subdomain.name
}

output "diffdock_endpoint_name" {
  description = "The connection endpoint for diffdock"
  value       = aws_route53_record.diffdock_subdomain.name
}

output "genmol_endpoint_name" {
  description = "The connection endpoint for genmol"
  value       = aws_route53_record.genmol_subdomain.name
}

output "framework_endpoint_name" {
  description = "The connection endpoint for bionemo framework"
  value       = aws_route53_record.jupyter_subdomain.name
}