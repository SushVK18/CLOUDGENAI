output "launch_template_id" {
  description = "The ID of the created EC2 launch template."
  value       = aws_launch_template.ec2_lt.id
}

output "launch_template_name" {
  description = "The name of the created EC2 launch template."
  value       = aws_launch_template.ec2_lt.name
}

output "key_pair_name" {
  description = "The name of the generated key pair."
  value       = aws_key_pair.generated_key.key_name
}

output "private_key" {
  description = "The generated private key material. Store this securely."
  value       = tls_private_key.key.private_key_pem
  sensitive   = true
}
