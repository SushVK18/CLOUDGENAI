# Generate a new RSA key pair (private key and corresponding public key).
resource "tls_private_key" "key" {
  algorithm = "RSA"
  rsa_bits  = 2048
}

# Create an AWS key pair using the generated public key.
resource "aws_key_pair" "generated_key" {
  key_name   = var.key_pair_name
  public_key = tls_private_key.key.public_key_openssh
}

# Compute the user data script to configure ECS on the instance.
locals {
  user_data = <<EOF
#!/bin/bash
echo ECS_CLUSTER=${var.ecs_cluster_name} >> /etc/ecs/ecs.config
echo ECS_ENABLE_GPU_SUPPORT=true >> /etc/ecs/ecs.config
EOF
}

# Create the EC2 Launch Template.
resource "aws_launch_template" "ec2_lt" {
  name = var.launch_template_name
  image_id = var.ami_id
  
  # Do not include instance type here; it will be defined by the Auto Scaling Group.
  key_name = aws_key_pair.generated_key.key_name

  iam_instance_profile {
    arn = var.iam_instance_profile_arn
  }

  # Block device mapping for the root EBS volume.
  block_device_mappings {
    device_name = var.ebs_device_name
    ebs {
      volume_size = var.ebs_volume_size
      volume_type = var.ebs_volume_type
    }
  }

  # You supply the existing security group IDs (e.g. from your load-balancer module).
  vpc_security_group_ids = var.security_group_ids

  # Tag the EC2 instance that will be created by this launch template.
  tag_specifications {
    resource_type = "instance"
    tags          = var.template_tags
  }

  # Tag the EBS volumes attached to the instance.
  tag_specifications {
    resource_type = "volume"
    tags          = var.resource_tags
  }

  # Provide your user data script after base64 encoding it.
  user_data = base64encode(local.user_data)
}
