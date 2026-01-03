resource "aws_autoscaling_group" "asg" {
  name                = var.asg_name
  vpc_zone_identifier = var.private_subnets

  desired_capacity = var.desired_capacity
  min_size         = var.min_size
  max_size         = var.max_size

  mixed_instances_policy {
    launch_template {
      launch_template_specification {
        launch_template_id = var.launch_template_id
        version            = "$Latest"
      }

      dynamic "override" {
        for_each = var.instance_types
        content {
          instance_type     = override.value
          weighted_capacity = "1"
        }
      }
    }
 
    instances_distribution {
      on_demand_percentage_above_base_capacity = var.on_demand_percentage_above_base_capacity
    }
  }

  protect_from_scale_in     = true

  tag {
    key                 = "Project"
    value               = var.project_name
    propagate_at_launch = true
  }

  tag {
    key                 = "Owner"
    value               = var.owner
    propagate_at_launch = true
  }

}
