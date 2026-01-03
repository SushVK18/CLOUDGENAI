data "aws_caller_identity" "current" {}
data "aws_region" "current" {}

##############################
# 1. ECS Cluster and Capacity Provider
##############################

resource "aws_ecs_cluster" "this" {
  name = var.cluster_name
}

resource "aws_ecs_capacity_provider" "this" {
  name = var.capacity_provider_name

  auto_scaling_group_provider {
    auto_scaling_group_arn = var.ecs_asg_arn
    managed_scaling {
      status                    = "ENABLED"
      target_capacity           = 80
      minimum_scaling_step_size = 1
      maximum_scaling_step_size = 1000
    }
    managed_termination_protection = "ENABLED"
  }
}

resource "aws_ecs_cluster_capacity_providers" "this" {
  cluster_name       = aws_ecs_cluster.this.name
  capacity_providers = [aws_ecs_capacity_provider.this.name]

  default_capacity_provider_strategy {
    capacity_provider = aws_ecs_capacity_provider.this.name
    weight            = 1
    base              = 1
  }
}

##############################
# 2. ECS Task Definitions
##############################

locals {
  # Define the common tasks.
  tasks = {
    "molmim-deploy" = {
      container_name         = "molmim"
      image                  = "${data.aws_caller_identity.current.account_id}.dkr.ecr.${data.aws_region.current.name}.amazonaws.com/molmim:1.0.0"
      container_port         = 8000
      host_port              = 8000
      cpu                    = 2048
      memory                 = 8192
      environment            = [
                                 { name = "CUDA_VISIBLE_DEVICES", value = "0" },
                                 { name = "LOCAL_NIM_CACHE",      value = "~/.cache/nim" },
                                 { name = "NGC_API_KEY",          value = var.ngc_api_key }
                               ]
      resource_requirements  = [ { type = "GPU", value = "1" } ]
      execution_role_arn     = var.ecs_task_execution_role_arn
      runtime_platform       = { cpu_architecture = "X86_64", operating_system_family = "LINUX" }
    },
    "diffdock-deploy" = {
      container_name         = "diffdock"
      image                  = "${data.aws_caller_identity.current.account_id}.dkr.ecr.${data.aws_region.current.name}.amazonaws.com/diffdock:latest"
      container_port         = 8000
      host_port              = 8000
      cpu                    = 4096
      memory                 = 12288
      environment            = [
                                 { name = "CUDA_VISIBLE_DEVICES", value = "0" },
                                 { name = "LOCAL_NIM_CACHE",      value = "~/.cache/nim" },
                                 { name = "NGC_API_KEY",          value = var.ngc_api_key }
                               ]
      resource_requirements  = [ { type = "GPU", value = "1" } ]
      execution_role_arn     = var.ecs_task_execution_role_arn
      runtime_platform       = { cpu_architecture = "X86_64", operating_system_family = "LINUX" }
    },
    "genmol-deploy" = {
      container_name         = "genmol"
      image                  = "${data.aws_caller_identity.current.account_id}.dkr.ecr.${data.aws_region.current.name}.amazonaws.com/genmol:latest"
      container_port         = 8000
      host_port              = 8000
      cpu                    = 2048
      memory                 = 8192
      environment            = [
                                 { name = "CUDA_VISIBLE_DEVICES", value = "0" },
                                 { name = "LOCAL_NIM_CACHE",      value = "~/.cache/nim" },
                                 { name = "NGC_API_KEY",          value = var.ngc_api_key }
                               ]
      resource_requirements  = [ { type = "GPU", value = "1" } ]
      execution_role_arn     = var.ecs_task_execution_role_arn
      runtime_platform       = { cpu_architecture = "X86_64", operating_system_family = "LINUX" }
    }
    # "framework-deploy" = {
    #   container_name         = "bionemo-framework"
    #   image                  = "${data.aws_caller_identity.current.account_id}.dkr.ecr.${data.aws_region.current.name}.amazonaws.com/bionemo-framework:1.10.1"
    #   container_port         = 8888
    #   host_port              = 8888
    #   cpu                    = 4096
    #   memory                 = 12288
    #   environment            = [
    #                              { name = "CUDA_VISIBLE_DEVICES", value = "0" },
    #                              { name = "LOCAL_NIM_CACHE",      value = "~/.cache/nim" },
    #                              { name = "NGC_API_KEY",          value = var.ngc_api_key }
    #                            ]
    #   resource_requirements  = [ { type = "GPU", value = "1" } ]
    #   execution_role_arn     = var.ecs_task_execution_role_arn
    #   runtime_platform       = { cpu_architecture = "X86_64", operating_system_family = "LINUX" }
    #   entry_point            = ["/bin/sh", "-c"]
    #   command                = [ "jupyter lab --allow-root --ip=* --port=8888 --no-browser --NotebookApp.token='' --NotebookApp.allow_origin='*' --ContentsManager.allow_hidden=True --notebook-dir=/workspace/bionemo" ]
    # }
  }
}

# Create a task definition for each task.
resource "aws_ecs_task_definition" "task" {
  for_each = local.tasks

  family                   = each.key
  network_mode             = "awsvpc"
  requires_compatibilities = ["EC2"]
  execution_role_arn       = each.value.execution_role_arn
  cpu                      = tostring(each.value.cpu)
  memory                   = tostring(each.value.memory)

  runtime_platform {
    cpu_architecture        = each.value.runtime_platform.cpu_architecture
    operating_system_family = each.value.runtime_platform.operating_system_family
  }

  container_definitions = jsonencode([
    merge(
      {
        name         = each.value.container_name,
        image        = each.value.image,
        essential    = true,
        portMappings = [
          {
            name          = "${each.value.container_name}-tf-port"
            containerPort = each.value.container_port,
            hostPort      = each.value.host_port,
            protocol      = "tcp",
            appProtocol   = "http"
          }
        ],
        environment          = each.value.environment,
        resourceRequirements = each.value.resource_requirements
        logConfiguration     = {
          logDriver = "awslogs"
          options = {
            awslogs-group         = "/ecs/${each.key}",
            awslogs-region        = "us-east-1"
            awslogs-stream-prefix = "ecs"
            awslogs-create-group  = "true"
            mode                  = "non-blocking"
            max-buffer-size       = "25m"
          }
        }
      },
      (lookup(each.value, "entry_point", null) != null ? { entryPoint = each.value.entry_point } : {} ),
      (lookup(each.value, "command", null) != null ? { command = each.value.command } : {} )
    )
  ])
}

##############################
# 3. ECS Services and Service Auto Scaling
##############################

# Create one ECS service per task.
resource "aws_ecs_service" "service" {
  for_each = local.tasks

  name            = "${each.key}-service"
  cluster         = aws_ecs_cluster.this.id
  task_definition = aws_ecs_task_definition.task[each.key].arn
  desired_count   = 1

  capacity_provider_strategy {
    capacity_provider = aws_ecs_capacity_provider.this.name
    weight            = 1
    base              = 1
  }

  network_configuration {
    subnets         = var.service_subnets
    security_groups = [var.service_security_group]
    assign_public_ip = false
  }

  load_balancer {
    target_group_arn = var.target_group_arns[each.key]
    container_name   = each.value.container_name
    container_port   = each.value.container_port
  }

  health_check_grace_period_seconds = 60

  deployment_controller {
    type = "ECS"
  }

  depends_on = [aws_ecs_cluster_capacity_providers.this]
}

# Configure service auto scaling targets and policies for each ECS service.
resource "aws_appautoscaling_target" "ecs_service" {
  for_each = local.tasks

  max_capacity       = 3
  min_capacity       = 0
  resource_id        = "service/${aws_ecs_cluster.this.name}/${aws_ecs_service.service[each.key].name}"
  scalable_dimension = "ecs:service:DesiredCount"
  service_namespace  = "ecs"
}

resource "aws_appautoscaling_policy" "ecs_service_policy" {
  for_each = local.tasks

  name               = "${each.key}-TargetTrackingPolicy"
  policy_type        = "TargetTrackingScaling"
  resource_id        = aws_appautoscaling_target.ecs_service[each.key].resource_id
  scalable_dimension = "ecs:service:DesiredCount"
  service_namespace  = "ecs"

  target_tracking_scaling_policy_configuration {
    target_value       = 85
    predefined_metric_specification {
      predefined_metric_type = "ECSServiceAverageCPUUtilization"
    }
    scale_in_cooldown  = 300
    scale_out_cooldown = 300
  }
}
