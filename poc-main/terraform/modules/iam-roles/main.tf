data "aws_caller_identity" "current" {}


##############################
# 1. ECS ContainerService Role
##############################
resource "aws_iam_role" "ecs_role" {
  name = var.role_name
  assume_role_policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
      {
        "Effect": "Allow",
        "Principal": {
          "Service": "ec2.amazonaws.com"
        },
        "Action": "sts:AssumeRole"
      }
    ]
  })
  tags = var.tags
}

resource "aws_iam_instance_profile" "ecs_instance_profile" {
  name = var.role_name
  role = aws_iam_role.ecs_role.name
}

# Attach AWS Managed Policies

resource "aws_iam_role_policy_attachment" "ecr_full_access" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryFullAccess"
}

resource "aws_iam_role_policy_attachment" "ecs_autoscale_role" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceAutoscaleRole"
}

resource "aws_iam_role_policy_attachment" "ecs_for_ec2_role" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_role_policy_attachment" "ecs_role_role" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceRole"
}

resource "aws_iam_role_policy_attachment" "ec2_full_access" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2FullAccess"
}

resource "aws_iam_role_policy_attachment" "s3_full_access" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3FullAccess"
}

resource "aws_iam_role_policy_attachment" "secrets_manager_read_write" {
  role       = aws_iam_role.ecs_role.name
  policy_arn = "arn:aws:iam::aws:policy/SecretsManagerReadWrite"
}

# Inline Policy: Cloudwatch Reduced Policy

resource "aws_iam_role_policy" "cloudwatch_reduced_policy" {
  name = "Cloudwatch-reduced-policy"
  role = aws_iam_role.ecs_role.id
  
  policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
      {
        "Action": [
          "logs:CreateLogGroup",
          "logs:CreateLogStream",
          "logs:PutLogEvents",
          "logs:DescribeLogStreams"
        ],
        "Resource": [
          "arn:aws:logs:*:${data.aws_caller_identity.current.account_id}:log-group:*"
        ],
        "Effect": "Allow"
      }
    ]
  })
}

# Inline Policy: SSM Reduced Policy

resource "aws_iam_role_policy" "ssm_reduced_policy" {
  name = "ssm-reduced-policy"
  role = aws_iam_role.ecs_role.id
  
  policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
      {
        "Action": [
          "ssm:DescribeAssociation",
          "ssm:GetDeployablePatchSnapshotForInstance",
          "ssm:GetDocument",
          "ssm:GetParameters",
          "ssm:ListAssociations",
          "ssm:ListInstanceAssociations",
          "ssm:PutInventory",
          "ssm:PutComplianceItems",
          "ssm:PutConfigurePackageResult",
          "ssm:UpdateAssociationStatus",
          "ssm:UpdateInstanceAssociationStatus",
          "ssm:UpdateInstanceInformation",
          "ssmmessages:CreateControlChannel",
          "ssmmessages:CreateDataChannel",
          "ssmmessages:OpenControlChannel",
          "ssmmessages:OpenDataChannel"
        ],
        "Resource": "*",
        "Effect": "Allow"
      },
      {
        "Action": [
          "ec2messages:AcknowledgeMessage",
          "ec2messages:DeleteMessage",
          "ec2messages:FailMessage",
          "ec2messages:GetEndpoint",
          "ec2messages:GetMessages",
          "ec2messages:SendReply"
        ],
        "Resource": "*",
        "Effect": "Allow"
      },
      {
        "Action": [
          "ec2:DescribeInstanceStatus"
        ],
        "Resource": "*",
        "Effect": "Allow"
      },
      {
        "Action": [
          "s3:PutObject"
        ],
        "Resource": "arn:aws:s3:::${data.aws_caller_identity.current.account_id}-ssm-patch-logs/*",
        "Effect": "Allow"
      }
    ]
  })
}


##############################
# 2. ECS Task Service Role
##############################
resource "aws_iam_role" "ecs_task_role" {
  name = "${var.role_name}-task"
  assume_role_policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
      {
        "Effect": "Allow",
        "Principal": {
          "Service": "ecs-tasks.amazonaws.com"
        },
        "Action": "sts:AssumeRole"
      }
    ]
  })
  tags = var.tags
}

resource "aws_iam_instance_profile" "ecs_task_instance_profile" {
  name = "${var.role_name}-task"
  role = aws_iam_role.ecs_task_role.name
}

# Attach AWS Managed Policies
resource "aws_iam_role_policy_attachment" "task_ecr_full_access" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryFullAccess"
}

resource "aws_iam_role_policy_attachment" "task_ecs_autoscale_role" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceAutoscaleRole"
}

resource "aws_iam_role_policy_attachment" "task_ecs_for_ec2_role" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_role_policy_attachment" "task_ecs_role_role" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceRole"
}

resource "aws_iam_role_policy_attachment" "task_ecs_full_access" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonECS_FullAccess"
}

resource "aws_iam_role_policy_attachment" "task_ecs_task_execution_role" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

resource "aws_iam_role_policy_attachment" "task_secrets_manager_read_write" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = "arn:aws:iam::aws:policy/SecretsManagerReadWrite"
}


##############################
# 3. Lambda Execution Role
##############################
resource "aws_iam_role" "lambda_role" {
  name = "lambda-execution-tf-role"
  assume_role_policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
      {
        "Effect": "Allow",
        "Principal": {
            "Service": "lambda.amazonaws.com"
        },
        "Action": "sts:AssumeRole"
      },
      {
        "Effect": "Allow",
        "Principal": {
            "AWS": "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/cloudops"
        },
        "Action": "sts:AssumeRole"
      }
    ]
  })
  tags = var.tags
}

resource "aws_iam_instance_profile" "lambda_role_instance_profile" {
  name = aws_iam_role.lambda_role.name
  role = aws_iam_role.lambda_role.name
}

# Attach AWS Managed Policies
resource "aws_iam_role_policy_attachment" "lambda_basic_exec_role" {
  role       = aws_iam_role.lambda_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole"
}