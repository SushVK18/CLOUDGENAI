provider "aws" {
  region = var.aws_region

  default_tags {
    tags = local.common_tags
  }
}

module "vpc" {
  source = "../../modules/vpc"

  vpc_name = "bionemo_tf_vpc"
}

module "load-balancer" {
  source = "../../modules/load-balancer"

  vpc_id = module.vpc.vpc_id
  public_subnets_ids  = module.vpc.public_subnets
  project_prefix = "bionemo-tf-poc"
  certificate_arn = var.certificate_arn
  route53_zone_id = var.route53_zone_id
  client_id = var.client_id
  client_secret = var.client_secret
}

module "s3-bucket" {
  source = "../../modules/s3-bucket"

  bucket_name = "bionemo-tf-s3"
}

module "secrets_manager" {
  source = "../../modules/secrets-manager"

  secret_name = "bionemo-tf-secrets"
  client_id = var.client_id
  client_secret = var.client_secret
  ngc_api_key = var.ngc_api_key
}

module "ecr_repos" {
  source = "../../modules/ecr"
  
  repositories = ["bionemo-framework-tf", "diffdock-tf", "genmol-tf", "molmim-tf"]
}

module "ecs_iam_roles" {
  source    = "../../modules/iam-roles"

  role_name = "Bionemo-EC2ContainerService-tf"
}

module "ec2_launch_template" {
  source = "../../modules/launch-template"

  launch_template_name  = "bionemo-tf-lt"
  ami_id = "ami-08096c2d3b6642595" # Bayer ECS optimized EC2 with GPU support
  key_pair_name = "bionemo-tf-keypair"
  ebs_device_name = "/dev/xvda"
  ebs_volume_size = 512
  ebs_volume_type = "gp3"
  iam_instance_profile_arn = module.ecs_iam_roles.instance_profile_arn
  
  security_group_ids = [module.load-balancer.alb_app_sg_id]
  
  template_tags = local.common_tags
  resource_tags = local.common_tags
  
  ecs_cluster_name = "bionemo-tf-inference"
}

module "ec2_asg" {
  source = "../../modules/ec2-asg"
  
  asg_name = "bionemo-tf-autoscaling"
  launch_template_id = module.ec2_launch_template.launch_template_id 
  instance_types = ["g5.xlarge"]  # You could add more types if desired
  private_subnets = module.vpc.private_subnets         # from your VPC module
  desired_capacity = 4
  min_size = 0
  max_size = 10

}

module "ecs_ec2_deploy" {
  source = "../../modules/ecs-ec2-deploy"

  cluster_name = "bionemo-tf-inference"
  launch_template_id = module.ec2_launch_template.launch_template_id
  ecs_asg_arn = module.ec2_asg.asg_arn
  capacity_provider_name = "bionemo-capacity-provider"
  ecs_task_execution_role_arn = module.ecs_iam_roles.task_execution_role_arn
  ngc_api_key = var.ngc_api_key
  service_subnets = module.vpc.private_subnets
  service_security_group = module.load-balancer.service_app_sg_id
  target_group_arns = {
    "molmim-deploy"   = module.load-balancer.molmim_target_group_arn,
    "diffdock-deploy" = module.load-balancer.diffdock_target_group_arn,
    "genmol-deploy"   = module.load-balancer.genmol_target_group_arn,
    "framework-deploy"= module.load-balancer.framework_target_group_arn
  }
}

module "api_gateway_lambda" {
  source = "../../modules/api-gateway-lambda"

  lambda_function_name = "azure-api-gateway-tf-authorizer"
  lambda_exection_role_arn = module.ecs_iam_roles.lambda_execution_role_arn
  lambda_zip_path = "../../modules/api-gateway-lambda/azure-apigw-authorizer-lambda.zip"
  audience = var.client_id

  api_name = "bionemo-poc-tf-gateway"

}
