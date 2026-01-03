resource "aws_ecr_repository" "ecr_repos" {
  for_each = { for repo in var.repositories : repo => repo }
  
  name = each.value

  image_tag_mutability = "MUTABLE"
  image_scanning_configuration {
    scan_on_push = true
  }

  tags = var.tags
}
