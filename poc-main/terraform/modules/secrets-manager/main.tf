resource "aws_secretsmanager_secret" "sm-secret" {
  name        = var.secret_name
  description = var.secret_description
  tags        = var.tags
}

resource "aws_secretsmanager_secret_version" "sm-secret-version" {
  secret_id     = aws_secretsmanager_secret.sm-secret.id
  secret_string = jsonencode({
    CLIENT_ID     = var.client_id,
    CLIENT_SECRET = var.client_secret,
    NGC_API_KEY   = var.ngc_api_key
  })
}
