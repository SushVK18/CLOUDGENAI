data "aws_caller_identity" "current" {}
data "aws_region" "current" {}

resource "aws_cloudwatch_log_group" "lambda_logs" {
  name              = "/aws/lambda/${var.lambda_function_name}"
  retention_in_days = 7
}

resource "aws_lambda_function" "authorizer" {
  function_name = var.lambda_function_name
  runtime       = "nodejs20.x"
  role          = var.lambda_exection_role_arn
  handler       = "index.handler"
  filename      = var.lambda_zip_path
  source_code_hash = filebase64sha256(var.lambda_zip_path)

  environment {
    variables = {
      algorithms  = var.algorithms
      audience    = var.audience
      issuersurls = var.issuersurls
    }
  }
}

resource "aws_lambda_permission" "api_invoke" {
  statement_id  = "AllowAPIGatewayInvoke"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.authorizer.function_name
  principal     = "apigateway.amazonaws.com"
}

resource "aws_apigatewayv2_api" "http_api" {
  name          = var.api_name
  protocol_type = "HTTP"
}

resource "aws_apigatewayv2_authorizer" "lambda_auth" {
  name                       = "AzureTokenAuthorizer"
  api_id                     = aws_apigatewayv2_api.http_api.id
  authorizer_type            = "REQUEST"
  # authorizer_uri             = "arn:aws:apigateway:${data.aws_region.current.name}:lambda:path/2015-03-31/functions/${aws_lambda_function.authorizer.arn}/invocations"
  authorizer_uri             = aws_lambda_function.authorizer.invoke_arn
  identity_sources           = ["$request.header.Authorization"]
  authorizer_payload_format_version = "2.0"
  enable_simple_responses    = true
}

locals {
  routes = [
    {
      name        = "diffdock"
      path        = "/diffdock/{proxy+}"
      target_url  = "https://diffdock-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com/{proxy}"
    },
    {
      name        = "genmol"
      path        = "/genmol/{proxy+}"
      target_url  = "https://genmol-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com/{proxy}"
    },
    {
      name        = "molmim"
      path        = "/molmim/{proxy+}"
      target_url  = "https://molmim-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com/{proxy}"
    }
  ]
}

resource "aws_apigatewayv2_integration" "routes" {
  for_each             = { for route in local.routes : route.name => route }
  api_id               = aws_apigatewayv2_api.http_api.id
  integration_type     = "HTTP_PROXY"
  integration_method   = "ANY"
  integration_uri      = each.value.target_url
  payload_format_version = "1.0"
}

resource "aws_apigatewayv2_route" "routes" {
  for_each         = { for route in local.routes : route.name => route }
  api_id           = aws_apigatewayv2_api.http_api.id
  route_key        = "ANY ${each.value.path}"
  target           = "integrations/${aws_apigatewayv2_integration.routes[each.key].id}"
  authorization_type = "CUSTOM"
  authorizer_id    = aws_apigatewayv2_authorizer.lambda_auth.id
}

resource "aws_apigatewayv2_stage" "default" {
  api_id      = aws_apigatewayv2_api.http_api.id
  name        = "$default"
  auto_deploy = true
}