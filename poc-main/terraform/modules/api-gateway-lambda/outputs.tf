output "api_endpoint" {
  description = "The endpoint of the HTTP API ."
  value = aws_apigatewayv2_api.http_api.api_endpoint
}

output "lambda_function_arn" {
  description = "The ARN of the lambda function."
  value = aws_lambda_function.authorizer.arn
}