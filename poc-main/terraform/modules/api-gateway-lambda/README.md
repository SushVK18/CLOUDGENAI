# Azure-Apigw-Authorizer-lambda
This lambda function is based on https://github.com/bayer-int/smart-aws-cp-Service-Catalog/blob/main/azure-apigw-authorizer-lambda and contains the files, which are used for the service catalog product 

## Goal
Authenticate against Azure AD in a secure and customizable manner.

## Overview
This product extends the deprecated version of the Azure authorizer.

### Prepare the Lambda Function for Deployment
Build and zip the lambda function if it does not exist yet before referring it in the terraform.
```
cd terraform/modules/api-gateway-lambda/azure-apigw-authorizer-lambda
npm install
zip - r ../azure-apigw-authorizer-lambda.zip .
```

### Product Features
Adds the below custom catalog product parameters to extend security and can be overridden by consumer:
* **Audience**: Token audience, set this for extra security. By default this is not set.
* **Algorithms**: Token signing algorithm. Default RS256.
* **IssuersUrls**: Base URL(s) of the Azure AD tenant. Use https://sts.windows.net/fcb2b37b-5da0-466b-9b83-0014b67a7c78/ for tokens version 1.0 and https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/v2.0 for tokens version 2.0. Default, https://sts.windows.net/fcb2b37b-5da0-466b-9b83-0014b67a7c78/

### Limitations
* **Payload Version**: This authorizer now supports both payload version 1.0 and 2.0. When using payload version 2.0, it is necessary to set the response mode to "IAM policy".