data "aws_caller_identity" "current" {}
data "aws_region" "current" {}

locals {
  documents_bucket_name = "${var.name_prefix}-${data.aws_caller_identity.current.account_id}-${data.aws_region.current.name}-documents"
  embedding_model_arn   = "arn:aws:bedrock:${data.aws_region.current.name}::foundation-model/${var.embedding_model_id}"
  generation_model_arn  = "arn:aws:bedrock:${data.aws_region.current.name}::foundation-model/${var.generation_model_id}"
}

resource "aws_s3_bucket" "documents" {
  bucket = local.documents_bucket_name
  tags   = var.tags
}

resource "aws_s3_bucket_public_access_block" "documents" {
  bucket                  = aws_s3_bucket.documents.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_server_side_encryption_configuration" "documents" {
  bucket = aws_s3_bucket.documents.id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_versioning" "documents" {
  bucket = aws_s3_bucket.documents.id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_iam_role" "bedrock_knowledge_base" {
  name = "${var.name_prefix}-bedrock-kb"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Principal = {
        Service = "bedrock.amazonaws.com"
      }
      Action = "sts:AssumeRole"
      Condition = {
        StringEquals = {
          "aws:SourceAccount" = data.aws_caller_identity.current.account_id
        }
        ArnLike = {
          "aws:SourceArn" = "arn:aws:bedrock:${data.aws_region.current.name}:${data.aws_caller_identity.current.account_id}:knowledge-base/*"
        }
      }
    }]
  })
  tags = var.tags
}

resource "aws_iam_role_policy" "bedrock_knowledge_base" {
  name = "${var.name_prefix}-bedrock-kb"
  role = aws_iam_role.bedrock_knowledge_base.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect   = "Allow"
        Action   = ["bedrock:InvokeModel"]
        Resource = local.embedding_model_arn
      },
      {
        Effect = "Allow"
        Action = ["s3:ListBucket"]
        Resource = [
          aws_s3_bucket.documents.arn
        ]
      },
      {
        Effect   = "Allow"
        Action   = ["s3:GetObject"]
        Resource = ["${aws_s3_bucket.documents.arn}/*"]
      },
      {
        Effect   = "Allow"
        Action   = ["aoss:APIAccessAll"]
        Resource = [aws_opensearchserverless_collection.vectors.arn]
      }
    ]
  })
}

resource "aws_opensearchserverless_security_policy" "encryption" {
  name = "${var.name_prefix}-encryption"
  type = "encryption"
  policy = jsonencode({
    Rules = [{
      Resource     = ["collection/${var.name_prefix}-vectors"]
      ResourceType = "collection"
    }]
    AWSOwnedKey = true
  })
}

resource "aws_opensearchserverless_security_policy" "network" {
  name = "${var.name_prefix}-network"
  type = "network"
  policy = jsonencode([{
    Description = "Allow Bedrock service access to the vector collection"
    Rules = [{
      Resource     = ["collection/${var.name_prefix}-vectors"]
      ResourceType = "collection"
    }]
    AllowFromPublic = false
    SourceServices  = ["bedrock.amazonaws.com"]
  }])
}

resource "aws_opensearchserverless_access_policy" "bedrock" {
  name = "${var.name_prefix}-bedrock"
  type = "data"
  policy = jsonencode([{
    Description = "Allow the Bedrock knowledge base to manage its vector index"
    Principal   = [aws_iam_role.bedrock_knowledge_base.arn]
    Rules = [
      {
        Resource = ["collection/${var.name_prefix}-vectors"]
        Permission = [
          "aoss:DescribeCollectionItems"
        ]
        ResourceType = "collection"
      },
      {
        Resource = ["index/${var.name_prefix}-vectors/*"]
        Permission = [
          "aoss:CreateIndex",
          "aoss:DescribeIndex",
          "aoss:ReadDocument",
          "aoss:UpdateIndex",
          "aoss:WriteDocument"
        ]
        ResourceType = "index"
      }
    ]
  }])
}

resource "aws_opensearchserverless_collection" "vectors" {
  name             = "${var.name_prefix}-vectors"
  type             = "VECTORSEARCH"
  standby_replicas = "DISABLED"
  tags             = var.tags

  depends_on = [
    aws_opensearchserverless_security_policy.encryption,
    aws_opensearchserverless_security_policy.network
  ]
}

resource "aws_bedrockagent_knowledge_base" "this" {
  name        = "${var.name_prefix}-knowledge-base"
  description = "RAG knowledge base backed by the project document bucket"
  role_arn    = aws_iam_role.bedrock_knowledge_base.arn
  tags        = var.tags

  knowledge_base_configuration {
    type = "VECTOR"

    vector_knowledge_base_configuration {
      embedding_model_arn = local.embedding_model_arn
    }
  }

  storage_configuration {
    type = "OPENSEARCH_SERVERLESS"

    opensearch_serverless_configuration {
      collection_arn    = aws_opensearchserverless_collection.vectors.arn
      vector_index_name = "${var.name_prefix}-index"

      field_mapping {
        vector_field   = "bedrock-knowledge-base-default-vector"
        text_field     = "AMAZON_BEDROCK_TEXT_CHUNK"
        metadata_field = "AMAZON_BEDROCK_METADATA"
      }
    }
  }

  depends_on = [
    aws_iam_role_policy.bedrock_knowledge_base,
    aws_opensearchserverless_access_policy.bedrock
  ]
}

resource "aws_bedrockagent_data_source" "documents" {
  knowledge_base_id    = aws_bedrockagent_knowledge_base.this.id
  name                 = "${var.name_prefix}-documents"
  data_deletion_policy = "RETAIN"

  data_source_configuration {
    type = "S3"

    s3_configuration {
      bucket_arn         = aws_s3_bucket.documents.arn
      inclusion_prefixes = ["documents/"]
    }
  }
}

resource "aws_iam_role" "application" {
  name = "${var.name_prefix}-application"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect = "Allow"
      Principal = {
        Service = "ecs-tasks.amazonaws.com"
      }
      Action = "sts:AssumeRole"
    }]
  })
  tags = var.tags
}

resource "aws_iam_role_policy" "application" {
  name = "${var.name_prefix}-bedrock-retrieval"
  role = aws_iam_role.application.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "bedrock:Retrieve",
          "bedrock:RetrieveAndGenerate"
        ]
        Resource = [aws_bedrockagent_knowledge_base.this.arn]
      },
      {
        Effect   = "Allow"
        Action   = ["bedrock:InvokeModel"]
        Resource = [local.generation_model_arn]
      }
    ]
  })
}
