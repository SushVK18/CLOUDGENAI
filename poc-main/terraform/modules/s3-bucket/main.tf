module "s3_bucket" {
  source = "terraform-aws-modules/s3-bucket/aws"
  version = "3.14.0"

  bucket = var.bucket_name
  acl    = "null"

  control_object_ownership = true
  object_ownership         = "BucketOwnerEnforced"

  # S3 bucket-level Public Access Block configuration
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true

  versioning = {
    enabled = true
  }

  tags = var.tags
}