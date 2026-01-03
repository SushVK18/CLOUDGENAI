data "aws_caller_identity" "current" {}

##############################
# 1. Security Groups
##############################

resource "aws_security_group" "alb_app_sg" {
  name   = "${var.project_prefix}-alb-app-sg"
  vpc_id = var.vpc_id

  ingress {
    description = "Allow HTTPS traffic from the internet"
    from_port   = 443
    to_port     = 443
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  egress {
    description = "Allow all outbound traffic"
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = var.tags
}

resource "aws_security_group" "service_app_sg" {
  name   = "${var.project_prefix}-service-app-sg"
  vpc_id = var.vpc_id

  ingress {
    description     = "Allow TCP 8888 from ALB security group"
    from_port       = 8888
    to_port         = 8888
    protocol        = "tcp"
    security_groups = [aws_security_group.alb_app_sg.id]
  }
  
  ingress {
    description     = "Allow TCP 8000 from ALB security group"
    from_port       = 8000
    to_port         = 8000
    protocol        = "tcp"
    security_groups = [aws_security_group.alb_app_sg.id]
  }

  egress {
    description = "Allow all outbound traffic"
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = var.tags
}


##############################
# 2. Target Groups (using for_each)
##############################

resource "aws_lb_target_group" "tg" {
  for_each        = var.target_group_settings
  name            = each.key
  port            = each.value.port
  protocol        = each.value.protocol
  target_type     = "ip"
  vpc_id          = var.vpc_id
  ip_address_type = "ipv4"
  protocol_version = each.value.protocol_version

  health_check {
    protocol = each.value.health_check_protocol
    path     = each.value.health_check_path
  }

  tags = var.tags
}


##############################
# 3. Application Load Balancer
##############################

resource "aws_lb" "this" {
  name               = "${var.project_prefix}-alb"
  internal           = false
  load_balancer_type = "application"
  security_groups    = [aws_security_group.alb_app_sg.id]
  enable_deletion_protection       = false
  enable_cross_zone_load_balancing = true
  ip_address_type    = "ipv4"
  subnets            = var.public_subnets_ids

  tags = var.tags
}


##############################
# 5. ALB Listener & Listener Rules
##############################

resource "aws_lb_listener" "https" {
  load_balancer_arn = aws_lb.this.arn
  port              = 443
  protocol          = "HTTPS"
  # ssl_policy        = "ELBSecurityPolicy-2016-08" 
  
  # Note: you need to have an SSL certificate ARN to attach here.
  # This is a placeholder and must be replaced.
  certificate_arn = var.certificate_arn

  # Provide a default action in case none of the rules match.
  default_action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.tg["app-tf-framework-tg"].arn
  }
}

# Listener Rule 4.1: Forward to app-tf-molmim-tg when host header matches
resource "aws_lb_listener_rule" "molmim" {
  listener_arn = aws_lb_listener.https.arn
  priority     = 10

  condition {
    host_header {
      values = ["molmim-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"]
    }
  }

  action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.tg["app-tf-molmim-tg"].arn
  }
}

# Listener Rule 4.2: Forward to app-tf-diffdock-tg when host header matches
resource "aws_lb_listener_rule" "diffdock" {
  listener_arn = aws_lb_listener.https.arn
  priority     = 20

  condition {
    host_header {
      values = ["diffdock-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"]
    }
  }

  action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.tg["app-tf-diffdock-tg"].arn
  }
}

# Listener Rule 4.3: Forward to app-tf-genmol-tg when host header matches
resource "aws_lb_listener_rule" "genmol" {
  listener_arn = aws_lb_listener.https.arn
  priority     = 30

  condition {
    host_header {
      values = ["genmol-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"]
    }
  }

  action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.tg["app-tf-genmol-tg"].arn
  }
}

# Listener Rule 4.4: Forward to app-tf-framework-tg with OIDC authentication when host header matches
resource "aws_lb_listener_rule" "framework" {
  listener_arn = aws_lb_listener.https.arn
  priority     = 40

  condition {
    host_header {
      values = ["jupyter-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"]
    }
  }

  action {
    type = "authenticate-oidc"

    authenticate_oidc {
      client_id              = var.client_id
      client_secret          = var.client_secret
      authorization_endpoint = "https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/oauth2/v2.0/authorize"
      issuer                 = "https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/v2.0"
      token_endpoint         = "https://login.microsoftonline.com/fcb2b37b-5da0-466b-9b83-0014b67a7c78/oauth2/v2.0/token"
      user_info_endpoint     = "https://graph.microsoft.com/oidc/userinfo"
      scope                  = "openid offline_access"
      session_timeout        = var.session_timeout
    }
  }

  action {
    type             = "forward"
    target_group_arn = aws_lb_target_group.tg["app-tf-framework-tg"].arn
  }
}


##############################
# 6. Route 53 DNS Records
##############################

# Create an alias record for the ALB inference in Route53.
resource "aws_route53_record" "alb_subdomain" {
  zone_id = var.route53_zone_id
  name    = "${var.project_prefix}.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"
  type    = "A"

  alias {
    name                   = "dualstack.${aws_lb.this.dns_name}"
    zone_id                = aws_lb.this.zone_id
    evaluate_target_health = true
  }
}

# Create an alias record for the ALB framework in Route53.
resource "aws_route53_record" "jupyter_subdomain" {
  zone_id = var.route53_zone_id
  name    = "jupyter-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"
  type    = "A"

  alias {
    name                   = "dualstack.${aws_lb.this.dns_name}"
    zone_id                = aws_lb.this.zone_id
    evaluate_target_health = true
  }
}

# Create a CNAME record for the service in Route53
resource "aws_route53_record" "molmim_subdomain" {
  zone_id = var.route53_zone_id
  name    = "molmim-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"
  type    = "CNAME"
  ttl     = 300

  records = [ aws_route53_record.alb_subdomain.name ]
}

resource "aws_route53_record" "diffdock_subdomain" {
  zone_id = var.route53_zone_id
  name    = "diffdock-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"
  type    = "CNAME"
  ttl     = 300

  records = [ aws_route53_record.alb_subdomain.name ]
}

resource "aws_route53_record" "genmol_subdomain" {
  zone_id = var.route53_zone_id
  name    = "genmol-tf.${data.aws_caller_identity.current.account_id}.cloud.bayer.com"
  type    = "CNAME"
  ttl     = 300

  records = [ aws_route53_record.alb_subdomain.name ]
}