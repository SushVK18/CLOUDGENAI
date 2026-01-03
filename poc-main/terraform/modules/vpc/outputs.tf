output "vpc_id" {
    value = module.vpc.vpc_id
    description = "vpc id"
}

output "private_subnets" {
    value = module.vpc.private_subnets
    description = "List of IDs of private subnets"
}

output "public_subnets" {
    value = module.vpc.public_subnets
    description = "List of IDs of public subnets"
}

output "public_subnets_cidr_blocks" {
    value = module.vpc.public_subnets_cidr_blocks
}

output "vpc_cidr_block" {
    value = module.vpc.vpc_cidr_block
}

output "private_subnets_cidr_blocks" {
    value = module.vpc.private_subnets_cidr_blocks
}