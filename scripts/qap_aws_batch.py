import re
import boto3

AWS_REGION='us-east-1'
PROJECT='abide'
KEY='/home/anibalsolon/.ssh/id_rsa.pub'

ec2 = boto3.client('ec2')
iam = boto3.client('iam')
batch = boto3.client('batch')


def tag_resource(client, resource, **kwargs):
    tags = [
        { 'Key': k, 'Value': v }
        for k, v in kwargs.items()
    ]
    return client.create_tags(
        Resources=[resource],
        Tags=tags
    )

def create():

    vpc = ec2.create_vpc(CidrBlock='10.0.0.0/16')
    tag_resource(ec2, vpc['Vpc']['VpcId'], Project=PROJECT)

    internet_gateway = ec2.create_internet_gateway()
    tag_resource(ec2, internet_gateway['InternetGateway']['InternetGatewayId'], Project=PROJECT)

    ec2.attach_internet_gateway(
        InternetGatewayId=internet_gateway['InternetGateway']['InternetGatewayId'],
        VpcId=vpc['Vpc']['VpcId']
    )

    subnet = ec2.create_subnet(
        VpcId=vpc['Vpc']['VpcId'],
        CidrBlock='10.0.1.0/24'
    )
    tag_resource(ec2, subnet['Subnet']['SubnetId'], Project=PROJECT)

    route_table = ec2.create_route_table(
        VpcId=vpc['Vpc']['VpcId']
    )
    tag_resource(ec2, route_table['RouteTable']['RouteTableId'], Project=PROJECT)

    ec2.associate_route_table(
        RouteTableId=route_table['RouteTable']['RouteTableId'],
        SubnetId=subnet['Subnet']['SubnetId']
    )

    route = ec2.create_route(
        DestinationCidrBlock='0.0.0.0/0',
        GatewayId=internet_gateway['InternetGateway']['InternetGatewayId'],
        RouteTableId=route_table['RouteTable']['RouteTableId'],
    )

    security_group = ec2.create_security_group(
        Description=PROJECT,
        GroupName=PROJECT,
        VpcId=vpc['Vpc']['VpcId']
    )
    tag_resource(ec2, security_group['GroupId'], Project=PROJECT)

    ec2.authorize_security_group_ingress(
        GroupId=security_group['GroupId'],
        CidrIp='0.0.0.0/0',
        FromPort=22,
        ToPort=22,
        IpProtocol='tcp',
    )


    iam.create_role(
        RoleName='%s_instance_role' % PROJECT,
        Path='/%s/' % PROJECT,
        AssumeRolePolicyDocument=re.sub('\s', '', """
        {
            "Version": "2012-10-17",
            "Statement": [
                {
                "Action": "sts:AssumeRole",
                "Effect": "Allow",
                "Principal": {
                    "Service": "ec2.amazonaws.com"
                }
                }
            ]
        }
        """)
    )

    iam.attach_role_policy(
        RoleName='%s_instance_role' % PROJECT,
        PolicyArn='arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role'
    )

    instance_role_policy = iam.create_policy(
        PolicyName='%s_instance_role_policy' % PROJECT,
        Path='/%s/' % PROJECT,
        PolicyDocument=re.sub('\s', '', """
        {
            "Version": "2012-10-17",
            "Statement": [
                {
                "Effect": "Allow",
                "Action": [
                    "s3:ListAllMyBuckets",
                    "s3:ListBucket",
                    "s3:PutObject",
                    "s3:GetObject",
                    "s3:DeleteObject"
                ],
                "Resource": ["arn:aws:s3:::*"]
                }
            ]
        }
        """)
    )

    iam.attach_role_policy(
        RoleName='%s_instance_role' % PROJECT,
        PolicyArn=instance_role_policy['Policy']['Arn']
    )

    instance_profile = iam.create_instance_profile(
        InstanceProfileName='%s_instance_profile' % PROJECT,
        Path='/%s/' % PROJECT,
    )

    iam.add_role_to_instance_profile(
        InstanceProfileName='%s_instance_profile' % PROJECT,
        RoleName='%s_instance_role' % PROJECT
    )

    batch_role = iam.create_role(
        RoleName='%s_batch_role' % PROJECT,
        Path='/%s/' % PROJECT,
        AssumeRolePolicyDocument=re.sub('\s', '', """
        {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Action": "sts:AssumeRole",
                    "Effect": "Allow",
                    "Principal": {
                    "Service": "batch.amazonaws.com"
                    }
                }
            ]
        }
        """)
    )

    iam.attach_role_policy(
        RoleName='%s_batch_role' % PROJECT,
        PolicyArn='arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole'
    )

    key_pair = ec2.import_key_pair(
        KeyName=PROJECT,
        PublicKeyMaterial=open(KEY, 'rb').read()
    )


    batch.create_compute_environment(
        computeEnvironmentName=PROJECT,
        type='MANAGED',
        state='ENABLED',
        serviceRole=batch_role['Role']['Arn'],
        computeResources={
            'type': 'EC2',
            'minvCpus': 0,
            'maxvCpus': 0,
            'desiredvCpus': 0,
            'instanceTypes': ['optimal'],
            'subnets': [subnet['Subnet']['SubnetId']],
            'securityGroupIds': [security_group['GroupId']],
            'ec2KeyPair': PROJECT,
            'instanceRole': instance_profile['InstanceProfile']['Arn'],
            'tags': {
                'Project': PROJECT
            },
        }
    )

    batch.create_job_queue(
        jobQueueName=PROJECT,
        state='ENABLED',
        priority=1,
        computeEnvironmentOrder=[
            {
                'order': 1,
                'computeEnvironment': PROJECT
            },
        ]
    )

    batch.register_job_definition(
        jobDefinitionName=PROJECT,
        type='container',
        parameters={
            'string': 'string'
        },
        containerProperties={
            "image": "bids/qap",
            "memory": 8192,
            "vcpus": 4
        }
    )


def destroy():

    if batch.describe_job_definitions(
        jobDefinitions=[PROJECT]
    )['jobDefinitions']:
        batch.deregister_job_definition(
            jobDefinition=PROJECT,
        )

    if batch.describe_job_queues(
        jobQueues=[PROJECT]
    )['jobQueues']:
        batch.delete_job_queue(
            jobQueue=PROJECT,
        )

    if batch.describe_compute_environments(
        computeEnvironments=[PROJECT]
    )['computeEnvironments']:

        batch.update_compute_environment(
            computeEnvironment=PROJECT,
            state='DISABLED',
        )

        batch.delete_compute_environment(
            computeEnvironment=PROJECT,
        )

    ec2.delete_key_pair(KeyName=PROJECT)

    for instance_profile in iam.list_instance_profiles(
        PathPrefix='/%s' % PROJECT
    )['InstanceProfiles']:

        for role in instance_profile['Roles']:

            iam.remove_role_from_instance_profile(
                InstanceProfileName=instance_profile['InstanceProfileName'],
                RoleName=role['RoleName']
            )

        iam.delete_instance_profile(
            InstanceProfileName=instance_profile['InstanceProfileName']
        )


    for role in iam.list_roles(
        PathPrefix='/%s' % PROJECT
    )['Roles']:

        for policy in iam.list_attached_role_policies(
            RoleName=role['RoleName']
        )['AttachedPolicies']:
        
            iam.detach_role_policy(
                RoleName=role['RoleName'],
                PolicyArn=policy['PolicyArn']
            )

            if not policy['PolicyArn'].startswith('arn:aws:iam::aws:'):
                iam.delete_policy(
                    PolicyArn=policy['PolicyArn']
                )

        iam.delete_role(
            RoleName=role['RoleName']
        )


    for security_group in ec2.describe_security_groups(
        Filters=[{ 'Name': 'tag:Project', 'Values': [PROJECT] }]
    )['SecurityGroups']:

        ec2.delete_security_group(
            GroupId=security_group['GroupId']
        )


    for route_table in ec2.describe_route_tables(
        Filters=[{ 'Name': 'tag:Project', 'Values': [PROJECT] }]
    )['RouteTables']:

        for association in route_table['Associations']:

            if association['Main']:
                continue

            ec2.disassociate_route_table(
                AssociationId=association['RouteTableAssociationId']
            )

        ec2.delete_route_table(
            RouteTableId=route_table['RouteTableId']
        )


    for subnet in ec2.describe_subnets(
        Filters=[{ 'Name': 'tag:Project', 'Values': [PROJECT] }]
    )['Subnets']:
        ec2.delete_subnet(
            SubnetId=subnet['SubnetId'],
        )


    for internet_gateway in ec2.describe_internet_gateways(
        Filters=[{ 'Name': 'tag:Project', 'Values': [PROJECT] }]
    )['InternetGateways']:

        for attachment in internet_gateway['Attachments']:

            ec2.detach_internet_gateway(
                InternetGatewayId=internet_gateway['InternetGatewayId'],
                VpcId=attachment['VpcId']
            )

        ec2.delete_internet_gateway(
            InternetGatewayId=internet_gateway['InternetGatewayId'],
        )


    for vpc in ec2.describe_vpcs(
        Filters=[{ 'Name': 'tag:Project', 'Values': [PROJECT] }]
    )['Vpcs']:

        ec2.delete_vpc(
            VpcId=vpc['VpcId'],
        )


create()
destroy()