#!/usr/bin/env python

import boto3
import argparse
import sys

AWS_ACCESS_KEY = "AKIAIVG5XVWFKGNZE47Q"
AWS_SECRET_KEY = "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"
AWS_REGION     = "us-west-1"

AWS_SESSION = boto3.session.Session(region_name = AWS_REGION, aws_access_key_id = AWS_ACCESS_KEY, aws_secret_access_key = AWS_SECRET_KEY)
session_ec2 = AWS_SESSION.resource('ec2')

def launch_pgp_instance():
    print "Launching i2.xlarge encryption/decryption instance"

    pgpInstance = session_ec2.create_instances(
        DryRun=False,
        ImageId = "ami-d1315fb1",
        MinCount = 1,
        MaxCount = 1,
        KeyName = "henryic@PIGEON",
        InstanceType='i2.xlarge',
        Placement={
            'AvailabilityZone': 'us-west-1b',
            'Tenancy': 'dedicated'
        },
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/sda1',
                'Ebs': {
                    'VolumeSize': 16,
                    'DeleteOnTermination': True,
                    'VolumeType': 'gp2'
                }
            },
            {
                'VirtualName': 'ephemeral0',
                'DeviceName': '/dev/sdb'
            }
        ],
        Monitoring={
            'Enabled': False
        },
        DisableApiTermination=False,
        InstanceInitiatedShutdownBehavior='stop',
        NetworkInterfaces=[
            {
                'DeviceIndex': 0,
                'SubnetId': 'subnet-1abb1243',
                'Description': 'Primary network interface',
                'Groups': [
                    'sg-da0951bf'
                ],
                'DeleteOnTermination': True,
                'AssociatePublicIpAddress': False
            }
        ],
        EbsOptimized=False
    )

    session_ec2.create_tags(Resources=[pgpInstance[0].id], Tags=[{'Key':'Name','Value':'VA_VPCPrivate_PGPInstance'}])

# end def launch_pgp_instance()

def stop_pgp_instance():
    print "Stop i2.xlarge encryption/decryption instance here"

def terminate_pgp_instance():
    print "Terminate i2.xlarge encryption/decryption instance here"


def launch_gateway_instance():
    print "Launching m3.medium gateway instance"

    sshGatewayInstance = session_ec2.create_instances(
        DryRun=False,
        ImageId = "ami-05cf2541",
        MinCount = 1,
        MaxCount = 1,
        KeyName = "henryic@PIGEON",
        InstanceType='m3.medium',
        Placement={
            'Tenancy': 'dedicated',
            'AvailabilityZone': 'us-west-1b'
        },
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/xvda',
                'Ebs': {
                    'VolumeSize': 8,
                    'DeleteOnTermination': True,
                    'VolumeType': 'gp2'
                }
            },
            {
                'VirtualName': 'ephemeral0',
                'DeviceName': '/dev/sdb'
            }
        ],
        Monitoring={
            'Enabled': False
        },
        DisableApiTermination=False,
        InstanceInitiatedShutdownBehavior='stop',
        NetworkInterfaces=[
            {
                'DeviceIndex': 0,
                'SubnetId': 'subnet-57fb550e',
                'Description': 'Public Subnet Interface',
                'PrivateIpAddress': '10.1.0.64',
                'Groups': [
                    'sg-bed7fddb',
                    'sg-22d0fe47'
                ],
                'DeleteOnTermination': True
            },
            {
                'DeviceIndex': 1,
                'SubnetId': 'subnet-1abb1243',
                'Description': 'Private Subnet Interface',
                'PrivateIpAddress': '10.1.1.64',
                'Groups': [
                    'sg-32d0fe57'
                ],
                'DeleteOnTermination': True
            },
        ],
        EbsOptimized=False
    )

    session_ec2.create_tags(Resources=[sshGatewayInstance[0].id], Tags=[{'Key':'Name','Value':'VA_VPCPublicSSHGateway'}])

    session_ec2_client.associate_address(
        DryRun = False,
        AllocationId = 'eipalloc-61ca0104',
        NetworkInterfaceId = sshGatewayInstance[0].network_interfaces[0].id,
        AllowReassociation = True
    )

    session_ec2.instances.filter(Filters=[{'Name':'tag:Name', 'Values':['VA_VPCPublicSSHGateway']}])

# end def launch_gateway_instance()

def stop_gateway_instance():
    print "Stop m3.medium gateway instance"

def terminate_gateway_instance():
    print "Terminate m3.medium gateway instance"


def launch_endpoint_instance():
    print "Launching m3.xlarge endpoint instance"
    dcEndpointInstance = session_ec2.create_instances(
        DryRun=False,
        ImageId = "ami-9ee593fe",
        MinCount = 1,
        MaxCount = 1,
        KeyName = "henryic@PIGEON",
        InstanceType='m3.xlarge',
        Placement={
            'Tenancy': 'dedicated',
            'AvailabilityZone': 'us-west-1b'
        },
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/xvda',
                'Ebs': {
                    'VolumeSize': 8,
                    'DeleteOnTermination': True,
                    'VolumeType': 'gp2'
                }
            },
            {
                'VirtualName': 'ephemeral0',
                'DeviceName': '/dev/sdb'
            }
        ],
        Monitoring={
            'Enabled': False
        },
        DisableApiTermination=False,
        InstanceInitiatedShutdownBehavior='stop',
        NetworkInterfaces=[
            {
                'DeviceIndex': 0,
                'SubnetId': 'subnet-1abb1243',
                'Description': 'Private Subnet Interface',
                'PrivateIpAddress': '10.1.1.65',
                'Groups': [
                    'sg-32d0fe57',
                    'sg-da0951bf'
                ],
                'DeleteOnTermination': True
            },
        ],
        EbsOptimized=False
    )

    session_ec2.create_tags(Resources=[dcEndpointInstance[0].id], Tags=[{'Key':'Name','Value':'VA_VPCPrivate_DCEndpoint'}])

# end def launch_endpoint_instance()

def stop_endpoint_instance():
    print "Stopping m3.xlarge endpoint instance"

    for instance in session_ec2.instances.filter(Filters=[{'Name':'tag:Name', 'Values':['VA_VPCPrivate_DCEndpoint']}, {'Name': 'instance-state-name', 'Values': ['running']}]):
        session_ec2.Instance(instance.id).stop()

def terminate_endpoint_instance():
    print "Terminating m3.xlarge endpoint instance"

    for instance in session_ec2.instances.filter(Filters=[{'Name':'tag:Name', 'Values':['VA_VPCPrivate_DCEndpoint']}, {'Name': 'instance-state-name', 'Values': ['running','stopped']}]):
        session_ec2.Instance(instance.id).terminate()


def main():
    arg_parser = argparse.ArgumentParser("Perform AWS operations for VA processing.")

    arg_parser.add_argument("function",  help = "work with a <gateway>, <endpoint>, or <pgp>-type instance")
    arg_parser.add_argument("action",    help = "<start>, <stop>, or <terminate> instance(s)")

    args = arg_parser.parse_args()

    if args.function == 'gateway':
        if args.action == 'start':
            launch_gateway_instance()
        elif args.action == 'stop':
            stop_gateway_instance()
        elif args.action == 'terminate':
            terminate_gateway_instance()
        else:
            sys.exit("Unrecognized action '%s'" % (args.action))
    elif args.function == 'endpoint':
        if args.action == 'start':
            launch_endpoint_instance()
        elif args.action == 'stop':
            stop_endpoint_instance()
        elif args.action == 'terminate':
            terminate_endpoint_instance()
        else:
            sys.exit("Unrecognized action '%s'" % (args.action))
    elif args.function == 'pgp':
        if args.action == 'start':
            launch_pgp_instance()
        elif args.action == 'stop':
            stop_pgp_instance()
        elif args.action == 'terminate':
            terminate_pgp_instance()
        else:
            sys.exit("Unrecognized action '%s'" % (args.action))
    else:
        sys.exit("Unrecognized function '%s'" % (args.function))


if __name__ == "__main__":
    main()


