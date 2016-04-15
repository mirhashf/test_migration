#!/usr/bin/env python

import boto3
import argparse
import sys

import json,subprocess

from aws_credentials import AWS_ACCESS_KEY, AWS_SECRET_KEY

AWS_REGION     = "us-west-1"

VA_ROUTE_TABLE = "rtb-0c27ab69"

AWS_IP_RANGE_URL = "https://ip-ranges.amazonaws.com/ip-ranges.json"

AWS_SESSION = boto3.session.Session(region_name = AWS_REGION, aws_access_key_id = AWS_ACCESS_KEY, aws_secret_access_key = AWS_SECRET_KEY)
session_ec2 = AWS_SESSION.resource('ec2')

def update_route_table():
    FIXED_ROUTES = ["10.0.0.0/24", "169.254.253.0/30", "10.1.0.0/16"]
    NAT_GATEWAY_ID = "nat-08f5a26287c1a2d14"

    # dodgy thing to get around urllib2 https handling error
    aws_ip_range_json = json.loads(subprocess.Popen(['curl', AWS_IP_RANGE_URL], stdout=subprocess.PIPE).communicate()[0])

    region_ips = [str(IP['ip_prefix']) for IP in aws_ip_range_json['prefixes'] if IP['region'] in [AWS_REGION, 'GLOBAL'] and IP['service'] in ['AMAZON']]
    route_table_dests = [ROUTE.destination_cidr_block for ROUTE in session_ec2.RouteTable(VA_ROUTE_TABLE).routes]

    routes_to_add = set(region_ips) - set(route_table_dests)
    routes_to_del = set(route_table_dests) - set(region_ips) - set(FIXED_ROUTES)

    for IP in routes_to_del:
        session_ec2.Route(route_table_id=VA_ROUTE_TABLE, destination_cidr_block=IP).delete()

    for IP in routes_to_add:
        session_ec2.RouteTable(VA_ROUTE_TABLE).create_route(DestinationCidrBlock=IP, NatGatewayId=NAT_GATEWAY_ID)


def launch_pgp_instance():
    print "Launching i2.xlarge encryption/decryption instance"

    pgp_seq_num = 0

    pgpInstance = session_ec2.create_instances(
        DryRun=False,
        ImageId = "ami-61166901",
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
                'PrivateIpAddress': '10.1.1.' + str(80 + pgp_seq_num),
                'Groups': [
                    'sg-da0951bf'
                ],
                'DeleteOnTermination': True,
                'AssociatePublicIpAddress': False
            }
        ],
        EbsOptimized=False
    )

    session_ec2.create_tags(Resources=[pgpInstance[0].id], Tags=[{'Key':'Name','Value':'VA_VPCPrivate_PGPInstance'}, {'Key':'SequenceNumber','Value':str(pgp_seq_num)}])

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
        EbsOptimized=False,
        UserData='#!/bin/bash\nsed -i \'s/nobootwait/nofail/\' /etc/fstab',
        IamInstanceProfile={
            'Arn':'arn:aws:iam::540758486469:instance-profile/VAVPCServer'
        }

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
        ImageId = "ami-64146b04",
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
        EbsOptimized=False,
        UserData='#!/bin/bash\nsed -i \'s/nobootwait/nofail/\' /etc/fstab\necho "interface=eth0" > /etc/dnsmasq.conf\nservice dnsmasq restart',
        IamInstanceProfile={
            'Arn':'arn:aws:iam::540758486469:instance-profile/VAVPCServer'
        }
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


def describe_instances(instance_function):
    if instance_function == 'all':
        functions_filter = []
    elif instance_function == 'gateway':
        functions_filter = [{'Name':'tag:Name', 'Values':['VA_VPCPublicSSHGateway']}]
    elif instance_function == 'endpoint':
        functions_filter = [{'Name':'tag:Name', 'Values':['VA_VPCPrivate_DCEndpoint']}]
    elif instance_function == 'pgp':
        functions_filter = [{'Name':'tag:Name', 'Values':['VA_VPCPrivate_PGPInstance']}]
    else:
        sys.exit("Unknown instance function: '%s'" % (instance_function))

    print "Instance\tState\tIP\tTags"
    for instance in session_ec2.instances.filter(Filters=functions_filter):
        print instance.id + "\t" + instance.state['Name'] + "\t" + instance.private_ip_address + "\t" +  ";".join(str(tag) for tag in instance.tags)


def main():
    arg_parser = argparse.ArgumentParser("Perform AWS operations for VA processing.")

    arg_parser.add_argument("function",  help = "work with a <gateway>, <endpoint>, or <pgp>-type instance")
    arg_parser.add_argument("action",    help = "<start>, <stop>, <terminate>, or <describe> [function]-type instance")

    args = arg_parser.parse_args()

    if args.action == 'describe':
        describe_instances(args.function)
    else:
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


