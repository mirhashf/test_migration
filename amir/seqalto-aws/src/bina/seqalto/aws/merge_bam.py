#!/usr/bin/python2.6
'''
Created on Sep 28, 2011
Script to spin up EC2 instances to download bam files and merge them

@author: amir@binatechnologies.com (Amirhossein Kiani)
'''
import boto.ec2, time, paramiko, sys
from config import config
from utils import execute_bash

create_instance =False 
ec2 = None
instances = []
number_of_instances = 1

WAIT_TIME_AFTER_INSTANCE_CREATION = 30

   
def s3_merge_bam(args):
    print "Args:" + args

    if create_instance:
        # connect to the us-east-1 region
        ec2 = boto.ec2.regions()[1].connect()
        image_id = 'ami-3f34f656'
        image_name = 'oneiric'
        new_reservation = ec2.run_instances(
            image_id=image_id,
            key_name='amir-bina',
            security_groups=['default'], instance_type='m2.4xlarge')
        
        instance = new_reservation.instances[0]
        instances.append(instance)
        
        # wait for instance to boot up 
        print "Spinning up instance for '%s' - %s. Waiting for it to boot up." % (image_id, image_name)
        while instance.state != 'running':
            sys.stdout.write(".")
            time.sleep(1)
            instance.update()
        
        # to allow for ssh server to load
        time.sleep(WAIT_TIME_AFTER_INSTANCE_CREATION)
        
        print ""
        print "Instance is running, ip: %s" % instance.ip_address
        print "Connecting to %s as user %s" % (instance.ip_address, 'ubuntu')
        ip_address = instance.ip_address
    
    ssh = paramiko.SSHClient()
    try:
        # connect to instance through ssh
        #ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        #ssh.connect(ip_address, username='ubuntu', key_filename=os.path.expanduser('~/Downloads/amir-bina.pem'))
        # ssh.connect(ip_address)
        ssh = ""
        execute_bash(ssh, "merge_bam.sh", break_on_stderr=True, print_output=True, config=config, args = args)
        print "Finisehd peacefully."
        
    finally:
        #terminate session and instance
        ssh.close()
        print "Closing session."
        #instance.terminate()
    
    
if __name__ == "__main__":
    print sys.argv
    if len(sys.argv) > 1:
        args = ""
        for arg in sys.argv[1:]:
            args = args +" "+ arg 
        s3_merge_bam(args) 
    else:
        print "Usage: seqalto_cloud.py s3://seqalto/reads1 s3://seqalto/reads2"
    
    
        
    
