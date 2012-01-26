#!/usr/bin/python2.7
'''
Created on Sep 28, 2011
Script to spin up EC2 instances to run SeqAlto, BWA and GATK.

@author: amir@binatechnologies.com (Amirhossein Kiani)
'''
import boto.ec2, time, paramiko, sys, os
from utils import execute_bash
ec2 = None
instances = []
number_of_instances = 1

WAIT_TIME_AFTER_INSTANCE_CREATION = 60

from config import config

def paired_end_on_aws(lane, sample, library,instance_type, memory, run_bwa=False, create_instance = True, test=False, ip_address=None, s3_out="seqalto"):
    config['s3_out_bucket']=s3_out
    config['lane'] = lane
    config['sample'] = sample
    config['library'] = library
    config['mem'] = memory
    aligner= ""
    if not run_bwa:
        aligner = "SeqAlto"
    else:
        aligner = "BWA"
                    
    print "Running " + aligner + " on lane: " + lane + " ..."
 
    image_id = 'ami-3f34f656'
    image_name = 'oneiric'
    
    if test:
        config['seqalto_genome']='chr21.fa_22.midx'
        config['bwa_genome_bundle']='bwa-index_test.tar.gz'
        config['bwa_genome']='chr21.fa'
    
    if create_instance:
        # connect to the us-east-1 region
        ec2 = boto.ec2.regions()[1].connect()
        new_reservation = ec2.run_instances(
            image_id=image_id,
            key_name='amir-bina',
            security_groups=['default'], instance_type=instance_type)
        
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
        print "Connecting to instance through SSH..."
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ip_address, username='ubuntu', key_filename=os.path.expanduser('~/Downloads/amir-bina.pem'))
        # ssh.connect(ip_address)
        
        execute_bash(ssh, "run_"+ aligner + ".sh", break_on_stderr=False, print_output=False, config=config)
        
        print "Finisehd peacefully."
    finally:
        #terminate session and instance
        ssh.close()
        print "Closing session."
        #instance.terminate()
    

    
if __name__ == "__main__":
    print sys.argv
    if len(sys.argv) > 7:
        lane = sys.argv[1]
        s3_out = sys.argv[2]
        sample = sys.argv[3]
        library = sys.argv[4]
        run_bwa = False
        create_instance = sys.argv[6]
        instance_type = sys.argv[7]
        memory = sys.argv[8]
	paired_end_on_aws(lane,sample,library,instance_type,memory, run_bwa,create_instance=create_instance,s3_out=s3_out)
    else:
        print "Usage: reads_to_bam.py LANE S3_OUT_BUCKET SAMPLE LIBRARY RUN_BWA CREATE_INSTANCE INSTANCE_TYPE MEMORY"
        print "Running TEST MODE!!!"
        lane = "small_pair"
        sample = "snyder"
        library = "library"
        run_bwa = False
        create_instance = False
        test_mode = True
        memory = "5"
        instance_type = "m1.large"
        s3_out = "seqalto-results"
       # paired_end_on_aws(lane,sample,library,instance_type,memory, run_bwa,create_instance=create_instance,test=test_mode,s3_out=s3_out, ip_address="ec2-107-20-46-235.compute-1.amazonaws.com")
    
        
    
