'''
Created on Sep 30, 2011

@author: amir
'''
from string import Template
import datetime

def execute_ssh_command(ssh, command, break_on_stderr=True, print_output=True, dir='', print_command=True):
    '''
    Execute ssh command on paramiko session and print errors/break if 
    unsuccessful.
    '''
    before = datetime.datetime.now()
    if(print_command):
        print "["+ str(before) +"] Executing: " + command
    else:
        print "["+ str(before) +"] Executing: hidden command."
    
    stdin, stdout, stderr = ssh.exec_command("cd " + dir + ";" + command)
    if break_on_stderr:
        error = stderr.readlines()
        if len(error) != 0:
            print "STDOUT:\n%s" % "".join(stdout.readlines())
            raise Exception("".join(error))     
    elif print_output:
            print "".join(stdout.readlines())
    else:
        return stdin, stdout, stderr


def execute_bash(ssh, bash_file_path, break_on_stderr=False, print_output=False, config={}, args=''):
    file = open(bash_file_path)
    content = ""
    while 1:
        line = file.readline()
        content = content + line
        if not line:
            break
        template = Template(content)
    execute_ssh_command(ssh, 'echo "'  + template.substitute(config).replace('"', '\\"') + '"> ~/script.sh ' + args + ' && chmod +x script.sh && nohup ./script.sh 1>std.out 2>std.err &' , break_on_stderr, print_output)
