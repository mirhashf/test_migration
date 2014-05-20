#!/bin/sh

nohup java -DLOG_HOME=/home/ec2-user/ap-activity-client/logs -jar /home/ec2-user/winstone/winstone.jar --warfile=ap-activity-client.war --httpPort=9000 > winstone.log 2>&1 &	
