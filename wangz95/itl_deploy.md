# ITL Deployment

## Overview

The objective of this document is describe how to do ITL deployment.

## Docker Deployment

### Configure Maven

On your laptop
```bash
scp -p t-rex:~/.m2/settings.xml ~/.m2/settings.xml

sed -i -e 's/t-rex/localhost/g' ~/.m2/settings.xml
```

### Create a Git tag for test

Clone the repo to your laptop
```bash
git clone git@github.com:BinaTechnologies/seqalto.git
```

Create a Git tag for test, the format of the tag is `test-GMS-YYYYMMDD(-[A-Z])*`.
```bash
cd seqalto

git pull

git tag -am test-GMS-20160520 test-GMS-20160520

git push origin test-GMS-20160520
```

### Set up SSH tunnel

On your laptop, forward local port 8081 to remote port 8081 on **t-rex**

```bash
ssh -L 8081:localhost:8081 t-rex -A
```
**NOTE:** This is required before build because we have a Nexus repository manager running on **t-rex**.

### Build with Maven

Build with **internal** profile and skip tests
```bash
mvn clean install -DskipTests -P internal
```
`test/integration/target/assembly/test-integration-launcher.jar` should be generated at this point.

### Deploy

Start a new tmux session on your laptop
```bash
tmux
```

Tear down current active ITL deployment running on target Docker containers
```bash
java -cp test/integration/target/assembly/test-integration-launcher.jar com.bina.seqalto.test.integration.IntegrationTest_Teardown --rigtype "BinaDevelopment" --rigargs "--gitURL git@github.com:BinaTechnologies/seqalto.git test-GMS-20160406 tehran-01-06 tehran-01-07 tehran-01-08 tehran-01-09" --username wangz95 --noprompt
```

Deploy specified Git tag (test-GMS-20160520 in this case) on target Docker containers
```bash
java -jar test/integration/target/assembly/test-integration-launcher.jar --rigtype "BinaDevelopment" --rigargs "--gitURL git@github.com:BinaTechnologies/seqalto.git test-GMS-20160520 tehran-01-06 tehran-01-07 tehran-01-08 tehran-01-09" --username wangz95 --noprompt --keepup
```

Add following `rigargs` if this is a ctDNA ITL deployment
```bash
--seedEnabled CtdnaLocalTestData
```

## AWS Deployment

### Set up SSH tunnel

On your laptop, forward local port 8081 to remote port 8081 on **t-rex**
```bash
ssh -L 8081:localhost:8081 t-rex -A
```

### Set up reverse SSH tunnel

Open a new tab in Terminal, set up reverse SSH tunnel as below
```bash
ssh -tt -R 8081:localhost:8081 admin@ext-itl-00.qe.bina.com ssh -R 8081:localhost:8081 build-00
```

### Abort current active ITL deployment on AWS host

Open another tab in Terminal, SSH to AWS host `ext-itl-00.qe.bina.com` as `admin` user
```bash
ssh admin@ext-itl-00.qe.bina.com -A
```
**NOTE:** This is possible once your public key is added to `/home/admin/.ssh/authorized_keys`.

Abort current active ITL deployment
```bash
tmux attach

Ctrl-C

exit
```

### Deploy

SSH to AWS host `ext-itl-00.qe.bina.com` as `admin` user
```bash
ssh admin@ext-itl-00.qe.bina.com -A
```

Start a new tmux session
```bash
tmux
```

Deploy specified Git tag (test-GMS-20160520 in this case)
```bash
java -jar ./aws-itl/test-integration-launcher.jar --rigtype BinaDevelopment --rigargs "--port 18933 --gitURL git@github.com:BinaTechnologies/seqalto.git test-GMS-20160520 172.17.0.1" --username admin --noprompt
```

Detach from the current tmux session
```
Ctrl-b

d
```

You can always come back using
```bash
tmux attach
```

## Bina Box Deployment

### Send out a notification

Send out an Email to people about coming Bina Box deployment.

### SSH to build-00

SSH to **t-rex** from your laptop
```bash
ssh t-rex -A
```

Start a new tmux session
```bash
tmux
```

SSH to **build-00** from **t-rex**
```bash
ssh build-00
```

Add following to `~/.ssh/config` on **build-00**
```
Host pune-0? lima-0? nara-0? maui-0? oslo-0? kiev hilo reno-0? science
  User binatech
```

You also need to add your public key on **build-00** to `~/.ssh/authorized_keys` on target Bina Boxes
* **pune-0?**
* **lima-0?**
* **reno-0?**

### Prepare repo for deployment

Clone the repo to **build-00**
```bash
git clone https://github.com/BinaTechnologies/seqalto.git
```

Checkout desired Git tag
```bash
git checkout test-GMS-20160520
```

Download 3rd-party packages using the Python script
```bash
./third-party.py
```

### Deploy
```bash
cd loomis2/sh

DANGER=1 CLEAN=1 ./release-auto-update pune lima reno
```
where `pune lime reno` are Bina Boxes that you want to deploy to.

