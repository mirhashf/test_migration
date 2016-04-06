Setting up Google Cloud CLI on VA Transfer Server
=================================================

Allowing DNS resolution from sneaker
------------------------------------

- Set up dnsmasq on private gateway instance (10.1.1.65), listen on interface eth0 ##
- Add private gateway instance (10.1.1.65) as nameserver on sneaker


Install Google Cloud SDK
------------------------

- On SSH gateway instance (10.1.0.64), follow the Google-provided instructions and run:

>    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)"
>    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | sudo tee /etc/apt/sources.list.d/google-cloud-sdk.list
>    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
>    sudo apt-get update && sudo apt-get install google-cloud-sdk

- Copy /var/cache/apt/archives/google-cloud-sdk_<version>_all.deb to sneaker
- Install google-cloud-sdk from deb
- Configure credentials using:
>    gcloud init --console-only

- Test using
>    gsutil ls gs://gbsc-gcp-project-mvp-received-from-bina


Transfering Data from S3 to Google Cloud Storage
------------------------------------------------

>    ssh ec2-user@10.1.1.80 '~/s3cmd/s3cmd-1.6.1/s3cmd ls s3://bina.va.uswest1.outputdata/${OUTPUTDIR}/${JOB_ID}' > filelist_${JOB_ID}
>
>    for F in $(cat filelist_${JOB_ID} | awk '{ print $4 }'); do
>        ssh ec2-user@10.1.1.80 "~/s3cmd/s3cmd-1.6.1/s3cmd get --no-progress ${F} -" | gsutil cp - gs://gbsc-gcp-project-mvp-received-from-bina/${OUTPUTDIR}/$(echo ${F} | cut -d"/" -f 5-)
>    done
