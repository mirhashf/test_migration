echo "Getting needed packages"
sudo apt-get update
sudo apt-get -yq install  samtools s3cmd bwa pigz python-boto openjdk-6-jre-headless

echo "Setting up s3cmd"

echo "
[default]
access_key = ${aws_access_key}
acl_public = False
bucket_location = US
cloudfront_host = cloudfront.amazonaws.com
cloudfront_resource = /2008-06-30/distribution
default_mime_type = binary/octet-stream
delete_removed = False
dry_run = False
encoding = UTF-8
encrypt = False
force = False
get_continue = False
gpg_command = /usr/bin/gpg
gpg_decrypt = %(gpg_command)s -d --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
gpg_encrypt = %(gpg_command)s -c --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
gpg_passphrase = amiroo
guess_mime_type = True
host_base = s3.amazonaws.com
host_bucket = %(bucket)s.s3.amazonaws.com
human_readable_sizes = False
list_md5 = False
preserve_attrs = True
progress_meter = True
proxy_host =
proxy_port = 0
recursive = False
recv_chunk = 4096
secret_key = ${aws_access_secret_key}
send_chunk = 4096
simpledb_host = sdb.amazonaws.com
skip_existing = False
urlencoding_mode = normal
use_https = False
verbosity = WARNING
" > ~/.s3cfg


echo "Setting up boto"

echo "
[Credentials]
aws_access_key_id = ${aws_access_key}
aws_secret_access_key = ${aws_access_secret_key}
" > ~/.boto


echo "Get access"

sudo chmod -R 777 /mnt


echo "Create working dir"

mkdir ${prefix}

cd ${prefix}


echo "Get first read file"
time s3cmd get s3://${s3_in_bucket}/${lane}_1.fq.gz reads1.gz


echo "Get second read file"
time s3cmd get s3://${s3_in_bucket}/${lane}_2.fq.gz reads2.gz


echo "Get genome file"
time s3cmd get s3://seqalto/${bwa_genome_bundle}

echo "Get utils bundle"
time s3cmd get s3://seqalto/apps.tar.gz

echo "Unzip genome file"
time tar xvzf ${bwa_genome_bundle}


echo "Alignment"

time bwa aln -q 20 annotation/${bwa_genome} -1 reads1.gz -t ${threads} 1>${lane}-1.sai
time bwa aln -q 20 annotation/${bwa_genome} -1 reads2.gz -t ${threads} 1>${lane}-2.sai
time bwa sampe -P annotation/${bwa_genome} ${lane}-1.sai ${lane}-2.sai reads1.gz  reads2.gz -r "@RG\tID:${lane}\tLB:${library}\tSM:{sample}" | samtools view -bS /dev/stdin > ${lane}.bam

echo "Sort Bam"
#time samtools sort -m 6000000000 ${lane}.bam ${lane}_sorted && samtools index ${lane}_sorted.bam
time java -Xms${mem}g -Xmx${mem}g -jar apps/picard-tools-1.55/SortSam.jar I=${lane}.bam O=${lane}_sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT  TMP_DIR=.

samtools index ${lane}_sorted.bam

echo "Get multipart uploader script"
time wget https://raw.github.com/chapmanb/cloudbiolinux/master/utils/s3_multipart_upload.py && chmod +x s3_multipart_upload.py


echo "Upload sorted bam index file to s3:"
time ./s3_multipart_upload.py ${lane}_sorted.bam.bai ${s3_out_bucket}


echo "Upload sorted bam  file to s3:"
time ./s3_multipart_upload.py ${lane}_sorted.bam ${s3_out_bucket}


echo "Upload stdout and stderr to s3 and then Goodbye! :)"
time ./s3_multipart_upload.py ~/std.out ${s3_out_bucket} ${lane}.o
time ./s3_multipart_upload.py ~/std.err ${s3_out_bucket} ${lane}.e

sudo /sbin/halt -f
