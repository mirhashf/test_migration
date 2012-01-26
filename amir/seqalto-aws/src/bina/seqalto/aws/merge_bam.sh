date
echo "Getting needed packages"
sudo apt-get -yq install samtools s3cmd pigz python-boto

date
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

date
echo "Setting up boto"

echo "
[Credentials]
aws_access_key_id = ${aws_access_key}
aws_secret_access_key = ${aws_access_secret_key}
" > ~/.boto

date
echo "Get multipart uploader script"
wget https://raw.github.com/chapmanb/cloudbiolinux/master/utils/s3_multipart_upload.py && chmod +x s3_multipart_upload.py


date
echo "Get access"

sudo chmod -R 777 /mnt

date
echo "Create working dir"

mkdir ${prefix}

cd ${prefix}


date
echo "Get BAM files"
c=0
for i in $$*
do 
./s3-mp-download.py -np 8 $$i lane_$$c.bam -f
let c=c+1  
done

samtools merge -f merged.bam lane_*.bam 1>~std.out 2>std.err
samtools index merged.bam 1>std.out 2>std.err

date
echo "Upload merged bam file to s3:"
./s3_multipart_upload.py merged.bam ${s3_out_bucket} merged.bam

echo "Upload merged bam index file to s3:"
./s3_multipart_upload.py merged.bam.bai ${s3_out_bucket} merged.bam.bai


sudo /sbin/halt -f
