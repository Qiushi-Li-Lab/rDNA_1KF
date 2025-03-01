#! /usr/bin/bash


# login
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=liqiushi192@mails.ucas.ac.cn' --data-urlencode 'password=asdf0147' -c cookies > /dev/null


# download curl
# curl 'https://genome-downloads.jgi.doe.gov

curl 'https://genome-downloads.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Fusredo2/download/_JAMO/5f8590e447675a20c84f99ed/pbio-2300.21224.ccs.fastq.gz' -b cookies > pbio-2300.21224.ccs.fastq.gz


md5sum * > md5summary.txt

echo "done!"



