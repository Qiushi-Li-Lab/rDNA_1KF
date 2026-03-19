#! /usr/bin/bash



##### prepare #####


# conda activate
source /data2/software/mambaforge/etc/profile.d/conda.sh

# bowtie2
conda activate sff2fastq

for i in `ls ./raw_data/*.sff`
do
        FaID=${i##*/}
	RID=${FaID%.*}

	sff2fastq raw_data/${FaID} -o raw_data/${RID}.fastq

done 

echo "sff2fastq done!"
