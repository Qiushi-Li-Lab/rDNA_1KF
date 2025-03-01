#! /usr/bin/bash

# parameters
# t, threads

while getopts "t:" opt
do

        case $opt in

                t)
                        Threads=$OPTARG
                        ;;


                ?)

                        echo "threads should be tapped!"
                        exit 1
                        ;;

                esac

done


##### prepare #####


# conda activate
source /data2/software/mambaforge/etc/profile.d/conda.sh

# bowtie2
conda activate bowtie2

if [ -d "./bowtie_ref" ]; then

	echo -e "\e[31m reference files has been exist ! \e[0m"

else

	mkdir bowtie_ref

fi

for i in `ls ./single_copy_gene/*fasta`
do
        FaID=${i##*/}
	RID=${FaID%.*}

        bowtie2-build -f ./single_copy_gene/${FaID} ./bowtie_ref/${RID}_ref

done 

for i in `ls ./multi_copy_gene/*fasta`
do
        FaID=${i##*/}
        RID=${FaID%.*}

        bowtie2-build -f ./multi_copy_gene/${FaID} ./bowtie_ref/${RID}_ref

done

# change
cd raw_data

# raw data file list
# Rscript ../../0.scripts/bowtie_findex.R 

# bowtie2 align 
# Fq_ind=`sed -n "1p" fq_rlt.txt`


# gzip if gz file exist
gznum=`ls ./*gz | wc -l`

if [ ${gznum} -gt 0 ]; then

	gzip -d *gz

        echo "all gz has been done!"

else

        echo "no .gz file"
fi

# bowtie2 fq index
fqnum=`ls ./*fastq | wc -l`

if [ ${fqnum} -gt 0 ]; then

        Rscript ../../0.scripts/bowtie_findex.R

        echo "bowtie2 fastq index outputed! "

else

        echo "fail !!!"
fi

# bowtie2 align
Fq_ind=`sed -n "1p" fq_rlt.txt`


for i in `ls ../single_copy_gene/*fasta`
do

	FaID=${i##*/}
        RID=${FaID%.*}
	
	bowtie2 --very-sensitive-local -p ${Threads} -x ../bowtie_ref/${RID}_ref -U ${Fq_ind} -S ${RID}.sam

done

for i in `ls ../multi_copy_gene/*fasta`
do

        FaID=${i##*/}
        RID=${FaID%.*}
        bowtie2 --very-sensitive-local -p ${Threads} -x ../bowtie_ref/${RID}_ref -U ${Fq_ind} -S ${RID}.sam

done

echo "bowtie2 done!"
