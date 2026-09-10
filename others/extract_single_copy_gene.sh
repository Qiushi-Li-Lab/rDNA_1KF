#! /bin/bash

# conda activate
source /data2/software/mambaforge/etc/profile.d/conda.sh

# fastqc
conda activate seqkit


for i in `ls *txt`
do
	
	ID=${i%.*}

	echo ${ID}

	ID1=${ID#*_}
	echo ${ID1}

	ID2=${ID%_*}

	echo ${ID2}

	seqkit grep -f *${ID1}*.txt ../*.fa -o ${ID2}_${ID1}.fasta
	
	# seqkit subseq --bed ${ID2}_${ID1}.bed -o ${ID2}_${ID1}.fasta ${ID2}_${ID1}_tmp.fasta


done

echo -e "\e[32m all process  has been done! please check!!! \e[0m"
