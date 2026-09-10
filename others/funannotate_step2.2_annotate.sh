#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

seq_filename=`ls *.fna`
sam_id=${seq_filename%%_*}

# funannotate
conda activate funannotate_1.8.15

if [ -d "./step_2/${sam_id}_predict" ]; then

	echo -e "\e[31m predict files has been exist ! step_3 will go go go ! \e[0m"

else
	
	echo -e "\e[31m no predict results ! \e[0m"
	exit 1

fi

# annotation
funannotate annotate -i ./step_2/${sam_id}_predict --cpus 30

# diamond 
conda activate diamond

cd ./step_2/${sam_id}_predict

# KOG
mkdir KOG
cd KOG

# diamond KOGs
diamond blastp --db /data2/liqs/database/KOG/KOG_db/KOG.diamond.dmnd --query ../predict_results/*.proteins.fa --out KOG_diamond.txt --threads 30

# done.
