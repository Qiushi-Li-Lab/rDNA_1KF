#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

# back to base
conda activate base

seq_filename=`ls *.fna`
sam_id=${seq_filename%%_*}


if [ -d "./step_3" ]; then

	echo -e "\e[31m step_3 file is here! \e[0m"
	exit 1

else
	
	mkdir step_3

fi

cd step_3

# cp files
cp ../step_2/${sam_id}_predict/annotate_results/*.annotations.txt ./
cp ../step_2/${sam_id}_predict/KOG/KOG_diamond.txt ./
mv *.annotations.txt ${sam_id}_annotations.txt

cp ../step_2/${sam_id}_predict/annotate_results/*.cds-transcripts.fa ./
mv *.cds-transcripts.fa ${sam_id}_cds-transcripts.fa

# single_copy_gene index
Rscript ../../../funannotate_extract.R ${sam_id}_annotations.txt

# anno saved
mkdir anno
mv ${sam_id}_annotations.txt ./anno
mv KOG_diamond.txt ./anno


# single_gene_tab saved
mkdir single_gene_tab
mv *_tab.txt ./single_gene_tab

# single_gene_extra
mkdir extract_single_copy_gene
mv *.txt ./extract_single_copy_gene
mv *.bed ./extract_single_copy_gene

cd extract_single_copy_gene

bash ../../../../extract_single_copy_gene.sh

rm *txt
# rm *bed
# rm *tmp*fasta
# rm *seqkit.fai

cd ../

mkdir seqkit_single_copy_gene
mkdir single_copy_gene

cd extract_single_copy_gene

# cd-hit
conda activate cd-hit

for i in `ls *fasta`
do

        cd-hit -i ${i} -o ../seqkit_single_copy_gene/${i} -c 0.9 -n 4 -T 5 &

done

wait

rm *clstr

cd ../seqkit_single_copy_gene

# seqkit
conda activate seqkit

for i in `ls *fasta`
do

        seqkit sort --by-length --reverse ${i} | seqkit head -n 1 > ../single_copy_gene/${i} &

done

wait

cd ../

# rm ./single_copy_gene/*ELF1*

cp -R single_copy_gene ../../

echo -e "\e[31m step_3 is done! \e[0m"

# done.
