#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

# funannotate
conda activate funannotate_1.8.15

seq_filename=`ls *.fna`
sam_id=${seq_filename%%_*}

if [ -d "./step_1" ]; then

	echo -e "\e[31m clean sort mask files has been exist ! \e[0m"

else

	mkdir step_1

fi

# clean
funannotate clean -i ${sam_id}_*.fna -o ./step_1/${sam_id}_clean.fna

# sort
funannotate sort -i ./step_1/${sam_id}_clean.fna -o ./step_1/${sam_id}_sort.fna -b ${sam_id}

# mask
funannotate mask -i ./step_1/${sam_id}_sort.fna -o ./step_1/${sam_id}_mask.fna --cpus 30 --debug

# done
echo "step_1 done!"
