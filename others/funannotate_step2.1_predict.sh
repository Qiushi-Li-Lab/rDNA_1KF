#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

# funannotate
conda activate funannotate_1.8.15

# parameters
# S, spc

while getopts "S:" opt
do

	case $opt in

		S)
			spc=$OPTARG
			;;

		?)
			echo "spc should be tapped!"
			exit 1
			;;
	esac

done

seq_filename=`ls *.fna`
sam_id=${seq_filename%%_*}

if [ -d "./step_2" ]; then

	echo -e "\e[31m predict files has been exist ! \e[0m"

else

	mkdir step_2

fi

echo ${spc}

# predict
funannotate predict -i ./step_1/${sam_id}_mask.fna -o ./step_2/${sam_id}_predict -s "${spc}" --name ${sam_id}_ --cpus 30

# done...



