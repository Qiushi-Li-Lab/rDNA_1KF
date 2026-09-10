#!/bin/bash

# mkfifo lock file

mkfifo ./parallel_lock_busco
exec 3<>./parallel_lock_busco
rm -f ./parallel_lock_busco

# parallel commands number = 3
for ((i=0; i<3; i++)); do
    echo >&3
done

# your command
for i in `cat ./proj_list.txt`
do
    # check key
    read -u 3
    {
        
	FaID=${i}
        cd ./rDNA_genome/${FaID}
	bash ../../BUSCO_run.sh
	# bash ../../BUSCO_clean.sh
	cd ../../../
	echo "${FaID} done!"

        # retrun key
        echo >&3
    } &  # 
done

# 
wait

# 
exec 3>&-
