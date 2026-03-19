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

geneNumbers=12

mainThreads=$Threads

subThreads=`expr $mainThreads / $geneNumbers`

if [ `expr $subThreads \* $geneNumbers` -lt $mainThreads ]
then
        subThreads=`expr $subThreads + 1`
else
        echo "pass this step"
fi

echo "$mainThreads"
echo "$subThreads"

##### prepare #####

bash ../0.scripts/rCNV_1_BS.sh -t ${mainThreads}

echo "step 1 done!"

bash ../0.scripts/rCNV_2_SBRD.sh -t ${subThreads}

echo "step 2 done!"

bash ../0.scripts/rCNV_3_calcu.sh

echo "step 3 done!"

rm ./raw_data/*.sam
rm ./raw_data/*.bam
rm ./bowtie_ref/*

echo "all done!"
