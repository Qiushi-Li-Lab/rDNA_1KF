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



# conda activate 
source /data2/software/mambaforge/etc/profile.d/conda.sh

# samtools
conda activate samtools

# raw data path
cd raw_data

# sam to bam
for i in `ls ../single_copy_gene/*fasta`
do

	FaID=${i##*/}
        RID=${FaID%.*}

	samtools view --threads ${Threads} -b ${RID}.sam > ${RID}.bam &

done

for i in `ls ../multi_copy_gene/*fasta`
do

        FaID=${i##*/}
        RID=${FaID%.*}
        
	samtools view --threads ${Threads} -b ${RID}.sam > ${RID}.bam &

done

wait

echo "all bam done!"

# bam sort
for i in `ls ../single_copy_gene/*fasta`
do

        FaID=${i##*/}
        RID=${FaID%.*}

        samtools sort --threads ${Threads} -T ${RID}_temp -o ${RID}_sorted.bam -O bam ${RID}.bam &

done

for i in `ls ../multi_copy_gene/*fasta`
do

        FaID=${i##*/}
        RID=${FaID%.*}

        samtools sort --threads ${Threads} -T ${RID}_temp -o ${RID}_sorted.bam -O bam ${RID}.bam &

done

wait

echo "all sort done!"


# back to project file
cd ..

# create depth_output
if [ -d "./depth_output" ]; then

        echo -e "\e[31m output files has been exist ! \e[0m"

else

        mkdir depth_output

fi

# raw data
cd raw_data

# depth
for i in `ls ../single_copy_gene/*fasta`
do
        FaID=${i##*/}
        RID=${FaID%.*}

        samtools depth -d 1000000 -q 20 --threads ${Threads} ${RID}_sorted.bam > ../depth_output/${RID}_depth_Q20.txt &
	samtools depth -d 1000000 --threads ${Threads} ${RID}_sorted.bam > ../depth_output/${RID}_depth_Q0.txt &
	# samtools depth --threads 10 ${RID}_sorted.bam > ../depth_output/${RID}_depth_Qnd.txt &

done

for i in `ls ../multi_copy_gene/*fasta`
do
        FaID=${i##*/}
        RID=${FaID%.*}

        samtools depth -d 1000000 -q 20 --threads ${Threads} ${RID}_sorted.bam > ../depth_output/${RID}_depth_Q20.txt &
	samtools depth -d 1000000 --threads ${Threads} ${RID}_sorted.bam > ../depth_output/${RID}_depth_Q0.txt &
        # samtools depth --threads 10 ${RID}_sorted.bam > ../depth_output/${RID}_depth_Qnd.txt &


done

wait

echo "all depth done!"

