#! /usr/bin/bash

mkdir CN_rlt

mkdir CN_rlt_tmp

# conda activate
source /data2/software/mambaforge/etc/profile.d/conda.sh

# return to base
conda activate base

# Q0
mv ./depth_output/*txt ./CN_rlt_tmp

cp ./CN_rlt_tmp/*Q0.txt ./depth_output

Rscript ../0.scripts/calcu_CN.R

mv ./CN_rlt/fungi_depth_totals.txt ./CN_rlt/Q0_fungi_depth_totals.txt

rm ./depth_output/*txt

echo "Q0 done.!"


# Q20
cp ./CN_rlt_tmp/*Q20.txt ./depth_output

Rscript ../0.scripts/calcu_CN.R

mv ./CN_rlt/fungi_depth_totals.txt ./CN_rlt/Q20_fungi_depth_totals.txt

rm ./depth_output/*txt

echo "Q20 done.!"

