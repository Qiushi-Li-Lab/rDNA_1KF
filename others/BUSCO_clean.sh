#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

# diamond 
conda activate busco_6.0.0


# busco
cd busco

# clean
cd ./output
rm -R ./logs
rm -R ./tmp

cd ../

# done.
