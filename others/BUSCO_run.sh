#! /usr/bin/bash

# conda
source /data2/software/mambaforge/etc/profile.d/conda.sh

# diamond 
conda activate busco_6.0.0


# busco
mkdir busco
cd busco

# genome quality
busco -i ../all_file/*_AssemblyScaffolds.fasta -o output -m genome -c 20 -l /data1/busco_db/lineages/fungi_odb12

# done.
