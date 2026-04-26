## Fungal RRN estimating pipeline

This pipeline is used for estimating fungal RRN (rDNA copy number) based on the sequencing depth of single-copy genes and rDNA genes. It can also be used for other multi-copy genes.

In addition, you will need to configure the Conda environment path to match your own setup and ensure the directory structure follows the example. Note that the pipeline currently needs to be run from a specific directory, which is not very flexible yet.

<br>

## Dependences
**Bowtie 2** version 2.5.3

**samtools** version 1.19.2

**R packages:** tidyverse and Biostrings

<br>

## Installation

To avoid software version conflicts, we recommend putting each dependency in its own conda environment.

First, make sure you have **Conda** or **Mamba** installed.

### Install the main dependencies
```bash
### Bowtie 2
conda create -n bowtie2
conda activate bowtie2
conda install -c conda-forge -c bioconda bowtie2
```

```bash
# switch back to base
conda activate base
```
```bash
### Samtools
conda create -n samtools
conda activate samtools
conda install -c conda-forge -c bioconda samtools
```bash

```bash
# switch back to base
conda activate base
```

```bash
### R environment
conda create -n R_env
conda activate R_env
conda install bioconda::bioconductor-biostrings
conda install conda-forge::r-tidyverse
```
<br>

### Prepare the input data

**Raw sequencing reads:** single-end or paired-end FASTQ files. Paired-end files should be named xxx_1.fastq / xxx_2.fastq (e.g. strain1_1.fastq and strain1_2.fastq), and single-end files can be named xxx.fastq.

**Single-copy gene sequences:** you can download them from MycoCosm or NCBI, or prepare your own based on genome annotation. File naming: `xxx_gene.fasta` (e.g. `strain1_RPB2.fasta`)

**rDNA gene sequences:** typically downloaded from NCBI or extracted from an assembled genome using ITSx or similar software.

<br>

### Project directory structure

The pipeline expects the following folder layout and the structure must be like this:

```bash
Fungal_RRN_test/
├── 0.script/    # all analysis scripts
└── 1.test/    # project directory
    ├── raw_data/    # your raw sequencing reads
    ├── single_copy_gene/    # single-copy gene sequences
    └── multi_copy_gene/    # rDNA gene sequences
```
<br>

### Running the pipeline

This pipeline is written in bash. Before you start, you need to manually update the Conda environment paths in the scripts (see line 31 in rCNV_1_BS_both.sh / rCNV_1_BS.sh, line 29 in rCNV_2_SBRD.sh, line 8 in rCNV_3_calcu.sh …).

Next, navigate to your project directory. Then run the appropriate command:

```bash
Single-end (unpaired) reads:
bash ../0.scripts/RCNV_pipeline.sh -t 30 # -t : number of threads
```
```bash
Paired-end reads:
bash ../0.scripts/RCNV_pipeline_both.sh -t 30 # -t : number of threads
```

The estimated results will be placed in the **CN_rlt** directory inside your project folder.
<br>

## An example

The example data was obtained from MycoCosm (sequencing project Neucr_trp3_1, https://mycocosm.jgi.doe.gov/Neucr_trp3_1/Neucr_trp3_1.home.html) and NCBI including genes and raw sequencing data. Because the raw sequencing data files are so large, we randomly extracted 20% of the reads for testing.


First, unzip the test file:
```bash
unzip calcu_fungal_rrn_test.zip
```

After unzipping, the test directory structure should look like this:
```bash
09.test/
├── raw_data/
├── single_copy_gene/
└── multi_copy_gene/
```

Then navigate to the project directory and run the pipeline:
```bash
cd 09.test/
bash ../0.scripts/RCNV_pipeline.sh -t 30
```

### Note on downloading data from MycoCosm:
Please note that MycoCosm enforces strict data usage policies. You must obtain the appropriate permissions in advance to avoid any potential issues. Raw sequencing data in MycoCosm is often stored on tape and not immediately available for download. You will need to register, request the data, and wait for it to be released.  
After that, download the project XML file and use **rCNV_0_export_downloadscript.R** to generate a bash download script.  
Finally, check the MD5 checksum to make sure the files are complete.


