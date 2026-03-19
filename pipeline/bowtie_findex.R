
# bowtie2 fileindex
# by qiushi-Li, IMCAS

# packages
library(stringr)

# all fastq
fq_list <- list.files(pattern = ".fastq")

# R1 fastq
R1_list <- fq_list[str_detect(fq_list, pattern = "_1.fastq")]

# R2 fastq
R2_list <- fq_list[str_detect(fq_list, pattern = "_2.fastq")]

# rlt
fq_rlt_tmp <- fq_list[!fq_list %in% c(R1_list, R2_list)]
fq_rlt <- str_c(fq_rlt_tmp, collapse = ",")

write.table(fq_rlt, "fq_rlt.txt", append = F, quote = F, row.names = F, col.names = F)


# done.
