
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
R1_rlt <- str_c(R1_list, collapse = ",")
R2_rlt <- str_c(R2_list, collapse = ",")

# write
write.table(R1_rlt, "R1_rlt.txt", append = F, quote = F, row.names = F, col.names = F)
write.table(R2_rlt, "R2_rlt.txt", append = F, quote = F, row.names = F, col.names = F)

# done.
