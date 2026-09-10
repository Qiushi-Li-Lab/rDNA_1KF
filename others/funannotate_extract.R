

##### for export single copy gene names --------------------
# Qiushi-Li, IM-CAS
# 2025.01.10


# packages
library(tidyverse)

# filename
funannotate_filename <- commandArgs(trailingOnly = T)


# read the data
funannotate_df <- read_delim(funannotate_filename, delim = "\t")
# funannotate_df <- read_delim("./anno/AMF35_annotations.txt", delim = "\t")

# strain_id
strain_id <- str_split_1(funannotate_filename, "\\_")[1]
# strain_id <- "AMF35"

# sub data
funannotate_df1 <- funannotate_df %>% 
  select(GeneID, TranscriptID, Name, Product, EC_number, BUSCO, PFAM, `GO Terms`, CAZyme)

# view(funannotate_df)

# KOG
#KOG_rlt <- read.table("../KOG/KOG_diamond.txt", header = F)
KOG_anno <- read.table("/data2/liqs/database/KOG/raw_data/kog.parsed.tab", header = F, sep = "\t")

KOG_rlt <- read.table("KOG_diamond.txt", header = F)
# KOG_anno <- read.table("kog.parsed.tab", header = F, sep = "\t")


# data process --------------------
# re-col name
colnames(KOG_rlt) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(KOG_anno) <- c("sseqid", "KOG_group", "KOG_Description")



KOG_rlt1 <- KOG_rlt %>% 
  left_join(KOG_anno, by = "sseqid")
# KOG_rlt1


# colnames(KOG_rlt1) <- c("TranscriptID", "KOG")
# funannotate_df1 <- funannotate_df1 %>% left_join(KOG_rlt1, by = c("TranscriptID" = "qseqid"))

# RPB1 ---------------------
FUN_RPB1 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0216")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(bitscore == max(bitscore))
FUN_RPB1

RPB1_filename <- str_c(strain_id, "_RPB1_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_RPB1, RPB1_filename, col.names = T, row.names = F, quote = F, sep = "\t")


RPB1_index <- FUN_RPB1 %>% select(qseqid) %>% distinct_all()
RPB1_indexname <- str_c(strain_id, "_RPB1", ".txt", sep = "")

# gene index
write.table(RPB1_index, RPB1_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# RPB2 ---------------------
FUN_RPB2 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0214")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(bitscore == max(bitscore))
FUN_RPB2

RPB2_filename <- str_c(strain_id, "_RPB2_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_RPB2, RPB2_filename, col.names = T, row.names = F, quote = F, sep = "\t")


RPB2_index <- FUN_RPB2 %>% select(qseqid) %>% distinct_all()
RPB2_indexname <- str_c(strain_id, "_RPB2", ".txt", sep = "")

# gene index
write.table(RPB2_index, RPB2_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# GAPDH ---------------------
FUN_GAPDH <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0657")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF00044")) %>% 
  filter(bitscore == max(bitscore))
FUN_GAPDH

GAPDH_filename <- str_c(strain_id, "_GAPDH_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_GAPDH, GAPDH_filename, col.names = T, row.names = F, quote = F, sep = "\t")


GAPDH_index <- FUN_GAPDH %>% select(qseqid) %>% distinct_all()
GAPDH_indexname <- str_c(strain_id, "_GAPDH", ".txt", sep = "")

# gene index
write.table(GAPDH_index, GAPDH_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# ELF1 ---------------------
FUN_ELF1 <- KOG_rlt1 %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF05129")) %>% 
  filter(bitscore == max(bitscore))
FUN_ELF1

ELF1_filename <- str_c(strain_id, "_ELF1_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_ELF1, ELF1_filename, col.names = T, row.names = F, quote = F, sep = "\t")


ELF1_index <- FUN_ELF1 %>% select(qseqid) %>% distinct_all()
ELF1_indexname <- str_c(strain_id, "_ELF1", ".txt", sep = "")

# gene index
write.table(ELF1_index, ELF1_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# GH63 ---------------------
FUN_GH63 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG2161")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF03200")) %>% 
  filter(bitscore == max(bitscore))
FUN_GH63

GH63_filename <- str_c(strain_id, "_GH63_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_GH63, GH63_filename, col.names = T, row.names = F, quote = F, sep = "\t")


GH63_index <- FUN_GH63 %>% select(qseqid) %>% distinct_all()
GH63_indexname <- str_c(strain_id, "_GH63", ".txt", sep = "")

# gene index
write.table(GH63_index, GH63_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# MCM7 ---------------------
FUN_MCM7 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0482")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF00493")) %>% 
  filter(bitscore == max(bitscore))
FUN_MCM7

MCM7_filename <- str_c(strain_id, "_MCM7_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_MCM7, MCM7_filename, col.names = T, row.names = F, quote = F, sep = "\t")


MCM7_index <- FUN_MCM7 %>% select(qseqid) %>% distinct_all()
MCM7_indexname <- str_c(strain_id, "_MCM7", ".txt", sep = "")

# gene index
write.table(MCM7_index, MCM7_indexname, col.names = F, row.names = F, quote = F, sep = "\t")

# G6PDH ---------------------
FUN_G6PDH <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0563")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF00479")) %>% 
  filter(bitscore == max(bitscore))
FUN_G6PDH

G6PDH_filename <- str_c(strain_id, "_G6PDH_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_G6PDH, G6PDH_filename, col.names = T, row.names = F, quote = F, sep = "\t")


G6PDH_index <- FUN_G6PDH %>% select(qseqid) %>% distinct_all()
G6PDH_indexname <- str_c(strain_id, "_G6PDH", ".txt", sep = "")

# gene index
write.table(G6PDH_index, G6PDH_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# MLS ---------------------
FUN_MLS <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG1261")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF01274")) %>% 
  filter(bitscore == max(bitscore))
FUN_MLS

MLS_filename <- str_c(strain_id, "_MLS_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_MLS, MLS_filename, col.names = T, row.names = F, quote = F, sep = "\t")


MLS_index <- FUN_MLS %>% select(qseqid) %>% distinct_all()
MLS_indexname <- str_c(strain_id, "_MLS", ".txt", sep = "")

# gene index
write.table(MLS_index, MLS_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# LYS2 ---------------------
FUN_LYS2 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG1178")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF00501")) %>% 
  filter(bitscore == max(bitscore))
FUN_LYS2

LYS2_filename <- str_c(strain_id, "_LYS2_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_LYS2, LYS2_filename, col.names = T, row.names = F, quote = F, sep = "\t")


LYS2_index <- FUN_LYS2 %>% select(qseqid) %>% distinct_all()
LYS2_indexname <- str_c(strain_id, "_LYS2", ".txt", sep = "")

# gene index
write.table(LYS2_index, LYS2_indexname, col.names = F, row.names = F, quote = F, sep = "\t")


# TOP2 ---------------------
FUN_TOP2 <- KOG_rlt1 %>% filter(str_detect(KOG_Description, "KOG0355")) %>% 
  left_join(funannotate_df1, by = c("qseqid" = "TranscriptID")) %>% 
  filter(str_detect(PFAM, "PF00204")) %>% 
  filter(bitscore == max(bitscore))
FUN_TOP2

TOP2_filename <- str_c(strain_id, "_TOP2_tab", ".txt", sep = "")

# gene table filtered
write.table(FUN_TOP2, TOP2_filename, col.names = T, row.names = F, quote = F, sep = "\t")


TOP2_index <- FUN_TOP2 %>% select(qseqid) %>% distinct_all()
TOP2_indexname <- str_c(strain_id, "_TOP2", ".txt", sep = "")

# gene index
write.table(TOP2_index, TOP2_indexname, col.names = F, row.names = F, quote = F, sep = "\t")

# done .......
