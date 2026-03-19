
# CNV calcu of Fungi
# by Qiushi-Li, IMCAS
# 2024.04.26

# the packages we needed in this pipeline
library(tidyverse)
library(Biostrings)


# gene index and path
single_gene_index <- list.files("./single_copy_gene/")
single_gene_index1 <- str_split(single_gene_index, pattern = "\\.", simplify = T)[, 1]
single_gene_path <- str_c("./single_copy_gene/", single_gene_index)

multi_gene_index <- list.files("./multi_copy_gene/")
multi_gene_index1 <- str_split(multi_gene_index, pattern = "\\.", simplify = T)[, 1]
multi_gene_path <- str_c("./multi_copy_gene/", multi_gene_index)


# depth table index and path
depth_table_ind <- list.files("./depth_output/")
depth_table_ind1 <- str_split(depth_table_ind, pattern = "\\.", simplify = T)[, 1]
depth_table_path <- str_c("./depth_output/", depth_table_ind)


# read the data
for(i in 1:length(single_gene_index1)) {
  
  assign(single_gene_index1[i], readDNAStringSet(single_gene_path[i])[[1]])
  
}

for(i in 1:length(multi_gene_index1)) {
  
  assign(multi_gene_index1[i], readDNAStringSet(multi_gene_path[i])[[1]])
  
}

for(i in 1:length(depth_table_ind1)) {
  
  assign(depth_table_ind1[i], read_table(depth_table_path[i], col_names = F))
  
}


# cbind
single_det <- str_split(single_gene_index1, pattern = "_", simplify = T)[, 2]
depth_table_s <- str_c(single_det, "_table")


for (i in 1:length(depth_table_s)) {
  
  X4 <- get(single_gene_index1[str_detect(single_gene_index1, single_det[i])])
  tmp_df <- get(depth_table_ind1[str_detect(depth_table_ind1, single_det[i])])
  
  min_length <- min(length(X4), length(tmp_df$X3))
  
  assign(depth_table_s[i], cbind(
    tmp_df[1:min_length, ],
    X4[1:min_length]
  ))
  
}

multi_det <- str_split(multi_gene_index1, pattern = "_", simplify = T)[, 2]
depth_table_m <- str_c(multi_det, "_table")

for (i in 1:length(depth_table_m)) {
  
  X4 <- get(multi_gene_index1[str_detect(multi_gene_index1, multi_det[i])])
  tmp_df <- get(depth_table_ind1[str_detect(depth_table_ind1, multi_det[i])])
  
  min_length <- min(length(X4), length(tmp_df$X3))
  
  assign(depth_table_m[i], cbind(
    tmp_df[1:min_length, ],
    X4[1:min_length]
  ))
  
}

# depth_table_s
# depth_table_m

## function to calculate GC within 100 bp of the target bp
sliding.bin.ave <- function(DNA_seq) {
  # this function calculates the GC % of each target bin, with bin size 100, and target bp centered
  n <- 50
  total <- length(DNA_seq)
  bins <- seq(from = 1, to = (total))
  result <- rep(NA, times = length(DNA_seq))
  thisGC <- as.integer(DNA_seq %in% c("G", "C"))
  for(i in 51:length(bins)){
    result[i] <- sum(thisGC[(bins[i] - n):(bins[i] + n)])
  }
  return(result)
}

for (i in 1:length(depth_table_s)) {
  
  X5 <- sliding.bin.ave(DNA_seq = get(depth_table_s[i])[, 4])
  
  assign(depth_table_s[i],
         cbind(
           get(depth_table_s[i]),
           X5
         )
  )
  
}

for (i in 1:length(depth_table_m)) {
  
  X5 <- sliding.bin.ave(DNA_seq = get(depth_table_m[i])[, 4])
  
  assign(depth_table_m[i],
         cbind(
           get(depth_table_m[i]),
           X5
         )
  )
  
}

# clean it up #1 ---------------
# delete row if there's an ambiguity code 
# depth_table_s
# depth_table_m

BP_list<- c("A", "C", "T", "G")

# single-copy
depth_table_s1 <- str_c(depth_table_s, "_1")
for(i in 1:length(depth_table_s1)) {
  
  tmp_df <- get(depth_table_s[i])
  
  assign(depth_table_s1[i],
         tmp_df[tmp_df$X4 %in% BP_list, ]
  )
  
}

# multi-copy
depth_table_m1 <- str_c(depth_table_m, "_1")
for(i in 1:length(depth_table_m1)) {
  
  tmp_df <- get(depth_table_m[i])
  
  assign(depth_table_m1[i],
         tmp_df[tmp_df$X4 %in% BP_list, ]
  )
  
}


# clean it up #2 --------------------
# delete row if there are depth values of zero, or NA's
# single-copy
depth_table_s2 <- str_c(depth_table_s, "_2")
for(i in 1:length(depth_table_s2)) {
  
  tmp_df <- get(depth_table_s1[i])
  
  assign(depth_table_s2[i],
         na.omit(tmp_df[!(tmp_df$X3==0), ])
  )
  
}

# multi-copy
depth_table_m2 <- str_c(depth_table_m, "_2")
for(i in 1:length(depth_table_m2)) {
  
  tmp_df <- get(depth_table_m1[i])
  
  assign(depth_table_m2[i],
         na.omit(tmp_df[!(tmp_df$X3==0), ])
  )
  
}

## change the names so that they match and can be merged ------------------
NEW.names <- c("sp.", "bp.pos", "long.depth","bp", "gc.bin")
# single-copy
for(i in 1:length(depth_table_s2)) {
  
  tmp_df <- get(depth_table_s2[i])
  colnames(tmp_df) <- NEW.names
  assign(depth_table_s2[i], tmp_df)
  
}

# multi-copy
for(i in 1:length(depth_table_m2)) {
  
  tmp_df <- get(depth_table_m2[i])
  colnames(tmp_df) <- NEW.names
  assign(depth_table_m2[i], tmp_df)
  
}

## chop first and last 50 pb to use for bin ave. calculation -----------------------
# single-copy
chop_table_s <- str_c(depth_table_s, "_chop")
for(i in 1:length(chop_table_s)) {
  
  tmp_df <- get(depth_table_s2[i])
  assign(chop_table_s[i], tmp_df[50:(nrow(tmp_df) - 50), ])
  
}

# multi-copy
chop_table_m <- str_c(depth_table_m, "_chop")
for(i in 1:length(chop_table_m)) {
  
  tmp_df <- get(depth_table_m2[i])
  assign(chop_table_m[i], tmp_df[50:(nrow(tmp_df) - 50), ])
  
}


# STOP HERE: check to see how good the data looks ----------------
# SCG's
pdf("single_chop_barplot.pdf")
for(i in 1:length(chop_table_s)) {
  
  tmp_df <- get(chop_table_s[i])
  barplot(tmp_df$long.depth, main = chop_table_s[i])
  
}
dev.off()

pdf("multi_chop_barplot.pdf")
for(i in 1:length(chop_table_m)) {
  
  tmp_df <- get(chop_table_m[i])
  barplot(tmp_df$long.depth, main = chop_table_m[i])
  
}
dev.off()

## are any of the "SC" genes likely to be "MC"? If they're more than 1 sd outside the median, don't include downstream.
# get means
# single_det
# multi_det
long_means_single <- str_c("long_mean_", single_det)
for(i in 1:length(long_means_single)) {
  
  tmp_df <- get(chop_table_s[i])
  assign(long_means_single[i], mean(tmp_df$long.depth))
  
}

# make a list of means
long_means <- vector()
for(i in 1:length(long_means_single)) {
  
  tmp_value <- get(long_means_single[i])
  long_means <- c(long_means, tmp_value)
  
}


# name the genes
names(long_means) <- long_means_single

# get cutoffs 
medianSCG <- median(long_means)
sdSCG <- sd(long_means)
upper <- medianSCG + sdSCG
lower <- medianSCG - sdSCG

# don't include genes in the next step that return 'FALSE', (above abline in graph). -------------
only_genes_in_range <- long_means < upper & long_means > lower
only_genes_in_range
# GAPDH and MLS

# gene in range
pdf("gene_in_range.pdf")
barplot(long_means, las = 2)
abline(a = upper, b = 0)
abline(a = lower, b = 0)
dev.off()


#combine all sc genes, # hash out genes above or below ab-line in graph
chop_table_s_in_range <- chop_table_s[only_genes_in_range]
all_SC_gene_tables <- data.frame()
for(i in 1:length(chop_table_s_in_range)) {
  
  tmp_df <- get(chop_table_s_in_range[i])
  all_SC_gene_tables <- rbind(all_SC_gene_tables, tmp_df)
  
}


# RD median for all bins --------------
all_SC_bins <- median(all_SC_gene_tables[, 3])

# get depth medians for each unique bin -------------
depth.ea.bin.all <- aggregate(data = all_SC_gene_tables,
                              all_SC_gene_tables[, 3] ~ all_SC_gene_tables[, 5],
                              FUN = median) 


# Function to adjust depth by GC bins --------------
cor.depth <- function(depth.df) {
  
  result  <- seq(from=1, to=(nrow(depth.df))) 
  x <- all_SC_bins # median RD over all single copy genes
  
  for(i in 1:nrow(depth.df)) {
    
    for (j in 1:nrow(depth.ea.bin.all)) {
      
      if (depth.ea.bin.all[j, 1] == depth.df[i, 5]) {
        result[i] <- depth.df[i, 3] * (x / (depth.ea.bin.all[j, 2]))
        
      }
      
    }
    
  }
  return(result)
}


# run function on each single copy region over the relevent bp's (100 chopped from either side), and get mean --------------
# hash out genes above or below ab-line in graph
chop_table_s_in_range
adj_gene <- str_c("adj_", str_split(chop_table_s_in_range, pattern = "_", simplify = T)[, 1])
for(i in 1:length(adj_gene)) {
  
  tmp_df <- get(chop_table_s_in_range[i])
  assign(adj_gene[i], mean(cor.depth(depth.df = tmp_df)))
  
}


# make a list of adjusted means ----------------------
# hash out genes above or below ab-line in graph
adjusted_means <- vector()
for(i in 1:length(adj_gene)){
  
  tmp_value <- get(adj_gene[i])
  adjusted_means <- c(adjusted_means, tmp_value)
  
}

gene_names <- str_c(str_split(chop_table_s_in_range, pattern = "_", simplify = T)[, 1], "=")
# gene_names


## to adjust multi-copy region of interest (ITS). ----------------------
cor.depth.multi_ITS <- function(depth.df) {
  # this function takes a dataframe (depth.df) with read depth at each bp, and the average GC content (%) of the 100 bps surrounding the bp of interest 
  result <- seq(from = 1, to = (nrow(depth.df)))
  # ave RD over ITS unadjusted bins
  x <- mean(ITS_table_chop$long.depth) 
  # median read depth for each unique GC bin (percent GC over 100 bp, centered on the bp of interest)
  depth.ea.bin.ITS <- aggregate(data = ITS_table_chop, ITS_table_chop[,3] ~ ITS_table_chop[,5], FUN = median)
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.ITS)){
      #for each unique GC % median
      if (depth.ea.bin.ITS[j, 1] == depth.df[i, 5]) {
        #for each pb position, multiply the depth by (the average over all positions devided by the median for that unique GC %)
        result[i] <- depth.df[i, 3] * (x / (depth.ea.bin.ITS[j, 2]))
      }
    }
  }
  return(result)
}

corrected.ITS.dep <- cor.depth.multi_ITS(depth.df = ITS_table_chop)
mean.ITS.corrected <- mean(corrected.ITS.dep)
mean.ITS.corrected


## to adjust multi-copy region of interest (LSU). -----------------------
cor.depth.multi_LSU <- function(depth.df){
  #this function takes a dataframe (depth.df) with read depth at each bp, and the average GC content (%) of the 100 bps surrounding the bp of interest
  result <- seq(from = 1, to = (nrow(depth.df)))
  #ave RD over LSU unadjusted bins
  x <- mean(LSU_table_chop$long.depth)
  #median read depth for each unique GC bin (percent GC over 100 bp, centered on the bp of interest)
  depth.ea.bin.LSU <- aggregate(data = LSU_table_chop, LSU_table_chop[, 3] ~ LSU_table_chop[, 5], FUN = median) 
  for(i in 1:nrow(depth.df)){
    for (j in 1:nrow(depth.ea.bin.LSU)){
      #for each unique GC % median
      if (depth.ea.bin.LSU[j, 1] == depth.df[i, 5]) {
        #for each pb position, multiply the depth by (the average over all positions devided by the median for that unique GC %)
        result[i] <- depth.df[i, 3] * (x / (depth.ea.bin.LSU[j, 2]))
      }
    }
  }
  return(result)
}

corrected.LSU.dep <- cor.depth.multi_LSU(depth.df = LSU_table_chop)
mean.LSU.corrected <- mean(corrected.LSU.dep)
mean.LSU.corrected


# check that you trust the data ------------------
pdf("check_trust.pdf")
barplot(ITS_table_chop$long.depth, main = "ITS before correction")
barplot(corrected.ITS.dep, main = "ITS after correction")

barplot(LSU_table_chop$long.depth, main = "LSU before correction")
barplot(corrected.LSU.dep, main = "LSU after correction")
dev.off()

# Average depth ---------------
SCG_ave <- mean(adjusted_means)
MC_ave <- (mean.LSU.corrected + mean.ITS.corrected) / 2

# Estimated Copy Number of rDNA cassette 
CN_est_ITS <- round(mean.ITS.corrected / SCG_ave)
CN_est <- round(MC_ave / SCG_ave)


# percent difference between ITS and LSU depth 
per_diff <- round(abs(mean.LSU.corrected - mean.ITS.corrected) / ((mean.LSU.corrected + mean.ITS.corrected) / 2)*100, digits = 3)

# generate output -----------
# individual totals
df <- data.frame(cbind(gene_names, adjusted_means))
write.table(df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

names.2 <- "ITS="
ITS.df <- data.frame(cbind(names.2, mean.ITS.corrected))
write.table(ITS.df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

names.3 <- "LSU="
LSU.df <- data.frame(cbind(names.3, mean.LSU.corrected))
write.table(LSU.df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#summary totals 
names.4 <- "SCG depth average="
SCG.df <- data.frame(cbind(names.4, round(SCG_ave)))
write.table(SCG.df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

names.5 <- "estimated rDNA CN (ITS only)="
CN.df.ITS <- data.frame(cbind(names.5, CN_est_ITS))
write.table(CN.df.ITS, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

names.6 <- "estimated rDNA CN (ITS / LSU)="
CN.df <- data.frame(cbind(names.6, CN_est))
write.table(CN.df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

names.7 <- "% diff ITS / LSU="
CN.df <- data.frame(cbind(names.7, per_diff))
write.table(CN.df, file = "./CN_rlt/fungi_depth_totals.txt",
            append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# done.