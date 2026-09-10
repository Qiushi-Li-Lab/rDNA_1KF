

###### step 5 mock commu ######
# by Qiushi-Li, IM-CAS


#####  packages we needed.
library(tidyverse) # R4DS

library(reshape2)

# read & writexl
library(readxl)
library(writexl)

# multi session calcu
library(furrr) 
library(future)
plan(multisession)

# numeric ecology
library(vegan)
library(compositions)


# radar plot
library(ggradar)
library(gghalves)

# net plot
# library(igraph)
# library(ggraph)
# library(tidygraph)

library(ggrepel)

# merge plot
library(patchwork)


# dPCR results
dPCR_rlts <- read_excel("./1.data/dPCR_rlts.xlsx", sheet = 1)
dPCR_rlts <- dPCR_rlts %>% mutate(
  strain_id = strain_ID,
  Taxa = str_split_i(Taxa, pattern = " ", 1)
) %>% dplyr::select(strain_id, Taxa)


###### part1 filtered ASV ref db ######
# ASVs table
all_table <- read_tsv("./1.data/fungal_mock_commu/0.all_ASVs/all_table.tsv", skip = 1)
colnames(all_table)[1] <- "OTU_ID"

# ASVs ref seqs 
all_ASVs_ref <- readDNAStringSet("./1.data/fungal_mock_commu/0.all_ASVs/all_ref_seqs.fasta")

# ASVs ref taxa
all_ASVs_taxa <- read_tsv("./1.data/fungal_mock_commu/0.all_ASVs/all_ref_seqs_taxa.tsv")
colnames(all_ASVs_taxa)[1] <- "OTU_ID"


all_ASVs_taxa1 <- all_ASVs_taxa %>%
  separate_wider_delim(
    Taxon, names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    delim = ";", too_few = "debug", ) %>%
  mutate(
    Phylum = replace_na(Phylum, replace = "p__unknown"),
    Class = replace_na(Class, replace = "c__unknown"),
    Order = replace_na(Order, replace = "o__unknown"),
    Family = replace_na(Family, replace = "f__unknown"),
    Genus = replace_na(Genus, replace = "g__unknown"),
    Species = replace_na(Species, replace = "s__unknown"),
  ) %>% select(-Taxon_ok, -Taxon_pieces, -Taxon_remainder)


all_ASVs_total_abd <-
  all_table %>%
  rowwise() %>%
  mutate(total_abd = sum(c_across(-OTU_ID))) %>%
  mutate(total_count = specnumber(c_across(-OTU_ID))) %>% 
  select(OTU_ID, total_abd, total_count) %>% arrange(desc(total_abd)) %>%
  left_join(all_ASVs_taxa1, by = "OTU_ID") %>%
  ungroup()
all_ASVs_total_abd

all_ASVs_total_abd$Kingdom %>% table()
# k__Alveolata         k__Fungi       k__Metazoa k__Viridiplantae       Unassigned 
# 1             4364               11              340             2771 

all_ASVs_total_abd_filtered <- all_ASVs_total_abd %>% filter(Kingdom == "k__Fungi")


# S01
all_ASVs_total_abd_filtered

ASVs_S01_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Pleurotus_ostreatus"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Pleurotus") %>% count(Species)

# S02
ASVs_S02_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Saccharomyces_cerevisiae"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Saccharomyces") %>% count(Species)

# S03
ASVs_S03_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Neurospora_crassa"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Neurospora") %>% count(Species)

# S04
ASVs_S04_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Podila_humilis"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Podila") %>% count(Species)

# S05
ASVs_S05_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Fusarium_proliferatum"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Fusarium") %>% count(Species)

# S07
ASVs_S07_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Gongronella_butleri"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Gongronella") %>% count(Species)
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Gongronella" & Species == "s__unknown")

# S08
ASVs_S08_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Rhizophagus_irregularis"))
ASVs_S08_yes_S08 <- ASVs_S08_yes %>% filter(str_detect(Species, pattern = "S08"))
ASVs_S08_yes_spc <- ASVs_S08_yes %>% filter(!OTU_ID %in% ASVs_S08_yes_S08$OTU_ID)

ASVs_S08_maybe <- all_ASVs_total_abd_filtered %>% filter(Genus == "g__Rhizophagus" & Species == "s__unknown")

ASVs_S08_all <- all_ASVs_total_abd_filtered %>% filter(Genus == "g__Rhizophagus")

# all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Rhizophagus_irregularis")) %>% count(Species)
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Rhizophagus") %>% count(Species)

# S10
ASVs_S10_yes <- all_ASVs_total_abd_filtered %>% filter(str_detect(Species, pattern = "Agaricus_bisporus"))
# all_ASVs_total_abd_filtered %>% filter(Genus == "g__Agaricus") %>% count(Species)

# usearch_global db
all_ASVs_ref

# S01
S01_ref <- all_ASVs_ref[ASVs_S01_yes$OTU_ID]
names(S01_ref) <- str_c("S01_", 1:length(S01_ref))
# S02
S02_ref <- all_ASVs_ref[ASVs_S02_yes$OTU_ID]
names(S02_ref) <- str_c("S02_", 1:length(S02_ref))
# S03
S03_ref <- all_ASVs_ref[ASVs_S03_yes$OTU_ID]
names(S03_ref) <- str_c("S03_", 1:length(S03_ref))
# S04
S04_ref <- all_ASVs_ref[ASVs_S04_yes$OTU_ID]
names(S04_ref) <- str_c("S04_", 1:length(S04_ref))
# S05
S05_ref <- all_ASVs_ref[ASVs_S05_yes$OTU_ID]
names(S05_ref) <- str_c("S05_", 1:length(S05_ref))
# S07
S07_ref <- all_ASVs_ref[ASVs_S07_yes$OTU_ID]
names(S07_ref) <- str_c("S07_", 1:length(S07_ref))
# S08
S08_ref_S08 <- all_ASVs_ref[ASVs_S08_yes_S08$OTU_ID]
names(S08_ref_S08) <- str_c("S08_", 1:length(S08_ref_S08))

S08_ref_spc <- all_ASVs_ref[ASVs_S08_yes_spc$OTU_ID]
names(S08_ref_spc) <- str_c("S08_spc_", 1:length(S08_ref_spc))

S08_ref_maybe <- all_ASVs_ref[ASVs_S08_maybe$OTU_ID]
names(S08_ref_maybe) <- str_c("S08_maybe_", 1:length(S08_ref_maybe))

# writeXStringSet(S08_ref_maybe, "S08_unknow.fasta")

# S10
S10_ref <- all_ASVs_ref[ASVs_S10_yes$OTU_ID]
names(S10_ref) <- str_c("S10_", 1:length(S10_ref))

mock_commu_ASVs_db <- c(
  S01_ref, S02_ref, S03_ref, S04_ref, S05_ref, S07_ref,
  S08_ref_S08, S08_ref_spc, S08_ref_maybe,
  S10_ref)

# ref_ASVs db
mock_commu_ASVs_db
# writeXStringSet(mock_commu_ASVs_db, "./1.data/mock_commu_ASVs_db_1.fasta")


# diversity of rDNA
# S01
ASVs_S01_abd1000 <- ASVs_S01_yes %>% filter(total_abd > 1000)
S01_div <- all_ASVs_ref[ASVs_S01_abd1000$OTU_ID]
names(S01_div) <- str_c("S01_div_", 1:length(S01_div))

ASVs_S01_abd1000 <- ASVs_S01_abd1000 %>% 
  mutate(Seq_ID = str_c("S01_div_", 1:length(S01_div)))


# S02
ASVs_S02_abd1000 <- ASVs_S02_yes %>% filter(total_abd > 1000)
S02_div <- all_ASVs_ref[ASVs_S02_abd1000$OTU_ID]
names(S02_div) <- str_c("S02_div_", 1:length(S02_div))

ASVs_S02_abd1000 <- ASVs_S02_abd1000 %>% 
  mutate(Seq_ID = str_c("S02_div_", 1:length(S02_div)))


# S03
ASVs_S03_abd1000 <- ASVs_S03_yes %>% filter(total_abd > 1000)
S03_div <- all_ASVs_ref[ASVs_S03_abd1000$OTU_ID]
names(S03_div) <- str_c("S03_div_", 1:length(S03_div))

ASVs_S03_abd1000 <- ASVs_S03_abd1000 %>% 
  mutate(Seq_ID = str_c("S03_div_", 1:length(S03_div)))

# S04
ASVs_S04_abd1000 <- ASVs_S04_yes %>% filter(total_abd > 1000)
S04_div <- all_ASVs_ref[ASVs_S04_abd1000$OTU_ID]
names(S04_div) <- str_c("S04_div_", 1:length(S04_div))

ASVs_S04_abd1000 <- ASVs_S04_abd1000 %>% 
  mutate(Seq_ID = str_c("S04_div_", 1:length(S04_div)))

# S05
ASVs_S05_abd1000 <- ASVs_S05_yes %>% filter(total_abd > 1000)
S05_div <- all_ASVs_ref[ASVs_S05_abd1000$OTU_ID]
names(S05_div) <- str_c("S05_div_", 1:length(S05_div))

ASVs_S05_abd1000 <- ASVs_S05_abd1000 %>% 
  mutate(Seq_ID = str_c("S05_div_", 1:length(S05_div)))

# S07
ASVs_S07_abd1000 <- ASVs_S07_yes %>% filter(total_abd > 1000)
S07_div <- all_ASVs_ref[ASVs_S07_abd1000$OTU_ID]
names(S07_div) <- str_c("S07_div_", 1:length(S07_div))

ASVs_S07_abd1000 <- ASVs_S07_abd1000 %>% 
  mutate(Seq_ID = str_c("S07_div_", 1:length(S07_div)))

# S08
ASVs_S08_abd1000 <- ASVs_S08_all %>% filter(total_abd > 1000)
S08_div <- all_ASVs_ref[ASVs_S08_abd1000$OTU_ID]
names(S08_div) <- str_c("S08_div_", 1:length(S08_div))

ASVs_S08_abd1000 <- ASVs_S08_abd1000 %>% 
  mutate(Seq_ID = str_c("S08_div_", 1:length(S08_div)))

# S10
ASVs_S10_abd1000 <- ASVs_S10_yes %>% filter(total_abd > 1000)
S10_div <- all_ASVs_ref[ASVs_S10_abd1000$OTU_ID]
names(S10_div) <- str_c("S10_div_", 1:length(S10_div))

ASVs_S10_abd1000 <- ASVs_S10_abd1000 %>% 
  mutate(Seq_ID = str_c("S10_div_", 1:length(S10_div)))


# writeXStringSet(S01_div, "./1.data/S01_div.fasta")
# writeXStringSet(S02_div, "./1.data/S02_div.fasta")
# writeXStringSet(S03_div, "./1.data/S03_div.fasta")
# writeXStringSet(S04_div, "./1.data/S04_div.fasta")
# writeXStringSet(S05_div, "./1.data/S05_div.fasta")
# writeXStringSet(S07_div, "./1.data/S07_div.fasta")
# writeXStringSet(S08_div, "./1.data/S08_div.fasta")
# writeXStringSet(S10_div, "./1.data/S10_div.fasta")

# done ...


####### part2 mock commu metadata #######
mock_commu_metadata <- read_excel("./1.data/fungal_mock_commu/mock_all_metadata.xlsx", sheet = 1)
# mock_commu_metadata %>% view()

process_vsearch_rlt <- function(SAMPLE_ID) {
  
  all_hit_file <- list.files("./1.data/fungal_mock_commu/1.vsearch_rlts_new_1/vsearch_rlts_new_1/")
  hit_file_name <- all_hit_file[str_detect(all_hit_file, pattern = SAMPLE_ID)][1]
  
  
  # read the hit file
  vsearch_hit_filepath <- str_glue("./1.data/fungal_mock_commu/1.vsearch_rlts_new_1/vsearch_rlts_new_1/{hit_file_name}")
  vsearch_hit <- read_table(vsearch_hit_filepath, col_names = F)
  
  
  all_log_file <- list.files("./1.data/fungal_mock_commu/1.vsearch_rlts_new_1/vsearch_log_new_1/")
  log_file_name <- all_log_file[str_detect(all_log_file, pattern = SAMPLE_ID)]
  
  # read the log file
  vsearch_log_filepath <- str_glue("./1.data/fungal_mock_commu/1.vsearch_rlts_new_1/vsearch_log_new_1/{log_file_name}")
  vsearch_log <- read_log(vsearch_log_filepath, skip = 5)
  
  total_reads <- vsearch_log[1, 7]
  
  
  # summary
  vsearch_hit1 <- vsearch_hit %>% 
    mutate(
      X2_1 = 
        if_else(
          str_count(X2, pattern = "_") == 1,
          str_split_i(X2, pattern = "_", 1),
          str_c(str_split_i(X2, pattern = "_", 1), str_split_i(X2, pattern = "_", 2), sep = "_")
        )
      ) %>% 
    count(X2_1) %>% 
    mutate(sample_reads = total_reads %>% as.numeric()) %>%
    mutate(n_per_0 = n / sample_reads) %>% dplyr::select(-sample_reads) %>% 
    mutate(n_per_1 = n / sum(n))
  
  # add exp n
  strain_list <- mock_commu_metadata %>% 
    filter(sample_id == SAMPLE_ID) %>%
    dplyr::select(strain_list) %>% pull() %>% 
    str_split_1(pattern = "_")
  
  strain_list_num <- length(strain_list)
  
  colnames(vsearch_hit1)[1] <- "strain_id"
  
  vsearch_hit2 <- vsearch_hit1 %>% 
    mutate(n_1 = if_else(strain_id %in% strain_list, n, 0), .after = n) %>% 
    mutate(n_per_2 = n_1 / sum(n_1)) %>% 
    mutate(n_per_exp = if_else(n_1 != 0, 1/strain_list_num, 0)) %>% 
    mutate(n_exp = sum(n_1)*n_per_exp)
  
  # summary
  vsearch_rlts <- vsearch_hit2 %>% mutate(
    sample_id = SAMPLE_ID, .before = strain_id
  )
  
  
  return(vsearch_rlts)
  
  
}


sample_list <- mock_commu_metadata$sample_id %>% as.character()
# sample_list


# map_dfr(sample_list, process_vsearch_rlt)


member_percentage_vsearch <- future_map_dfr(sample_list, process_vsearch_rlt, .progress = T)
# write_xlsx(member_percentage_vsearch, "vsearch_rlts_new_20260802.xlsx")

# member_percentage_vsearch$strain_id %>% table()


# member_percentage_vsearch %>% filter(sample_id == "D3M9a")
# member_percentage_vsearch %>% filter(sample_id == "D3M9b")
# member_percentage_vsearch %>% filter(sample_id == "D3M9c")
#  
#  
# member_percentage_vsearch %>% filter(sample_id == "S3M9a")
# member_percentage_vsearch %>% filter(sample_id == "S3M9b")
# member_percentage_vsearch %>% filter(sample_id == "S3M9c")
#  
#  
# member_percentage_vsearch %>% filter(sample_id == "M3M9a")
# member_percentage_vsearch %>% filter(sample_id == "M3M9b")
# member_percentage_vsearch %>% filter(sample_id == "M3M9c")



###### calcu bias #######
member_bias <- 
  member_percentage_vsearch %>% filter(n_per_2 != 0) %>% 
  group_by(sample_id) %>% 
  mutate(
    D_m = mean(log(n_per_2 / n_per_exp))
  ) %>% 
  mutate(
    bias = log(n_per_2 / n_per_exp),
    bias_clr = bias - D_m
  ) %>% 
  mutate(
    group = str_sub(sample_id, 1, 1)
  )



radar_plot_data_process <- function(DF) {
  
  # levels
  DF$group <- factor(DF$group, levels = c("D", "S", "M"))
  
  DF_sum <- DF %>% group_by(strain_id, group) %>% 
    summarise(bias_clr_m = mean(bias_clr))
  
  DF_sum_rad <- 
    DF_sum %>% pivot_wider(
      names_from = strain_id, values_from = bias_clr_m
    ) %>% mutate(
      group = case_when(
        group == "D" ~ "Equal DNA",
        group == "S" ~ "Equal single copy gene",
        group == "M" ~ "Equal rDNA",
        .default = group
      )
    )
  
  DF_sum_rad$group <- 
    factor(DF_sum_rad$group, levels = c("Equal DNA", "Equal single copy gene", "Equal rDNA"))
    
  return(DF_sum_rad)
  
}


all_commu_bias_radar <- radar_plot_data_process(member_bias)

all_commu_bias_radar %>% dplyr::select(-group) %>% range()
# -1.6550352  0.7347183

##### all commu bias #####
# fig4_all_commu_bias_radar <- 
#   ggradar(all_commu_bias_radar,
#           values.radar = c("-1.7", "0", "1.2"),
#           grid.label.size = 10,
#           #plot.extent.x.sf = 2,
#           #plot.extent.y.sf = 2,
#           grid.min = -1.7,
#           grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
#           grid.max = 1.2,
#           group.line.width = 1, 
#           group.point.size = 5,
#           group.colours = c("#FF007F", "#4455CC", "#008A1A"),
#           axis.label.size = 6,
#           # axis.labels = dPCR_rlts1$Taxa,
#           axis.label.offset = 1.1,
#           legend.title = "Treatment",
#           plot.title = "all commu bias",
#           background.circle.colour = "grey95") +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5, vjust = 0.5)
#   )
# fig4_all_commu_bias_radar



###### AMF all ASV in one ######
member_percentage_vsearch_combAMF <- 
  member_percentage_vsearch %>% 
  dplyr::select(sample_id, strain_id, n, n_1, n_per_exp) %>% 
  mutate(
    strain_id = str_split_i(strain_id, pattern = "_", 1)
  ) %>% group_by(sample_id, strain_id) %>% 
  summarise(
    n = sum(n),
    n_1 = sum(n_1)
  ) %>% 
  mutate(
    n_1 = if_else(n_1 == 0, 0, n)
  ) %>% 
  mutate(
    n_per_2 = n_1 / sum(n_1)
  ) %>% ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(
    n_per_exp = if_else(n_per_2 != 0, 1/sum(n_per_2 != 0), 0)
  ) %>% ungroup()

# member_percentage_vsearch_combAMF %>% filter(str_detect(strain_id, "S08"))

member_bias_combAMF <- 
  member_percentage_vsearch_combAMF %>% filter(n_per_2 != 0) %>% 
  group_by(sample_id) %>% 
  mutate(
    D_m = mean(log(n_per_2 / n_per_exp))
  ) %>% 
  mutate(
    bias = log(n_per_2 / n_per_exp),
    bias_clr = bias - D_m
  ) %>% 
  mutate(
    group = str_sub(sample_id, 1, 1)
  )

all_commu_bias_combAMF_radar <- radar_plot_data_process(member_bias_combAMF)
all_commu_bias_combAMF_radar %>% dplyr::select(-group) %>% range()


sort_strain <- all_commu_bias_combAMF_radar %>% filter(group == "Equal rDNA") %>% pivot_longer(
  names_to = "strain_id", values_to = "bias", -group
) %>% arrange(desc(bias)) %>% select(strain_id) %>% pull()


all_commu_bias_combAMF_radar1 <- all_commu_bias_combAMF_radar %>% 
  select(group, all_of(sort_strain))

fig5d_axis_label <- data.frame(
  strain_id = colnames(all_commu_bias_combAMF_radar1)[-1]
) %>% left_join(dPCR_rlts %>% select(strain_id, Taxa), by = "strain_id") %>% 
  mutate(
    Taxa = str_split_i(Taxa, pattern = " ", 1)
  )

###### 
fig4d_all_commu_bias_combAMF_radar <- 
  ggradar(all_commu_bias_combAMF_radar1,
          values.radar = c("-1.7", "0", "1.2"),
          grid.label.size = 10,
          #plot.extent.x.sf = 2,
          #plot.extent.y.sf = 2,
          grid.min = -1.7,
          grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
          grid.max = 1.2,
          group.line.width = 1, 
          group.point.size = 5,
          # group.colours = c("red", "blue", "darkgreen"),
          group.colours = c("#FF007F", "#4455CC", "#F5B700"),
          axis.label.size = 6,
          axis.labels = fig5d_axis_label$Taxa,
          axis.label.offset = 1.1,
          legend.title = "Method",
          # plot.title = "Mock communities including S08",
          background.circle.colour = "grey95") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5)
  )
fig4d_all_commu_bias_combAMF_radar


tm <- now() %>% str_split_i(pattern = " ", 1)
fig4d_pdf <- str_c("fig4d_", "all_bias_", tm, ".pdf", sep = "")
fig4d_jpg <- str_c("fig4d_", "all_bias_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4d_pdf <- str_c(fig_path, fig4d_pdf)
fig_fullpath_4d_jpg <- str_c(fig_path, fig4d_jpg)

ggsave(fig_fullpath_4d_pdf, fig4d_all_commu_bias_combAMF_radar, width = 7.96, height = 7.21)
ggsave(fig_fullpath_4d_jpg, fig4d_all_commu_bias_combAMF_radar, width = 7.96, height = 7.21)


###### use combAMF ######
fig4_all_commu_bias_combAMF_radar


###### AMF mock vs no AMF mock ######
# mock_commu_metadata

sample_AMF <- mock_commu_metadata %>% filter(
  str_detect(strain_list, pattern = "S08")
) %>% dplyr::select(sample_id) %>% pull()

sample_others <- mock_commu_metadata %>% filter(
  !str_detect(strain_list, pattern = "S08")
) %>% dplyr::select(sample_id) %>% pull()


######## AMF bias
# member_bias_combAMF

member_bias_combAMF_AMF <- member_bias_combAMF %>% 
  filter(sample_id %in% sample_AMF)

AMF_commu_radar <- radar_plot_data_process(member_bias_combAMF_AMF)


# fig4_AMF_commu_radar <- 
#   ggradar(AMF_commu_radar,
#           values.radar = c("-1.7", "0", "1.2"),
#           grid.label.size = 10,
#           #plot.extent.x.sf = 2,
#           #plot.extent.y.sf = 2,
#           grid.min = -1.7,
#           grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
#           grid.max = 1.2,
#           group.line.width = 1, 
#           group.point.size = 5,
#           # group.colours = c("red", "blue", "darkgreen"),
#           group.colours = c("#FF007F", "#4455CC", "#008A1A"),
#           axis.label.size = 6,
#           # axis.labels = dPCR_rlts1$Taxa,
#           axis.label.offset = 1.1,
#           legend.title = "Treatment",
#           plot.title = "Bias in mock communities with S08",
#           background.circle.colour = "grey95") +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5, vjust = 0.5)
#   )
# fig4_AMF_commu_radar


##### other bias
member_bias_combAMF_others <- member_bias_combAMF %>% 
  filter(sample_id %in% sample_others)

other_commu_radar <- radar_plot_data_process(member_bias_combAMF_others)

other_commu_radar1 <- other_commu_radar %>% 
  select(group, S07, S03, S05, S02, S01, S04)

# dPCR_rlts

fig5c_axis_label <- data.frame(
  strain_id = colnames(other_commu_radar1)[-1]
) %>% left_join(dPCR_rlts %>% select(strain_id, Taxa), by = "strain_id") %>% 
  mutate(Taxa = str_split_i(Taxa, pattern = " ", 1))

fig4c_other_commu_radar <- 
  ggradar(other_commu_radar1,
          values.radar = c("-1.7", "0", "1.2"),
          grid.label.size = 10,
          #plot.extent.x.sf = 2,
          #plot.extent.y.sf = 2,
          grid.min = -1.7,
          grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
          grid.max = 1.2,
          group.line.width = 1, 
          group.point.size = 5,
          # group.colours = c("red", "blue", "darkgreen"),
          group.colours = c("#FF007F", "#4455CC", "#F5B700"),
          axis.label.size = 6,
          axis.labels = fig5c_axis_label$Taxa,
          axis.label.offset = 1.1,
          legend.title = "Method",
          # plot.title = "Mock communities excluding S08",
          background.circle.colour = "grey95") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5)
  )
fig4c_other_commu_radar
# fig4_radar_noAMF_AMF <- (fig4_other_commu_radar + fig4_AMF_commu_radar) / guide_area() + plot_layout(guides = "collect")

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4c_pdf <- str_c("fig4c_", "bias_without_S08_", tm, ".pdf", sep = "")
fig4c_jpg <- str_c("fig4c_", "bias_without_S08_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4c_pdf <- str_c(fig_path, fig4c_pdf)
fig_fullpath_4c_jpg <- str_c(fig_path, fig4c_jpg)

ggsave(fig_fullpath_4c_pdf, fig4c_other_commu_radar, width = 7.96, height = 7.21)
ggsave(fig_fullpath_4c_jpg, fig4c_other_commu_radar, width = 7.96, height = 7.21)



###### which one is robust? ######
strain_pool_check <- c("S01", "S02", "S03", "S04", "S05", "S07")

check_g2 <- data.frame(
  strain_number = 2,
  strain_list = combn(strain_pool_check, 2) %>% as.data.frame() %>% map_chr(~ str_c(., collapse = "_"))
)
# check_g2


check_bias_calcu <- function(STRAIN_LIST){
  
  strain_1 <- str_split_i(STRAIN_LIST, pattern = "_", 1)
  strain_2 <- str_split_i(STRAIN_LIST, pattern = "_", 2)
  
  # strain_1 and strain_2
  sample_list_tmp <- 
    mock_commu_metadata %>% 
    filter(str_detect(strain_list, strain_1)) %>% 
    filter(str_detect(strain_list, strain_2)) %>% 
    dplyr::select(sample_id) %>% pull()
  # sample_list_tmp
  
  member_per_tmp <- 
    member_percentage_vsearch_combAMF %>% filter(sample_id %in% sample_list_tmp) %>% 
    filter(strain_id == strain_1 | strain_id == strain_2) %>% 
    dplyr::select(sample_id, strain_id, n) %>% 
    group_by(sample_id) %>% 
    mutate(
      n_per = n/sum(n)
    ) %>% 
    mutate(
      n_per_exp = 0.5
    ) %>% 
    mutate(
      bias_clr = log(n_per / n_per_exp) - mean(log(n_per / n_per_exp))
    )
  
  return(member_per_tmp)
  
}


check_bias <- future_map_dfr(check_g2$strain_list, check_bias_calcu, .progress = T)


check_bias1 <- check_bias %>% mutate(
  group = str_sub(sample_id, 1, 1)
)

# check_bias1

check_bias1_radar <- radar_plot_data_process(check_bias1)

# fig4_check_radar <- 
#   ggradar(check_bias1_radar,
#           values.radar = c("-1.7", "0", "1.2"),
#           grid.label.size = 10,
#           #plot.extent.x.sf = 2,
#           #plot.extent.y.sf = 2,
#           grid.min = -1.7,
#           grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
#           grid.max = 1.2,
#           group.line.width = 1, 
#           group.point.size = 5,
#           group.colours = c("red", "blue", "darkgreen"),
#           axis.label.size = 6,
#           # axis.labels = dPCR_rlts1$Taxa,
#           axis.label.offset = 1.1,
#           legend.title = "Treatment",
#           plot.title = "Check bias in all commu",
#           background.circle.colour = "grey95") +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold")
#   )
# fig4_check_radar


# no S10 and S08
check_bias_calcu_without <- function(STRAIN_LIST){
  
  strain_1 <- str_split_i(STRAIN_LIST, pattern = "_", 1)
  strain_2 <- str_split_i(STRAIN_LIST, pattern = "_", 2)
  
  # strain_1 and strain_2
  sample_list_tmp <- 
    mock_commu_metadata %>% 
    filter(str_detect(strain_list, strain_1)) %>% 
    filter(str_detect(strain_list, strain_2)) %>% 
    filter(!str_detect(strain_list, "S10")) %>%
    filter(!str_detect(strain_list, "S08")) %>% 
    dplyr::select(sample_id) %>% pull()
  # sample_list_tmp
  
  member_per_tmp <- 
    member_percentage_vsearch_combAMF %>% filter(sample_id %in% sample_list_tmp) %>% 
    filter(strain_id == strain_1 | strain_id == strain_2) %>% 
    dplyr::select(sample_id, strain_id, n) %>% 
    group_by(sample_id) %>% 
    mutate(
      n_per = n/sum(n)
    ) %>% 
    mutate(
      n_per_exp = 0.5
    ) %>% 
    mutate(
      bias_clr = log(n_per / n_per_exp) - mean(log(n_per / n_per_exp))
    )
  
  return(member_per_tmp)
  
}

# check_bias_calcu


check_bias_without <- future_map_dfr(check_g2$strain_list, check_bias_calcu_without, .progress = T)
# check_bias_without

check_bias_without1 <- check_bias_without %>% mutate(
  group = str_sub(sample_id, 1, 1)
)

# radar
check_without1_radar <- radar_plot_data_process(check_bias_without1)


# fig4_check_no_S01_S08_radar <- 
#   ggradar(check_without1_radar,
#           values.radar = c("-1.7", "0", "1.2"),
#           grid.label.size = 10,
#           #plot.extent.x.sf = 2,
#           #plot.extent.y.sf = 2,
#           grid.min = -1.7,
#           grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
#           grid.max = 1.2,
#           group.line.width = 1, 
#           group.point.size = 5,
#           group.colours = c("red", "blue", "darkgreen"),
#           axis.label.size = 6,
#           # axis.labels = dPCR_rlts1$Taxa,
#           axis.label.offset = 1.1,
#           legend.title = "Treatment",
#           plot.title = "Check bias in mock commu without S01 and S08",
#           background.circle.colour = "grey95") +
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold")
#   )
# fig4_check_no_S01_S08_radar
# 
# fig4_check_radar + fig4_check_no_S01_S08_radar


###### bias among S02, S05, S08, S10 ######
strain_pool_adj <- c("S02", "S05", "S08", "S10")

adj_g2 <- data.frame(
  strain_number = 2,
  strain_list = combn(strain_pool_adj, 2) %>% as.data.frame() %>% map_chr(~ str_c(., collapse = "_"))
)

adj_g2a <- 
  adj_g2 %>% filter(
    !(str_detect(strain_list, pattern = "S02")  &  str_detect(strain_list, pattern = "S05"))
  ) %>% filter(
    !(str_detect(strain_list, pattern = "S08")  &  str_detect(strain_list, pattern = "S10"))
  )


adj_bias_calcu <- function(STRAIN_LIST){
  
  strain_1 <- str_split_i(STRAIN_LIST, pattern = "_", 1)
  strain_2 <- str_split_i(STRAIN_LIST, pattern = "_", 2)
  
  # strain_1 and strain_2
  sample_list_tmp <- 
    mock_commu_metadata %>% 
    filter(str_detect(strain_list, strain_1)) %>% 
    filter(str_detect(strain_list, strain_2)) %>% 
    dplyr::select(sample_id) %>% pull()
  # sample_list_tmp
  
  member_per_tmp <- 
    member_percentage_vsearch_combAMF %>% filter(sample_id %in% sample_list_tmp) %>% 
    filter(strain_id == strain_1 | strain_id == strain_2) %>% 
    dplyr::select(sample_id, strain_id, n) %>% 
    group_by(sample_id) %>% 
    mutate(
      n_per = n/sum(n)
    ) %>% 
    mutate(
      n_per_exp = 0.5
    ) %>% 
    mutate(
      bias_clr = log(n_per / n_per_exp) - mean(log(n_per / n_per_exp))
    )
  
  return(member_per_tmp)
  
}


adj_bias <- future_map_dfr(adj_g2a$strain_list, adj_bias_calcu, .progress = T)
# adj_bias

adj_bias1 <- adj_bias %>% mutate(
  group = str_sub(sample_id, 1, 1)
)


# radar_plot_data_process(adj_bias1)


# adj_bias1 %>% filter(group == "M") %>% 
#   filter(strain_id == "S08" | strain_id == "S10") %>% 
#   ggplot() +
#   geom_boxplot(aes(strain_id, bias_clr)) +
#   geom_point(aes(strain_id, bias_clr))

adj_value_S08_S10 <- 
  adj_bias1 %>% filter(group == "M") %>% 
  filter(strain_id == "S08" | strain_id == "S10") %>% 
  group_by(strain_id) %>% 
  summarise(
    bias_clr_m = mean(bias_clr),
    n_per_m = mean(n_per),
    n_per_m_paired = 1 - n_per_m
  ) %>% 
  mutate(
    adj_value = n_per_m / 0.5
  ) %>% dplyr::select(strain_id, adj_value)



member_percentage_vsearch_adj_exp <- 
  member_percentage_vsearch_combAMF %>% 
  group_by(sample_id) %>% 
  mutate(group = str_sub(sample_id, 1, 1)) %>% 
  mutate(
    n_per_exp = if_else(
      (group == "M" & strain_id == "S08"),
      n_per_exp*adj_value_S08_S10$adj_value[1],
      n_per_exp
    )
  ) %>% 
  mutate(
    n_per_exp = if_else(
      (group == "M" & strain_id == "S10"),
      n_per_exp*adj_value_S08_S10$adj_value[2],
      n_per_exp
    )
  ) %>% left_join(
    mock_commu_metadata, by = "sample_id"
  ) %>% 
  mutate(
    n_per_exp = if_else(n_per_exp != 0 & strain_id != "S08" & str_detect(strain_list, "S08") & !str_detect(strain_list, "S10") & group == "M",
                        (1 - (1/strain_num*adj_value_S08_S10$adj_value[1])) / (strain_num - 1),
                        n_per_exp)
  ) %>% 
  mutate(
    n_per_exp = if_else(n_per_exp != 0 & strain_id != "S10" & !str_detect(strain_list, "S08") & str_detect(strain_list, "S10") & group == "M",
                        (1 - (1/strain_num*adj_value_S08_S10$adj_value[2])) / (strain_num - 1),
                        n_per_exp)
  ) %>% 
  mutate(
    n_per_exp = if_else(n_per_exp != 0 & strain_id != "S10" & strain_id != "S08" & str_detect(strain_list, "S08") & str_detect(strain_list, "S10")  & group == "M",
                        (1 - (1/strain_num*adj_value_S08_S10$adj_value[1] +
                                1/strain_num*adj_value_S08_S10$adj_value[2])) / (strain_num - 2),
                        n_per_exp)
  ) %>% filter(n_per_2 != 0) %>% 
  mutate(
    bias_clr = log(n_per_2 / n_per_exp) - mean(log(n_per_2 / n_per_exp))
  )

member_percentage_vsearch_adj_exp


# ggplot(member_percentage_vsearch_adj_exp %>% 
#          filter(group == "S" & strain_id == "S08") %>% 
#          mutate(strain_list = str_remove(strain_list, pattern = "S08")), aes(strain_list, bias_clr)) +
#   geom_point()



# radar 
bias_adj_radar <- radar_plot_data_process(member_percentage_vsearch_adj_exp)

bias_adj_radar1 <- bias_adj_radar %>% 
  select(group, all_of(sort_strain))

fig5e_axis_label <- data.frame(
  strain_id = colnames(bias_adj_radar1)[-1]
) %>% left_join(dPCR_rlts %>% select(strain_id, Taxa), by = "strain_id") %>% 
  mutate(
    Taxa = str_split_i(Taxa, pattern = " ", 1)
  )

fig4e_bias_adj_radar <- 
  ggradar(bias_adj_radar1,
          values.radar = c("-1.7", "0", "1.2"),
          grid.label.size = 10,
          #plot.extent.x.sf = 2,
          #plot.extent.y.sf = 2,
          grid.min = -1.7,
          grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
          grid.max = 1.2,
          group.line.width = 1, 
          group.point.size = 5,
          # group.colours = c("red", "blue", "darkgreen"),
          group.colours = c("#FF007F", "#4455CC", "#F5B700"),
          axis.label.size = 6,
          axis.labels = fig5e_axis_label$Taxa,
          axis.label.offset = 1.1,
          # plot.title = "Mock communities with S08 correction",
          legend.title = "Method",
          background.circle.colour = "grey95") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5)
  )
fig4e_bias_adj_radar

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4e_pdf <- str_c("fig4e_", "bias_corrected_", tm, ".pdf", sep = "")
fig4e_jpg <- str_c("fig4e_", "bias_corrected_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4e_pdf <- str_c(fig_path, fig4e_pdf)
fig_fullpath_4e_jpg <- str_c(fig_path, fig4e_jpg)

ggsave(fig_fullpath_4e_pdf, fig4e_bias_adj_radar, width = 7.96, height = 7.21)
ggsave(fig_fullpath_4e_jpg, fig4e_bias_adj_radar, width = 7.96, height = 7.21)


# fig4_all_commu_bias_combAMF_radar + fig4_bias_adj_radar

# dPCR_rlts


####### ASVs div ########


# UC
S01_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S01_div.uc", delim = "\t", col_names = F)
S02_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S02_div.uc", delim = "\t", col_names = F)
S03_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S03_div.uc", delim = "\t", col_names = F)
S04_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S04_div.uc", delim = "\t", col_names = F)
S05_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S05_div.uc", delim = "\t", col_names = F)
S07_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S07_div.uc", delim = "\t", col_names = F)
S08_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S08_div.uc", delim = "\t", col_names = F)
S10_div_UC <- read_delim("./1.data/fungal_mock_commu/2.emboss_aln/S10_div.uc", delim = "\t", col_names = F)


# parse ur data
parse_ur_data <- function(UC, ABD){
  
  # data
  uc_tmp <- UC
  abd_tmp <- ABD %>% dplyr::select(Seq_ID, total_abd)
  
  # cluster
  C_list <- uc_tmp %>% filter(X1 == "C") %>% dplyr::select(X2, X9)
  H_list <- uc_tmp %>% filter(X1 == "H") %>% dplyr::select(X2, X9)
  
  C_list_H <- C_list %>% filter(X2 %in% H_list$X2)
  C_list_noH <- C_list %>% filter(!X2 %in% H_list$X2)
  
  # abd summary H
  H_abd <- H_list %>% 
    left_join(abd_tmp, by = c("X9" = "Seq_ID")) %>% 
    group_by(X2) %>% summarise(total_abd = sum(total_abd))
  
  # abd comb
  C_list_H_abd <- C_list_H %>% left_join(H_abd, by = "X2")
  C_list_noH_abd <- C_list_noH %>% left_join(abd_tmp, by = c("X9" = "Seq_ID"))
  
  # comb rlts
  abd_comb_tmp <- rbind(C_list_H_abd, C_list_noH_abd)
  abd_rlts <- abd_comb_tmp %>% mutate(
    X2 = str_c("clr_", X2)
  )
  colnames(abd_rlts) <- c("clr", "Seq_ID", "reads_num")
  
  # return
  return(abd_rlts)
  
}


#### net plot
# node data
S01_node <- parse_ur_data(UC = S01_div_UC, ABD = ASVs_S01_abd1000)
S02_node <- parse_ur_data(UC = S02_div_UC, ABD = ASVs_S02_abd1000)
S03_node <- parse_ur_data(UC = S03_div_UC, ABD = ASVs_S03_abd1000)
S04_node <- parse_ur_data(UC = S04_div_UC, ABD = ASVs_S04_abd1000)
S05_node <- parse_ur_data(UC = S05_div_UC, ABD = ASVs_S05_abd1000)
S07_node <- parse_ur_data(UC = S07_div_UC, ABD = ASVs_S07_abd1000)
S08_node <- parse_ur_data(UC = S08_div_UC, ABD = ASVs_S08_abd1000)


S08_clr <- S08_div[S08_node$Seq_ID]
names(S08_clr) <- str_c("S08_", S08_node$clr)
# S08_clr

S08_node_taxa <- S08_node %>%
  left_join(ASVs_S08_abd1000 %>% filter(Seq_ID %in% S08_node$Seq_ID), by = "Seq_ID")


# writeXStringSet(S08_clr, "S08_clr.fasta")
# write_xlsx(S08_node_taxa, "S08_clr_taxa.xlsx")

S10_node <- parse_ur_data(UC = S10_div_UC, ABD = ASVs_S10_abd1000)


all_node_data <- rbind(S01_node, S02_node, S03_node, S04_node,
                       S05_node, S07_node, S08_node, S10_node)

node_range <- range(all_node_data$reads_num)



# vsearch_rlts
S01_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S01_v_aln.txt", delim = "\t", col_names = F)
S02_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S02_v_aln.txt", delim = "\t", col_names = F)
S03_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S03_v_aln.txt", delim = "\t", col_names = F)
S04_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S04_v_aln.txt", delim = "\t", col_names = F)
S05_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S05_v_aln.txt", delim = "\t", col_names = F)
S07_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S07_v_aln.txt", delim = "\t", col_names = F)
S08_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S08_v_aln.txt", delim = "\t", col_names = F)
S10_v_sime <- read_delim("./1.data/fungal_mock_commu/2.vsearch_aln/S10_v_aln.txt", delim = "\t", col_names = F)

all_v_data <- rbind(S01_v_sime, S02_v_sime, S03_v_sime, S04_v_sime,
                    S05_v_sime, S07_v_sime, S08_v_sime, S10_v_sime) %>% 
  mutate(X3 = X3 / 100)

v_range <- range(all_v_data$X3)


parse_vsearch_paired_aln_rlts <- function(ALN, NODE) {
  
  # ALN data
  ALN_tmp <- ALN %>% dplyr::select(1:3)
  
  colnames(ALN_tmp) <- c("seq1", "seq2", "pident")
  
  # matrix
  MAT_tmp <- 
    ALN_tmp %>%
    bind_rows(ALN_tmp %>% dplyr::rename(seq1 = seq2, seq2 = seq1)) %>%
    complete(seq1, seq2) %>%
    pivot_wider(names_from = seq2, values_from = pident) %>%
    column_to_rownames("seq1") %>%
    as.matrix()
  
  diag(MAT_tmp) <- 100
  ALN_tmp_dis <- 100 - MAT_tmp
  
  # PCoA analysis
  TMP_pcoa_result <- cmdscale(ALN_tmp_dis, k = 2, eig = TRUE)
  
  # tmp
  # return(TMP_pcoa_result)
  
  TMP_points <- as.data.frame(TMP_pcoa_result$points)
  
  colnames(TMP_points) <- c("PCoA1", "PCoA2")
  
  TMP_points1 <- TMP_points %>% mutate(
    Seq_ID = rownames(.)
  ) %>% left_join(
    NODE, by = "Seq_ID")
  
  # eig
  all_eig <- TMP_pcoa_result$eig
  pos_eig <- all_eig[all_eig > 0]
  
  TMP_var_explained <- round(pos_eig / sum(pos_eig) * 100, 2)[1:2]
  
  TMP_seg_comb <- combn(1:nrow(TMP_points), 2)
  
  # segments data
  TMP_segments <- data.frame(
    from = rownames(TMP_points)[TMP_seg_comb[1, ]],
    to = rownames(TMP_points)[TMP_seg_comb[2, ]],
    x = TMP_points$PCoA1[TMP_seg_comb[1, ]],
    y = TMP_points$PCoA2[TMP_seg_comb[1, ]],
    xend = TMP_points$PCoA1[TMP_seg_comb[2, ]],
    yend = TMP_points$PCoA2[TMP_seg_comb[2, ]]
  )
  
  # with mid
  TMP_segments_withmid <- TMP_segments %>% 
    mutate(
      mid_x = (x + xend) / 2,
      mid_y = (y + yend) / 2
    ) %>% 
    mutate(
      FT_key = str_c(
        pmin(from, to),
        pmax(from, to), sep = "|")
    ) %>% 
    left_join(
      ALN_tmp %>% mutate(
        FT_key = str_c(
          pmin(seq1, seq2),
          pmax(seq1, seq2), sep = "|")
      ) %>% dplyr::select(FT_key, pident), by = "FT_key"
    ) %>% 
    mutate(
      pident = pident / 100
    )
  
  PCoA1_range <- range(TMP_points1$PCoA1)
  PCoA2_range <- range(TMP_points1$PCoA2)
  
  spc_lab_tmp <- c((PCoA1_range[2] - PCoA1_range[1]) / 2 + PCoA1_range[1],
                   (PCoA2_range[2] - PCoA2_range[1]) / 2 + PCoA2_range[1])
  
  # rlts
  pcoa_identity_rlts <- list(
    main_df = TMP_points1, 
    sub_seg = TMP_segments_withmid,
    axis_eig = TMP_var_explained,
    spc_lab = spc_lab_tmp
  )
  
  return(pcoa_identity_rlts)
  
}


# S01
S01_pcoa_data <- parse_vsearch_paired_aln_rlts(ALN = S01_v_sime, NODE = S01_node)

S01_pcoa_data$main_df$ASVs_ID <- str_replace(S01_pcoa_data$main_df$Seq_ID, "S01_div", "ASV")

S01_pcoa_plot <- 
  ggplot() +
  geom_segment(
    data = S01_pcoa_data$sub_seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pident), 
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = S01_pcoa_data$sub_seg,
    aes(x = mid_x, y = mid_y, label = pident),
    size = 3,
    max.overlaps = 5
  ) +
  geom_point(
    data = S01_pcoa_data$main_df, 
    aes(PCoA1, PCoA2, size = reads_num),
    colour = "lightblue", alpha = 0.8) +
  geom_text_repel(
    data = S01_pcoa_data$main_df,
    aes(PCoA1, PCoA2, label = ASVs_ID),
    colour = "black", size = 5, fontface = "bold"
  ) +
  # geom_text(
  #   aes(x = 0, y = -7, label = "Num. of cluster: 6")
  # ) +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_colour_gradient2(limits = v_range, , name = "Identity",
                         low = "#FF007F", high = "grey90", mid = "grey99",
                         midpoint = (v_range[1] + v_range[2]) / 2) +
  labs(x = str_c("PCoA1 (", S01_pcoa_data$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", S01_pcoa_data$axis_eig[2], "%)"),
       title = "Pleurotus ostreatus") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1)
S01_pcoa_plot


# S02
S02_pcoa_data <- parse_vsearch_paired_aln_rlts(ALN = S02_v_sime, NODE = S02_node)

S02_pcoa_data$main_df$ASVs_ID <- str_replace(S02_pcoa_data$main_df$Seq_ID, "S02_div", "ASV")

S02_pcoa_plot <- 
  ggplot() +
  geom_segment(
    data = S02_pcoa_data$sub_seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pident), 
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = S02_pcoa_data$sub_seg,
    aes(x = mid_x, y = mid_y, label = pident),
    size = 3,
    max.overlaps = 5
  ) +
  geom_point(
    data = S02_pcoa_data$main_df, 
    aes(PCoA1, PCoA2, size = reads_num),
    colour = "lightblue", alpha = 0.8) +
  geom_text_repel(
    data = S02_pcoa_data$main_df,
    aes(PCoA1, PCoA2, label = ASVs_ID),
    colour = "black", size = 5, fontface = "bold"
  ) +
  # geom_text(
  #   aes(x = 0, y = -7, label = "Num. of cluster: 6")
  # ) +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_colour_gradient2(limits = v_range, , name = "Identity",
                         low = "#FF007F", high = "grey90", mid = "grey99",
                         midpoint = (v_range[1] + v_range[2]) / 2) +
  labs(x = str_c("PCoA1 (", S02_pcoa_data$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", S02_pcoa_data$axis_eig[2], "%)"),
       title = "Saccharomyces cerevisiae") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1)
S02_pcoa_plot


# S03
S03_pcoa_data <- parse_vsearch_paired_aln_rlts(ALN = S03_v_sime, NODE = S03_node)

S03_pcoa_data$main_df$ASVs_ID <- str_replace(S03_pcoa_data$main_df$Seq_ID, "S03_div", "ASV")

S03_pcoa_plot <- 
  ggplot() +
  geom_segment(
    data = S03_pcoa_data$sub_seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pident), 
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = S03_pcoa_data$sub_seg,
    aes(x = mid_x, y = mid_y, label = pident),
    size = 3,
    max.overlaps = 5
  ) +
  geom_point(
    data = S03_pcoa_data$main_df, 
    aes(PCoA1, PCoA2, size = reads_num),
    colour = "lightblue", alpha = 0.8) +
  geom_text_repel(
    data = S03_pcoa_data$main_df,
    aes(PCoA1, PCoA2, label = ASVs_ID),
    colour = "black", size = 5, fontface = "bold"
  ) +
  # geom_text(
  #   aes(x = 0, y = -7, label = "Num. of cluster: 6")
  # ) +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_colour_gradient2(limits = v_range, , name = "Identity",
                         low = "#FF007F", high = "grey90", mid = "grey99",
                         midpoint = (v_range[1] + v_range[2]) / 2) +
  labs(x = str_c("PCoA1 (", S03_pcoa_data$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", S03_pcoa_data$axis_eig[2], "%)"),
       title = "Neurospora crassa") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1)
S03_pcoa_plot



# S04
S04_node

S04_plot_data <- S04_node %>% 
  mutate(
    x = 1,
    y = 1
  )


S04_plot_data$ASVs_ID <- str_replace(S04_plot_data$Seq_ID, "S04_div", "ASV")

S04_net_p <- 
  ggplot(S04_plot_data) +
  geom_point(aes(x, y, size = reads_num), colour = "lightblue", alpha = 0.8, show.legend = F) +
  geom_text_repel(aes(x, y, label = ASVs_ID), show.legend = F, colour = "black", size = 5, fontface = "bold") +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL, title = "Podila humilis") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1,
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    axis.ticks = element_blank(),
    axis.text = element_blank())
S04_net_p

# S05
S05_node

S05_plot_data <- S05_node %>% 
  mutate(
    x = 1,
    y = 1
  )

S05_plot_data$ASVs_ID <- str_replace(S05_plot_data$Seq_ID, "S05_div", "ASV")


S05_net_p <- 
  ggplot(S05_plot_data) +
  geom_point(aes(x, y, size = reads_num), colour = "lightblue", alpha = 0.8, show.legend = F) +
  geom_text_repel(aes(x, y, label = ASVs_ID), show.legend = F, colour = "black", size = 5, fontface = "bold") +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL, title = "Fusarium proliferatum") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1,
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    axis.ticks = element_blank(),
    axis.text = element_blank())
S05_net_p


# S07
S07_pcoa_data <- parse_vsearch_paired_aln_rlts(ALN = S07_v_sime, NODE = S07_node)

S07_pcoa_data$main_df$ASVs_ID <- str_replace(S07_pcoa_data$main_df$Seq_ID, "S07_div", "ASV")

S07_pcoa_plot <- 
  ggplot() +
  geom_segment(
    data = S07_pcoa_data$sub_seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pident), 
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = S07_pcoa_data$sub_seg,
    aes(x = mid_x, y = mid_y, label = pident),
    size = 3,
    max.overlaps = 5
  ) +
  geom_point(
    data = S07_pcoa_data$main_df, 
    aes(PCoA1, PCoA2, size = reads_num),
    colour = "lightblue", alpha = 0.8) +
  geom_text_repel(
    data = S07_pcoa_data$main_df,
    aes(PCoA1, PCoA2, label = ASVs_ID),
    colour = "black", size = 5, fontface = "bold"
  ) +
  # geom_text(
  #   aes(x = 0, y = -7, label = "Num. of cluster: 6")
  # ) +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_colour_gradient2(limits = v_range, , name = "Identity",
                         low = "#FF007F", high = "grey90", mid = "grey99",
                         midpoint = (v_range[1] + v_range[2]) / 2) +
  labs(x = str_c("PCoA1 (", S07_pcoa_data$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", S07_pcoa_data$axis_eig[2], "%)"),
       title = "Gongronella butleri") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1)
S07_pcoa_plot


# S08
S08_pcoa_data <- parse_vsearch_paired_aln_rlts(ALN = S08_v_sime, NODE = S08_node)

S08_pcoa_data$main_df$ASVs_ID <- str_replace(S08_pcoa_data$main_df$Seq_ID, "S08_div", "ASV")

S08_pcoa_plot <- 
  ggplot() +
  geom_segment(
    data = S08_pcoa_data$sub_seg,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pident), 
    linewidth = 0.6
  ) +
  geom_text_repel(
    data = S08_pcoa_data$sub_seg,
    aes(x = mid_x, y = mid_y, label = pident),
    size = 3,
    max.overlaps = 5
  ) +
  geom_point(
    data = S08_pcoa_data$main_df, 
    aes(PCoA1, PCoA2, size = reads_num),
    colour = "lightblue", alpha = 0.8) +
  geom_text_repel(
    data = S08_pcoa_data$main_df,
    aes(PCoA1, PCoA2, label = ASVs_ID),
    colour = "black", size = 5, fontface = "bold"
  ) +
  # geom_text(
  #   aes(x = 0, y = -7, label = "Num. of cluster: 6")
  # ) +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_colour_gradient2(limits = v_range, , name = "Identity",
                         low = "#FF007F", high = "grey90", mid = "grey99",
                         midpoint = (v_range[1] + v_range[2]) / 2) +
  labs(x = str_c("PCoA1 (", S08_pcoa_data$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", S08_pcoa_data$axis_eig[2], "%)"),
       title = "Rhizophagus irregularis") +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1)
S08_pcoa_plot


# S10
S10_node

S10_plot_data <- S10_node %>% 
  mutate(
    x = 1,
    y = 1
  )

S10_plot_data$ASVs_ID <- str_replace(S10_plot_data$Seq_ID, "S10_div", "ASV")

S10_net_p <- 
  ggplot(S10_plot_data) +
  geom_point(aes(x, y, size = reads_num), colour = "lightblue", alpha = 0.8, show.legend = F) +
  geom_text_repel(aes(x, y, label = ASVs_ID), show.legend = F, colour = "black", size = 5, fontface = "bold") +
  scale_size(range = c(2, 12), limits = node_range, name = "Reads number",
             breaks = seq(100000, 1800000, 300000)) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  scale_y_continuous(limits = c(0.5, 1.5)) +
  labs(x = NULL, y = NULL, title = "Agaricus bisporus") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 13, colour = "black", face = "bold"),
    aspect.ratio = 1,
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    axis.ticks = element_blank(),
    axis.text = element_blank())
S10_net_p



fig4b_pcoa <- 
  S08_pcoa_plot + S07_pcoa_plot + S03_pcoa_plot + S05_net_p +
  S02_pcoa_plot +  S01_pcoa_plot + S04_net_p + S10_net_p + 
  plot_layout(guides = "collect", nrow = 1)
fig4b_pcoa


tm <- now() %>% str_split_i(pattern = " ", 1)
fig4b_pdf <- str_c("fig4b_", "member_pcoa_", tm, ".pdf", sep = "")
fig4b_jpg <- str_c("fig4b_", "member_pcoa_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4b_pdf <- str_c(fig_path, fig4b_pdf)
fig_fullpath_4b_jpg <- str_c(fig_path, fig4b_jpg)

ggsave(fig_fullpath_4b_pdf, fig4b_pcoa, width = 30, height = 9.18, limitsize = FALSE)
ggsave(fig_fullpath_4b_jpg, fig4b_pcoa, width = 30, height = 9.18, limitsize = FALSE)



###### amplicon vs metatranscriptom ######

# FRRN
FRRN_cla_tab %>% view()
FRRN_phy_tab %>% view()

FRRN_phy_tab %>% filter(
  phy %in% (FRRN_rlt_taxa_spl_grp %>% filter(
    grp == "Main_taxa"
  ) %>% mutate(
    main_phy = str_sub(phy, start = 3)
  ) %>% select(main_phy) %>% pull())
) %>% select(rDNA_m) %>% pull() %>% mean()

FRRN_cla_mean <- mean(FRRN_cla_tab$rDNA_m)


# metagenome and metatrans
EPICON_MG_MT <- read_tsv("./1.data/EPICON/taxa_table_c_abundance.tsv")
colnames(EPICON_MG_MT)

# amplicon
EPICON_fungal_otutab1 %>% dim()

colnames(EPICON_fungal_otutab1)

amplicon_sample_list <- data.frame(
  Sample_ID = colnames(EPICON_fungal_otutab1),
  Methods = "Amplicon"
)



amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11Z21")
)
## TP11Z21a ~ TP11Z21

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R3")
)
amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R03")
)
## TP11R03z ~ TP11R3

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R22")
)
## TP11R22z ~ TP11R22

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R1")
)
amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R01")
)
## TP11R01z ~ TP11R1

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP15R1")
)
## TP15R1t ~ TP15R1

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP02R11")
)
## TP02R11t ~ TP02R11

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP03R2")
)
## TP03R2r ~ TP03R2

amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R2")
)
amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11R02")
)
## TP11R02z ~ TP11R2
amplicon_sample_list %>% filter(
  str_detect(Sample_ID, pattern = "TP11Z21")
)
## TP11Z21a ~ TP11Z21



# metadata
EPICON_metadata <- read_csv("./1.data/EPICON/group.csv")

EPICON_metadata1 <- EPICON_metadata %>% 
  mutate(
    Methods = str_split_i(Group, pattern = "_", 1)
  ) %>% 
  mutate(
    Sample_ID = str_split_i(Group2, pattern = "_", 1)
  ) %>% 
  mutate(
    Sample_ID_1 = case_when(
      Sample_ID == "TP11Z21" ~ "TP11Z21a",
      Sample_ID == "TP11R3" ~ "TP11R03z",
      Sample_ID == "TP11R22" ~ "TP11R22z",
      Sample_ID == "TP11R1" ~ "TP11R01z",
      Sample_ID == "TP15R1" ~ "TP15R1t",
      Sample_ID == "TP02R11" ~ "TP02R11t",
      Sample_ID == "TP03R2" ~ "TP03R2r",
      Sample_ID == "TP11R2" ~ "TP11R02z",
      Sample_ID == "TP11Z21" ~ "TP11Z21a",
      .default = Sample_ID
    )
  ) %>% 
  left_join(
    amplicon_sample_list, by = c("Sample_ID_1" = "Sample_ID") 
  )



######## amp vs metatranscirptom #####
EPICON_metadata_amp_met <- EPICON_metadata1 %>% filter(
  Methods.x == "Metatranscriptome" & Methods.y == "Amplicon"
) %>% distinct(Sample_ID, .keep_all = T)

EPICON_metadata_amp_metRoot <- 
  EPICON_metadata_amp_met %>% filter(
    str_detect(Group, pattern = "Root")
  ) %>% mutate(
    TIME = str_split_i(Group, pattern = "_", 4)
  ) %>% mutate(
    Treatment = str_split_i(Group, patter = "_", 3)
  )



# amp
EPICON_amp1 <- EPICON_fungal_otutab1_1 %>% 
  dplyr::select(OTU_ID, all_of(EPICON_metadata_amp_metRoot$Sample_ID_1))

EPICON_amp1 %>% dim()

# met
EPICON_met <- EPICON_MG_MT %>% 
  dplyr::select(Taxon, all_of(EPICON_metadata_amp_metRoot$ID))

EPICON_met %>% dim()

# rename
colnames(EPICON_amp1)
colnames(EPICON_met) <- colnames(EPICON_amp1)

Fungi_phy <- c("p__Ascomycota", "p__Basidiomycota", "p__Blastocladiomycota", "p__Chytridiomycota",
               "p__Mucoromycota", "p__Olpidiomycota", "p__Zoopagomycota")


EPICON_met_taxa <- EPICON_met %>% select(
  OTU_ID
) %>% separate_wider_delim(
  OTU_ID, names = c("Kingdom", "Phylum", "Class"),
  delim = ";", too_few = "debug") %>%
  select(-OTU_ID_ok, -OTU_ID_pieces, -OTU_ID_remainder) %>% 
  filter(
    Kingdom == "k__Eukaryota"
  ) %>% filter(Phylum %in% Fungi_phy)


EPICON_met %>% select(
  OTU_ID
) %>% separate_wider_delim(
  OTU_ID, names = c("Kingdom", "Phylum", "Class"),
  delim = ";", too_few = "debug") %>%
  select(-OTU_ID_ok, -OTU_ID_pieces, -OTU_ID_remainder) %>% 
  filter(
    Kingdom == "k__Eukaryota"
  ) %>% filter(Phylum == "p__NA")


# summary
EPICON_amp1_cla <- EPICON_amp1 %>% left_join(
  EPICON_fungi_taxa1 %>% select(ID, Subphylum, Class), by = c("OTU_ID" = "ID")) %>% 
  select(Class, all_of(EPICON_metadata_amp_metRoot$Sample_ID_1)) %>% 
  group_by(Class) %>% 
  summarise(across(all_of(EPICON_metadata_amp_metRoot$Sample_ID_1), ~ sum(.x))) %>% 
  filter(Class != "unclassified") %>% 
  filter(!is.na(Class))


EPICON_met_cla <- EPICON_met %>% 
  filter(OTU_ID %in% EPICON_met_taxa$OTU_ID) %>% 
  mutate(
    Class = str_split_i(OTU_ID, pattern = ";", 3)
  ) %>% filter(!is.na(Class)) %>% 
  mutate(
    Class = str_sub(Class, start = 4)
  ) %>% select(Class, all_of(EPICON_metadata_amp_metRoot$Sample_ID_1)) %>% 
  group_by(Class) %>% 
  summarise(across(all_of(EPICON_metadata_amp_metRoot$Sample_ID_1), ~ sum(.x)))


# RRN adj tab
EPICON_amp1_Class_tax <-   
  EPICON_amp1 %>% left_join(
  EPICON_fungi_taxa1 %>% 
  select(ID, Subphylum, Class), by = c("OTU_ID" = "ID")) %>% 
  select(Subphylum, Class) %>% 
  distinct_all() %>% drop_na()
colnames(EPICON_amp1_Class_tax) <- c("Phylum", "Class")

EPICON_amp1_Class_RRN_tab_ok <- 
  EPICON_amp1_Class_tax %>% 
  left_join(
    FRRN_cla_tab, by = c("Class" = "cla")
  ) %>% filter(!is.na(rDNA_m))

EPICON_amp1_Class_RRN_tab_na <- 
  EPICON_amp1_Class_tax %>% 
  left_join(
    FRRN_cla_tab, by = c("Class" = "cla")
  ) %>% filter(is.na(rDNA_m)) %>% select(-rDNA_m) %>% 
  left_join(
    FRRN_phy_tab, by = c("Phylum" = "phy")
  ) %>% drop_na()

EPICON_amp_RRN_tab <- rbind(EPICON_amp1_Class_RRN_tab_ok,
                            EPICON_amp1_Class_RRN_tab_na)

EPICON_amp_RRN_mean <- EPICON_amp_RRN_tab$rDNA_m %>% mean()

# wider to longer
# amp
EPICON_amp1_cla
EPICON_amp1_cla_l <- 
  EPICON_amp1_cla %>% pivot_longer(
    names_to = "Sample_ID", values_to = "Reads_num_amp", -Class
  )


# met
EPICON_met_cla
EPICON_met_cla_l <- 
  EPICON_met_cla %>% pivot_longer(
    names_to = "Sample_ID", values_to = "Reads_num_met", -Class
  ) %>% 
  group_by(Sample_ID) %>% mutate(
    R_p_met = Reads_num_met / sum(Reads_num_met)
  ) %>% mutate(
    MERGE_index = str_c(Class, Sample_ID, sep = "_")
  ) %>% ungroup()

# amp adj
EPICON_amp1_cla_l_adj <- 
  EPICON_amp1_cla_l %>% left_join(
    EPICON_amp_RRN_tab %>% select(Class, rDNA_m), by = "Class") %>% 
  mutate(
    rDNA_m = if_else(is.na(rDNA_m), EPICON_amp_RRN_mean, rDNA_m)
  ) %>% mutate(
    Reads_num_amp_adj = Reads_num_amp / rDNA_m
  ) %>% group_by(
    Sample_ID
  ) %>% mutate(
    R_p = Reads_num_amp / sum(Reads_num_amp),
    R_p_adj = Reads_num_amp_adj / sum(Reads_num_amp_adj)
  ) %>% mutate(
    MERGE_index = str_c(Class, Sample_ID, sep = "_")
  )

EPICON_amp1_cla_adj <- EPICON_amp1_cla_l_adj %>% select(Class, Sample_ID, Reads_num_amp_adj) %>% 
  pivot_wider(
    names_from = Sample_ID, values_from = Reads_num_amp_adj
  )

# Glomer no-adj or adj
EPICON_amp_met_cla_adj <- 
  EPICON_amp1_cla_l_adj %>%
  filter(Class == "Glomeromycetes") %>% 
  left_join(EPICON_met_cla_l %>% 
              filter(Class == "Glomeromycetes") %>% 
              select(MERGE_index, Reads_num_met, R_p_met), by = "MERGE_index")


# EPICON_amp_met_cla_adj %>% filter(R_p_met > 0.5)
# percentage
p_amp_met1 <- 
  ggplot(EPICON_amp_met_cla_adj,
       aes(R_p, R_p_met)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point() +
  geom_abline(slope = 1)

p_amp_met2 <- 
  ggplot(EPICON_amp_met_cla_adj,
       aes(R_p_adj, R_p_met)) +
  geom_smooth(method = "lm") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point() +
  geom_abline(slope = 1)

p_amp_met1 + p_amp_met2
# ggsave("bias_percentage_amp_met.pdf")

# composition
EPICON_amp1_cla_l_adj
EPICON_met_cla_l


# within and all
amp1_met_within_Class <- EPICON_amp1_cla$Class[EPICON_amp1_cla$Class %in% EPICON_met_cla$Class]
all_Class <- c(EPICON_amp1_cla$Class, EPICON_met_cla$Class) %>% unique()


amp1_unadj_total_arrange <- 
  EPICON_amp1_cla_l_adj %>% group_by(Class) %>% 
  summarise(
    total_p_sum = sum(Reads_num_amp)
  ) %>% arrange(desc(total_p_sum))
amp1_unadj_total_arrange


# check Saccharomycetes
EPICON_amp1_cla_l_adj %>% 
  filter(
    Class == "Saccharomycetes"
  )



amp1_unadj_within_Class <- 
  amp1_unadj_total_arrange %>% filter(
    Class %in% amp1_met_within_Class
  ) %>% select(Class) %>% pull()
# top
amp1_unadj_top_Class <- amp1_unadj_within_Class[1:10]


amp1_adj_total_arrange <- 
  EPICON_amp1_cla_l_adj %>% group_by(Class) %>% 
  summarise(
    total_p_sum = sum(Reads_num_amp_adj)
  ) %>% arrange(desc(total_p_sum))
amp1_adj_total_arrange

amp1_adj_within_Class <- 
  amp1_adj_total_arrange %>% filter(
    Class %in% amp1_met_within_Class
  ) %>% select(Class) %>% pull()
# top
amp1_adj_top_Class <- amp1_adj_within_Class[1:10]


met_total_arrange <- 
  EPICON_met_cla_l %>% group_by(Class) %>% 
  summarise(
    total_p_sum = sum(Reads_num_met)
  ) %>% arrange(desc(total_p_sum))
met_total_arrange

met_within_Class <- 
  met_total_arrange %>% filter(
    Class %in% amp1_met_within_Class
  ) %>% select(Class) %>% pull()
# top
met_top_Class <- met_within_Class[1:10]

# check Saccharomycetes
EPICON_met_cla_l %>% 
  filter(
    Class == "Saccharomycetes"
  ) %>% left_join(
    EPICON_metadata_amp_metRoot, by = c("Sample_ID" = "Sample_ID_1")
  ) %>% view()



amp1_met_Class_top <- c(amp1_unadj_top_Class,
                        amp1_adj_top_Class,
                        met_top_Class) %>% unique()


EPICON_comp_amp_met_part1 <- 
  EPICON_amp1_cla_l_adj %>% 
  select(Class, Sample_ID, Reads_num_amp,Reads_num_amp_adj) %>% 
  mutate(
    Class1 = if_else(Class %in% amp1_met_Class_top, Class, "Others")
  ) %>% left_join(
    EPICON_metadata_amp_metRoot %>% select(Sample_ID_1, TIME, Treatment),
    by = c("Sample_ID" = "Sample_ID_1")) %>% 
  group_by(Sample_ID) %>% 
  mutate(
    per_amp = Reads_num_amp / sum(Reads_num_amp),
    per_amp_adj = Reads_num_amp_adj / sum(Reads_num_amp_adj)
  ) %>% ungroup() %>% 
  group_by(TIME, Treatment, Class1) %>% 
  summarise(
    p_amp_unadj_m = mean(per_amp),
    p_amp_adj_m = mean(per_amp_adj)
  ) %>% 
  pivot_longer(
    names_to = "Methods", values_to = "abd", 4:5
  )

EPICON_comp_amp_met_part2 <- 
  EPICON_met_cla_l %>% select(Class, Sample_ID, Reads_num_met) %>% 
  mutate(
    Class1 = if_else(Class %in% amp1_met_Class_top, Class, "Others")
  ) %>% left_join(
    EPICON_metadata_amp_metRoot %>% select(Sample_ID_1, TIME, Treatment),
    by = c("Sample_ID" = "Sample_ID_1")) %>% 
  group_by(Sample_ID) %>% 
  mutate(
    per_met = Reads_num_met / sum(Reads_num_met)
  ) %>% ungroup() %>% 
  group_by(TIME, Treatment, Class1) %>% 
  summarise(
    per_met_m = mean(per_met)
  ) %>% 
  pivot_longer(
    names_to = "Methods", values_to = "abd", 4
  )

# EPICON_comp_amp_met_part2 %>% filter(TIME == "TP02")


EPICON_met_cla_l %>% select(Class, Sample_ID, Reads_num_met) %>% 
  mutate(
    Class1 = if_else(Class %in% amp1_met_Class_top, Class, "Others")
  ) %>% left_join(
    EPICON_metadata_amp_metRoot %>% select(Sample_ID_1, TIME, Treatment),
    by = c("Sample_ID" = "Sample_ID_1")) %>% 
  group_by(Sample_ID) %>% 
  mutate(
    per_met = Reads_num_met / sum(Reads_num_met)
  ) %>% ungroup() %>% 
  filter(
    TIME == "TP02"
  ) %>% 
  filter(
    Treatment == "Control"
  ) %>% filter(
    Class == "Saccharomycetes"
  )



# merge
EPICON_comp_amp_met <- 
  rbind(EPICON_comp_amp_met_part1, EPICON_comp_amp_met_part2)


EPICON_comp_amp_met$Class1 %>% unique()

EPICON_comp_amp_met_lev <- 
  EPICON_comp_amp_met_part2 %>% group_by(Class1) %>% 
  summarise(
    total_abd = sum(abd)
  ) %>% arrange(desc(total_abd)) %>% select(Class1) %>% pull()
EPICON_comp_amp_met_lev

# Class level
EPICON_comp_amp_met_lev1 <- c(EPICON_comp_amp_met_lev[-13], "Others")



EPICON_comp_amp_met$Class1 <- factor(EPICON_comp_amp_met$Class1, 
                                     levels = EPICON_comp_amp_met_lev1)

EPICON_comp_amp_met1 <- 
  EPICON_comp_amp_met %>% 
  mutate(
    TIME = str_replace(TIME, pattern = "TP0", "")
  ) %>% 
  mutate(
    TIME = str_replace(TIME, pattern = "TP", "")
  )

EPICON_comp_amp_met1$TIME <- factor(EPICON_comp_amp_met1$TIME,
                                    levels = c(2:17))


my_colors_soft <- c(
  "#6CA6CD",
  "#66CDAA",
  "#FF007F",
  "#EEAD0E",
  "#B0A8D1",
  "#7EC8E3",
  "#F0A080",
  "tomato",
  "violet",
  "yellowgreen",
  "peachpuff",
  "pink",
  "peru",
  "#A9A9A9"
)

fig4f_EPICON_comp_amp_met_control <- 
  ggplot(EPICON_comp_amp_met1 %>% filter(Treatment == "Control"), 
         aes(Methods, abd, fill = Class1)) +
  geom_col(position = "fill", width = 0.8) +
  scale_x_discrete(limits = c("p_amp_unadj_m", "per_met_m", "p_amp_adj_m"),
                   labels = c("uncorrected",
                              "metaT",
                              "rrn corrected")) +
  facet_grid(. ~ TIME %>% as.factor()) +
  scale_fill_manual(values = my_colors_soft) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  labs(y = "Relative abundance") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 15,face = "bold"), 
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold.italic"),
        axis.text = element_text(size = 10, face = "bold", colour = "black"),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        title = element_text(size = 15, face = "bold"),
        legend.position = "right",
        axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1)) +
  guides(fill = guide_legend(title = "", byrow = T, ncol = 1, title.position = "left", title.hjust = 0.5))

fig4f_EPICON_comp_amp_met_control


tm <- now() %>% str_split_i(pattern = " ", 1)
fig4f_pdf <- str_c("fig4f_", "EPICON_comp_amp_met_control_", tm, ".pdf", sep = "")
fig4f_jpg <- str_c("fig4f_", "EPICON_comp_amp_met_control_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4f_pdf <- str_c(fig_path, fig4f_pdf)
fig_fullpath_4f_jpg <- str_c(fig_path, fig4f_jpg)

ggsave(fig_fullpath_4f_pdf, fig4f_EPICON_comp_amp_met_control, width = 15.3, height = 4.46)
ggsave(fig_fullpath_4f_jpg, fig4f_EPICON_comp_amp_met_control, width = 15.3, height = 4.46)



# prapare class level abd table cla ~ Cla
EPICON_amp1_Cla <- data.frame(
  Class = all_Class
) %>% left_join(EPICON_amp1_cla, by = "Class") %>% 
  mutate(across(-Class, ~ replace_na(.x, 0)))
EPICON_amp1_Cla %>% dim()


EPICON_amp1_Cla_adj <- data.frame(
  Class = all_Class
) %>% left_join(EPICON_amp1_cla_adj, by = "Class") %>% 
  mutate(across(-Class, ~ replace_na(.x, 0)))
EPICON_amp1_Cla_adj %>% dim()


EPICON_met_Cla <- data.frame(
  Class = all_Class
) %>% left_join(EPICON_met_cla, by = "Class") %>% 
  mutate(across(-Class, ~ replace_na(.x, 0)))
EPICON_met_Cla %>% dim()

# CLR transform
EPICON_amp1_Cla %>% dim()
EPICON_amp1_Cla_adj %>% dim()
EPICON_met_Cla %>% dim()


EPICON_amp1_Cla_df <- EPICON_amp1_Cla %>%
  select(-Class) %>% as.data.frame()
rownames(EPICON_amp1_Cla_df) <- EPICON_amp1_Cla$Class

EPICON_amp1_Cla_adj_df <- EPICON_amp1_Cla_adj %>%
  select(-Class) %>% as.data.frame()
rownames(EPICON_amp1_Cla_adj_df) <- EPICON_amp1_Cla_adj$Class

EPICON_met_Cla_df <- EPICON_met_Cla %>%
  select(-Class) %>% as.data.frame()
rownames(EPICON_met_Cla_df) <- EPICON_met_Cla$Class


EPICON_amp1_Cla_clr <- clr(EPICON_amp1_Cla_df %>% t()) %>%
  as.data.frame() %>% 
  select(all_of(amp1_met_Class_top)) %>%
  mutate(Sample_ID = rownames(.)) %>% 
  pivot_longer(names_to = "Class", values_to = "non_adj_clr", all_of(amp1_met_Class_top)) %>% 
  mutate(MERGE_INDEX = str_c(Sample_ID, Class, sep = "_"))


# EPICON_amp1_Cla_clr %>% filter(is.na(non_adj_clr))

EPICON_amp1_Cla_adj_clr <- clr(EPICON_amp1_Cla_adj_df %>% t()) %>% 
  as.data.frame() %>% 
  select(all_of(amp1_met_Class_top)) %>%
  mutate(Sample_ID = rownames(.)) %>% 
  pivot_longer(names_to = "Class", values_to = "adj_clr", all_of(amp1_met_Class_top)) %>% 
  mutate(MERGE_INDEX = str_c(Sample_ID, Class, sep = "_")) %>% 
  select(MERGE_INDEX, adj_clr)

# EPICON_amp1_Cla_adj_clr %>% filter(is.na(adj_clr))


EPICON_met_Cla_clr <- clr(EPICON_met_Cla_df %>% t()) %>%
  as.data.frame() %>% 
  select(all_of(amp1_met_Class_top)) %>%
  mutate(Sample_ID = rownames(.)) %>% 
  pivot_longer(names_to = "Class", values_to = "met_clr", all_of(amp1_met_Class_top)) %>% 
  mutate(MERGE_INDEX = str_c(Sample_ID, Class, sep = "_")) %>% 
  select(MERGE_INDEX, met_clr)



EPICON_amp_met_bias <- 
  EPICON_amp1_Cla_clr %>% 
  left_join(EPICON_amp1_Cla_adj_clr) %>% 
  left_join(EPICON_met_Cla_clr) %>% 
  mutate(
    amp_non_adj_bias = non_adj_clr - met_clr,
    amp_adj_bias = adj_clr - met_clr
  )

EPICON_amp_met_bias1 <- 
  EPICON_amp_met_bias %>% 
  select(Sample_ID, Class, amp_non_adj_bias, amp_adj_bias) %>% 
  group_by(Class) %>% 
  summarise(
    non_adj_bias_m = mean(amp_non_adj_bias),
    adj_bias_m = mean(amp_adj_bias)
  ) %>% pivot_longer(
    names_to = "group", values_to = "bias_m", -Class
  )
EPICON_amp_met_bias1$bias_m %>% range()

EPICON_amp_met_bias_radar <- 
  EPICON_amp_met_bias1 %>% 
  pivot_wider(
    names_from = Class, values_from = bias_m
  )

EPICON_amp_met_bias_radar$group <- c("RRN uncorrected", "RRN corrected")
EPICON_amp_met_bias_radar$group <- factor(EPICON_amp_met_bias_radar$group,
                                          levels = c("RRN uncorrected", "RRN corrected"))


amp1_met_Class_main <- c(amp1_unadj_top_Class[1:5],
                        amp1_adj_top_Class[1:5],
                        met_top_Class[1:5]) %>% unique()




amp1_met_Class_main_sort <- EPICON_amp_met_bias1 %>% filter(
  Class %in% amp1_met_Class_main
) %>% filter(group == "adj_bias_m") %>% arrange(desc(bias_m)) %>% select(Class) %>% pull()



EPICON_amp_met_bias_radar_main_Class <- EPICON_amp_met_bias_radar %>% 
  select(group, all_of(amp1_met_Class_main_sort))

EPICON_amp_met_bias_radar_main_Class %>% select(-group) %>% range()


fig4h_EPICON_amp_met_radar_mainC <- 
  ggradar(EPICON_amp_met_bias_radar_main_Class,
          values.radar = c("-4.3", "0", "0.5"),
          grid.label.size = 10,
          #plot.extent.x.sf = 2,
          #plot.extent.y.sf = 2,
          grid.min = -4.3,
          grid.mid = 0, gridline.mid.linetype = 2, gridline.mid.colour = "black",
          grid.max = 0.5,
          group.line.width = 1, 
          group.point.size = 5,
          group.colours = c("#4455CC", "#FF007F"),
          axis.label.size = 6,
          # axis.labels = dPCR_rlts1$Taxa,
          axis.label.offset = 1.1,
          # plot.title = "Bias between metabarcoding and metatranscriptom",
          legend.title = "Method",
          background.circle.colour = "grey95") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, vjust = 0.5)
  )
fig4h_EPICON_amp_met_radar_mainC

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4h_pdf <- str_c("fig4h_", "EPICON_amp_met_radar_mainC_", tm, ".pdf", sep = "")
fig4h_jpg <- str_c("fig4h_", "EPICON_amp_met_radar_mainC_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4h_pdf <- str_c(fig_path, fig4h_pdf)
fig_fullpath_4h_jpg <- str_c(fig_path, fig4h_jpg)

ggsave(fig_fullpath_4h_pdf, fig4h_EPICON_amp_met_radar_mainC, width = 7.96, height = 7.21)
ggsave(fig_fullpath_4h_jpg, fig4h_EPICON_amp_met_radar_mainC, width = 7.96, height = 7.21)

# Saving 6.69 x 5.75 in image ?


EPICON_amp_met_bias_Glo <- 
  EPICON_amp_met_bias %>% filter(Class == "Glomeromycetes") %>% 
  select(Sample_ID, amp_non_adj_bias, amp_adj_bias) %>% 
  left_join(
    EPICON_metadata_amp_metRoot %>% select(Sample_ID_1, TIME, Treatment), 
    by = c("Sample_ID" = "Sample_ID_1") 
  ) %>% 
  pivot_longer(
    names_to = "Methods", values_to = "bias", 2:3
  )


EPICON_amp_met_bias_Glo1 <- 
  EPICON_amp_met_bias_Glo %>% 
  mutate(
    TIME = str_replace(TIME, pattern = "TP0", "")
  ) %>% 
  mutate(
    TIME = str_replace(TIME, pattern = "TP", "")
  )


EPICON_amp_met_bias_Glo1$TIME <- 
  factor(EPICON_amp_met_bias_Glo1$TIME, levels = c(2:17))




EPICON_amp_met_bias_Glo1$Methods <- factor(EPICON_amp_met_bias_Glo$Methods, 
                                          levels = c("amp_non_adj_bias", "amp_adj_bias"))


EPICON_amp_met_bias_Glo1$TIME %>% table()

EPICON_bias_AMF <- EPICON_amp_met_bias_Glo1 %>%
  filter(Treatment == "Control")

EPICON_bias_AMF$TIME %>% table()


EPICON_bias_AMF$TIME <- EPICON_bias_AMF$TIME %>% as.vector %>% as.numeric()

EPICON_bias_AMF$TIME %>% table()

EPICON_bias_AMF_t.test <- 
  EPICON_bias_AMF %>% select(Sample_ID, Methods, bias) %>%
  pivot_wider(names_from = Methods, values_from = bias)


bias_AMF_t.test <- t.test(EPICON_bias_AMF_t.test$amp_non_adj_bias, EPICON_bias_AMF_t.test$amp_adj_bias, paired = T)

bias_AMF_t.test$p.value

# help("geom_half_violin")

EPICON_bias_AMF %>% group_by(Methods) %>% 
  summarise(
    bias_m = mean(bias)
  ) %>% 
  mutate(
    fold_m = exp(bias_m)
  ) %>% 
  mutate(
    fold_m1 = 1/fold_m
  )


###### AMF dev ######
fig4i_amf_dev <- 
  ggplot(EPICON_bias_AMF) +
  geom_half_violin(aes(Methods, bias), side = c("l", "r"), colour = NA, fill = "grey90") +
  geom_path(aes(Methods, bias, group = Sample_ID), colour = "grey90") +
  geom_point(aes(Methods, bias, colour = TIME), size = 2.8, alpha = 0.6) +
  scale_x_discrete(limits = c("amp_non_adj_bias", "amp_adj_bias"),
                   labels = c("rrn\nuncorrected", "rrn\ncorrected")) +
  scale_colour_gradient2(breaks = seq(2, 17, 1),
                         labels = seq(2, 17, 1),
                         # low = "#4455CC", high = "#FF007F",
                         low = "blue", high = "red",
                         mid = "grey90", midpoint = 9.5) +
  #scale_colour_manual(limits = c("amp_non_adj_bias", "amp_adj_bias"),
  #                    labels = c("rrn uncorrected", "rrn corrected"),
  #                    values = c("#4455CC", "#FF007F")) +
  geom_hline(yintercept = 0, linetype = 2, colour = "black") +
  annotate(geom = "text", x = 1.5, y = 3.3, 
           label = expression("t" == "-160.99;" ~~ "df" == "1"), size = 3.8) +
  annotate(geom = "text", x = 1.5, y = 2.8, 
           label = expression(italic(p) == "2.99e-102"), size = 3.8) +
  #scale_alpha(breaks = seq(2, 17, 1),
  #            labels = seq(2, 17, 1),
  #            range = c(0.05, 0.7)) +
  guides(colour = guide_legend(title = "Week")) +
  labs(y = "Devation of Glomeromycetes") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 15,face = "bold"), 
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 15, face = "bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right")
fig4i_amf_dev


tm <- now() %>% str_split_i(pattern = " ", 1)
fig4i_pdf <- str_c("fig4i_", "EPICON_amf_dev_", tm, ".pdf", sep = "")
fig4i_jpg <- str_c("fig4i_", "EPICON_amf_dev_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4i_pdf <- str_c(fig_path, fig4i_pdf)
fig_fullpath_4i_jpg <- str_c(fig_path, fig4i_jpg)

ggsave(fig_fullpath_4i_pdf, fig4i_amf_dev, width = 3.08, height = 5.52)
ggsave(fig_fullpath_4i_jpg, fig4i_amf_dev, width = 3.08, height = 5.52)


###### diff methods PCoA #####
# data
EPICON_amp1_Cla %>% dim()
EPICON_amp1_Cla_adj %>% dim()
EPICON_met_Cla

# sub data
EPICON_amp1_Cla_trt <- EPICON_amp1_Cla %>% select(Class, all_of(EPICON_metadata_amp_metRoot %>%
                                                                  filter(Treatment == "Control") %>%
                                                                  select(Sample_ID_1) %>% pull()))
EPICON_amp1_Cla_adj_trt <- EPICON_amp1_Cla_adj %>% select(Class, all_of(EPICON_metadata_amp_metRoot %>%
                                                                          filter(Treatment == "Control") %>%
                                                                          select(Sample_ID_1) %>% pull()))

colnames(EPICON_amp1_Cla_trt)[-1] <- str_c(colnames(EPICON_amp1_Cla_trt)[-1], "Unc", sep = "_")
colnames(EPICON_amp1_Cla_adj_trt)[-1] <- str_c(colnames(EPICON_amp1_Cla_adj_trt)[-1], "C", sep = "_")


EPICON_met_Cla_trt <- EPICON_met_Cla %>% select(Class, all_of(EPICON_metadata_amp_metRoot %>%
                                                                filter(Treatment == "Control") %>%
                                                                select(Sample_ID_1) %>% pull()))

colnames(EPICON_met_Cla_trt)[-1] <- str_c(colnames(EPICON_met_Cla_trt)[-1], "metaT", sep = "_")


EPICON_UCT <- 
  EPICON_amp1_Cla_trt %>% 
  left_join(EPICON_amp1_Cla_adj_trt) %>% 
  left_join(EPICON_met_Cla_trt)


EPICON_UCT_df <- EPICON_UCT %>%
  select(-Class) %>% as.data.frame()
rownames(EPICON_UCT_df) <- EPICON_UCT$Class


EPICON_UCT_df_t <- EPICON_UCT_df %>% t() %>% as.data.frame()

# EPICON_pcoa_UCT_clr <- clr() %>%
#   as.data.frame()


# metadata
EPICON_UCT_METADATA <- data.frame(
  Sample_ID = rownames(EPICON_UCT_df_t)
) %>% mutate(
  Methods = str_split_i(Sample_ID, pattern = "_", 2)
) %>% mutate(
  sample_id = str_split_i(Sample_ID, pattern = "_", 1)
) %>% left_join(
  EPICON_metadata_amp_metRoot %>% select(Sample_ID_1, TIME), by = c("sample_id" = "Sample_ID_1")
) %>% mutate(
  TIME = str_replace(TIME, pattern = "TP0", "")
) %>% mutate(
    TIME = str_replace(TIME, pattern = "TP", "")
  )

EPICON_UCT_METADATA$Methods <- as.factor(EPICON_UCT_METADATA$Methods)
EPICON_UCT_METADATA$TIME <- EPICON_UCT_METADATA$TIME %>% as.vector() %>% as.numeric()

rownames(EPICON_UCT_METADATA) <- EPICON_UCT_METADATA$Sample_ID

# PCoA
parse_PcoA <- function(DF, METADATA, PRE_METHOD) {
  
  # calcu euclidean dist
  DF_rel <- decostand(DF, method = PRE_METHOD)
  
  DF_dist <- vegdist(DF_rel, method = "bray")
  
  # PCoA analysis
  TMP_pcoa_result <- wcmdscale(DF_dist, k = 2, eig = TRUE)
  
  # tmp
  # return(TMP_pcoa_result)
  
  TMP_points <- as.data.frame(TMP_pcoa_result$points)
  
  colnames(TMP_points) <- c("PCoA1", "PCoA2")
  
  TMP_points1 <- TMP_points %>% mutate(
    Sample_ID = rownames(.)
  ) %>% left_join(
    METADATA, by = "Sample_ID")
  
  # eig
  all_eig <- TMP_pcoa_result$eig
  pos_eig <- all_eig[all_eig > 0]
  
  TMP_var_explained <- round(pos_eig / sum(pos_eig) * 100, 2)[1:2]
  
  
  # rlts
  pcoa_rlts <- list(
    main_df = TMP_points1,
    dist = DF_dist,
    axis_eig = TMP_var_explained
  )
  
  return(pcoa_rlts)
  
}

# UCT pcoa
EPICON_UCT_pcoa <- parse_PcoA(DF = EPICON_UCT_df_t, METADATA = EPICON_UCT_METADATA, PRE_METHOD = "total")
EPICON_UCT_pcoa_adonis2 <- 
  adonis2(EPICON_UCT_pcoa$dist ~ Methods, data = EPICON_UCT_pcoa$main_df %>% select(Sample_ID, Methods),
          permutations = 999)



# UT pcoa
EPICON_UCT_df_t_UT <- EPICON_UCT_df_t %>% mutate(Sample_ID = rownames(.)) %>% 
  filter(str_detect(Sample_ID, pattern = "_Unc") | str_detect(Sample_ID, pattern = "_metaT")) %>% 
  select(-Sample_ID)

EPICON_UT_pcoa <- parse_PcoA(DF = EPICON_UCT_df_t_UT , METADATA = EPICON_UCT_METADATA, PRE_METHOD = "total")
adonis2(EPICON_UT_pcoa$dist ~ Methods, data = EPICON_UT_pcoa$main_df %>% select(Sample_ID, Methods, TIME))


# CT pcoa
EPICON_UCT_df_t_CT <- EPICON_UCT_df_t %>% mutate(Sample_ID = rownames(.)) %>% 
  filter(str_detect(Sample_ID, pattern = "_C") | str_detect(Sample_ID, pattern = "_metaT")) %>% 
  select(-Sample_ID)

EPICON_CT_pcoa <- parse_PcoA(DF = EPICON_UCT_df_t_CT, METADATA = EPICON_UCT_METADATA, PRE_METHOD = "total")
adonis2(EPICON_CT_pcoa$dist ~ Methods, data = EPICON_CT_pcoa$main_df %>% select(Sample_ID, Methods))

# EPICON_UCT_pcoa$main_df$TIME <- EPICON_UCT_pcoa$main_df$TIME %>% as.vector() %>% as.numeric()

# EPICON_UCT_pcoa$main_df$TIME %>% str()

# EPICON_UCT_pcoa$main_df %>% str()


fig4g_EPICON_UCT_pcoa_plot <- 
  ggplot() +
  geom_path(
    data = EPICON_UCT_pcoa$main_df,
    aes(PCoA1, PCoA2, group = sample_id),
    colour = "grey90"
  ) +
  geom_point(
    data = EPICON_UCT_pcoa$main_df, 
    aes(PCoA1, PCoA2, colour = Methods, alpha = TIME),
    size = 3.2) +
  annotate(geom = "text", x = 0.27, y = 0.25, 
           label = expression("un corrected vs metaT:" ~~ R^2 == 0.129 ~~ italic(p) == 0.001),
           size = 3.8, parse = T, hjust = 0) +
  annotate(geom = "text", x = 0.27, y = 0.21, 
           label = expression("rrn corrected vs metaT:" ~~ R^2 == 0.0918 ~~ italic(p) == 0.001),
           size = 3.8, parse = T, hjust = 0) +
  scale_alpha(breaks = c(seq(2, 17, 3)),
              labels = c(seq(2, 17, 3)),
              range = c(0.2, 0.7)) +
  scale_colour_manual(
    limits = c("Unc", "C", "metaT"),
    labels = c("uncorrected", "rrn corrected", "metaT"),
    values = c("#4455CC", "#FF007F", "#F5B700")
  ) +
  guides(alpha = guide_legend(title = "Week", byrow = T, nrow = 1, direction = "horizontal"),
         colour = guide_legend(title = "Method", byrow = T, nrow = 1, direction = "horizontal")) +
  labs(x = str_c("PCoA1 (", EPICON_UCT_pcoa$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", EPICON_UCT_pcoa$axis_eig[2], "%)")) +
  # scale_x_continuous(limits = c(-0.2, 0.3)) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    ),
    legend.title = element_text(size = 12, colour = "black", face = "bold"),
    legend.text = element_text(size = 9, colour = "black", face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.13),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.spacing.y = unit(0, "cm"),
    legend.margin = margin(2, 2, 2, 2, unit = "pt")
  )
fig4g_EPICON_UCT_pcoa_plot


tm <- now() %>% str_split_i(pattern = " ", 1)
fig4g_pdf <- str_c("fig4g_", "EPICON_UCT_pcoa_", tm, ".pdf", sep = "")
fig4g_jpg <- str_c("fig4g_", "EPICON_UCT_pcoa_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4g_pdf <- str_c(fig_path, fig4g_pdf)
fig_fullpath_4g_jpg <- str_c(fig_path, fig4g_jpg)

ggsave(fig_fullpath_4g_pdf, fig4g_EPICON_UCT_pcoa_plot, width = 7.55, height = 4.65)
ggsave(fig_fullpath_4g_jpg, fig4g_EPICON_UCT_pcoa_plot, width = 7.55, height = 4.65)


fig4_design <- c("AAAAAAAAAABBBBBBCC")

fig4ghi <- EPICON_UCT_pcoa_plot + plot_spacer() + fig4_amf_dev + plot_layout(design = fig4_design)
fig4ghi

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4ghi_pdf <- str_c("fig4g_", "ghi_", tm, ".pdf", sep = "")
fig4ghi_jpg <- str_c("fig4g_", "ghi_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_4ghi_pdf <- str_c(fig_path, fig4ghi_pdf)
fig_fullpath_4ghi_jpg <- str_c(fig_path, fig4ghi_jpg)

ggsave(fig_fullpath_4ghi_pdf, fig4ghi, width = 15.3, height = 5.8)
ggsave(fig_fullpath_4ghi_jpg, fig4ghi, width = 15.3, height = 5.8)


ggsave("EPICON_PCoA_AMF_dev.pdf", fig4ghi, width = 15.3, height = 5.8)
# Saving 15.3 x 5.41 in image



















###### check Sacc and Muco ######

# Sacc
EPICON_UCT_METADATA


EPICON_UCT_df_Class <- EPICON_UCT_df %>% mutate(
  Class = rownames(.), .before = TP02R12_Unc
) %>% select(Class, all_of(
  EPICON_UCT_METADATA %>% filter(Methods != "C") %>% select(Sample_ID) %>% pull()
))


EPICON_UCT_df_Class$Class %>% unique()


EPICON_UCT_df_Class_Sacc <- EPICON_UCT_df_Class %>% 
  mutate(
    Class1 = if_else(Class == "Saccharomycetes", Class, "Others")
  ) %>% select(-Class)

EPICON_UCT_df_Class_Sacc1 <- 
  EPICON_UCT_df_Class_Sacc %>% 
  pivot_longer(
    names_to = "Sample_ID",values_to = "abd", -Class1
  ) %>% group_by(
    Class1, Sample_ID
  ) %>% 
  summarise(
    abd_sum = sum(abd)
  ) %>% ungroup() %>% group_by(Sample_ID) %>% 
  mutate(
    abd_per = abd_sum / sum(abd_sum)
  ) %>% 
  mutate(
    Methods = str_split_i(Sample_ID, pattern = "_", 2)
  )
  
EPICON_UCT_df_Class_Sacc1


ggplot(EPICON_UCT_df_Class_Sacc1 %>% filter(Class1 == "Saccharomycetes"), 
       aes(Methods, abd_per)) +
  geom_point() +
  theme_bw()



ggplot(EPICON_UCT_df_Class_Sacc1 %>% filter(Class1 == "Saccharomycetes") %>% 
         filter(Methods == "metaT"), 
       aes(abd_per)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw()


EPICON_UCT_df_Class_Sacc1 %>% filter(Class1 == "Saccharomycetes") %>% 
  filter(Methods == "metaT") %>% arrange(desc(abd_per)) %>% view()


metaT_Sacc <- 
  EPICON_UCT_df_Class_Sacc1 %>%
  filter(Methods == "metaT") %>% 
  left_join(EPICON_UCT_METADATA %>% select(Sample_ID, TIME))

metaT_Sacc$Class1 <- factor(metaT_Sacc$Class1, levels = c("Saccharomycetes", "Others"))

Unc_Sacc <- 
  EPICON_UCT_df_Class_Sacc1 %>%
  filter(Methods == "Unc") %>% 
  left_join(EPICON_UCT_METADATA %>% select(Sample_ID, TIME))

Unc_Sacc$Class1 <- factor(Unc_Sacc$Class1, levels = c("Saccharomycetes", "Others"))

p_metaT_Sacc <- 
  ggplot(metaT_Sacc, aes(Sample_ID, abd_per, fill = Class1)) +
  geom_col() +
  scale_fill_manual(values = c("red", "grey90")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    )
  )


p_Unc_Sacc <- 
  ggplot(Unc_Sacc, aes(Sample_ID, abd_per, fill = Class1)) +
  geom_col() +
  scale_fill_manual(values = c("red", "grey90")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    )
  )

p_Unc_Sacc / p_metaT_Sacc


# Mucc

EPICON_UCT_df_Class$Class %>% unique()


EPICON_UCT_df_Class_Muco <- EPICON_UCT_df_Class %>% 
  mutate(
    Class1 = if_else(Class == "Mucoromycetes", Class, "Others")
  ) %>% select(-Class)

EPICON_UCT_df_Class_Muco1 <- 
  EPICON_UCT_df_Class_Muco %>% 
  pivot_longer(
    names_to = "Sample_ID",values_to = "abd", -Class1
  ) %>% group_by(
    Class1, Sample_ID
  ) %>% 
  summarise(
    abd_sum = sum(abd)
  ) %>% ungroup() %>% group_by(Sample_ID) %>% 
  mutate(
    abd_per = abd_sum / sum(abd_sum)
  ) %>% 
  mutate(
    Methods = str_split_i(Sample_ID, pattern = "_", 2)
  )

EPICON_UCT_df_Class_Muco1


ggplot(EPICON_UCT_df_Class_Muco1 %>% filter(Class1 == "Mucoromycetes"), 
       aes(Methods, abd_per)) +
  geom_point() +
  theme_bw()



ggplot(EPICON_UCT_df_Class_Muco1 %>% filter(Class1 == "Mucoromycetes") %>% 
         filter(Methods == "metaT"), 
       aes(abd_per)) +
  geom_histogram(binwidth = 0.01) +
  theme_bw()


EPICON_UCT_df_Class_Muco1 %>% filter(Class1 == "Mucoromycetes") %>% 
  filter(Methods == "metaT") %>% arrange(desc(abd_per)) %>% view()


metaT_Muco <- 
  EPICON_UCT_df_Class_Muco1 %>%
  filter(Methods == "metaT") %>% 
  left_join(EPICON_UCT_METADATA %>% select(Sample_ID, TIME))


Unc_Muco <- 
  EPICON_UCT_df_Class_Muco1 %>%
  filter(Methods == "Unc") %>% 
  left_join(EPICON_UCT_METADATA %>% select(Sample_ID, TIME))

p_metaT_Muco <- 
  ggplot(metaT_Muco, aes(Sample_ID, abd_per, fill = Class1)) +
  geom_col() +
  scale_fill_manual(values = c("red", "grey90")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    )
  )


p_Unc_Muco <- 
  ggplot(Unc_Muco, aes(Sample_ID, abd_per, fill = Class1)) +
  geom_col() +
  scale_fill_manual(values = c("red", "grey90")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    )
  )

p_Unc_Muco / p_metaT_Muco



# fig5g_amf_box1 <- 
#   ggplot(EPICON_amp_met_bias_Glo1 %>% filter(Treatment == "Control"), aes(TIME, bias, colour = Methods)) +
#   geom_jitter(position = position_dodge(width = 1)) +
#   geom_boxplot(position = position_dodge(width = 1)) +
#   facet_grid(. ~ TIME, scales = "free") +
#   scale_colour_manual(limits = c("amp_non_adj_bias", "amp_adj_bias"),
#                       labels = c("rrn uncorrected", "rrn corrected"),
#                       values = c("#4455CC", "#FF007F")) +
#   geom_hline(yintercept = 0, linetype = 3) +
#   labs(y = "Devation") +
#   guides(colour = guide_legend(title = "", byrow = T, nrow = 1)) +
#   theme_bw() +
#   theme(axis.title.x = element_blank(), 
#         strip.text = element_text(size = 15,face = "bold"), 
#         panel.spacing = unit(0, "lines"),
#         legend.title = element_blank(),
#         legend.background = element_rect(fill = alpha("grey", 0.3), color = NA),
#         legend.text = element_text(colour = "black", size = 10, face = "bold"),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
#         axis.title = element_text(size = 15, face = "bold", colour = "black"),
#         title = element_text(size = 15, face = "bold"),
#         legend.position = "inside",
#         legend.position.inside = c(0.5, 0.1))
# fig5g_amf_box1
# 
# 
# # ggsave("p_amf_boxplot_test.pdf")
# 
# fig5fg_EPICON_comp_dev <- fig5f_EPICON_comp_amp_met_control / fig5g_amf_box1
# ggsave("fig5fg_EPICON_comp_dev1.pdf", fig5fg_EPICON_comp_dev,
#        width = 15.2, height = 8.07)




###### mock commu PCoA #####

all_table


filtered_ASVs <- 
  rbind(
    ASVs_S01_yes,
    ASVs_S02_yes,
    ASVs_S03_yes,
    ASVs_S04_yes,
    ASVs_S05_yes,
    ASVs_S07_yes,
    ASVs_S08_all,
    ASVs_S10_yes
  )


all_table1 <- all_table %>% filter(
  OTU_ID %in% filtered_ASVs$OTU_ID
)

all_table1_df <- all_table1 %>% select(-OTU_ID) %>% as.data.frame()
rownames(all_table1_df) <- all_table1$OTU_ID

all_table1_df %>% apply(sum, MARGIN = 2) %>% min()


all_table1_df_r <- rrarefy(all_table1_df %>% t(), 6927)



all_table_metadata <- mock_commu_metadata %>% 
  select(sample_id, strain_list) %>% 
  mutate(
    group = str_sub(sample_id, 1, 1)
  )

colnames(all_table_metadata)[1] <- "Sample_ID"

# PCoA
all_table_pcoa <- 
  parse_PcoA(DF = all_table1_df_r, METADATA = all_table_metadata, PRE_METHOD = "hellinger")

all_table1_df_r %>% view()

mock_commu_pcoa_plot <- 
  ggplot() +
  geom_point(
    data = all_table_pcoa$main_df, 
    aes(PCoA1, PCoA2, colour = group),
    size = 3.2) +
  scale_colour_manual(
    limits = c("D", "S", "M"),
    labels = c("Equal_D", "Equal_S", "Equal_M"),
    values = c("#4455CC", "#FF007F", "#F5B700")
  ) +
  labs(x = str_c("PCoA1 (", EPICON_UCT_pcoa$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", EPICON_UCT_pcoa$axis_eig[2], "%)")) +
  # scale_x_continuous(limits = c(-0.2, 0.3)) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    )
  )
mock_commu_pcoa_plot
# ggsave("EPICON_UCT_pcoa.pdf")

mock_commu_vrlts_table <- member_percentage_vsearch_combAMF %>% 
  select(sample_id, strain_id, n_1) %>% 
  pivot_wider(
    names_from = strain_id, values_from = n_1
  )


mock_commu_vrlts_df <- mock_commu_vrlts_table %>% 
  select(-sample_id) %>% as.data.frame()

rownames(mock_commu_vrlts_df) <- mock_commu_vrlts_table$sample_id


mock_commu_vrlts_df %>% t() %>% apply(MARGIN = 2, sum) %>% min()


mock_commu_vrlts_df_r <- rrarefy(mock_commu_vrlts_df, 11424)

mock_commu_vrlts_df_r

# PCoA
mock_commu_vrlts_pcoa <- 
  parse_PcoA(DF = mock_commu_vrlts_df_r, METADATA = all_table_metadata, PRE_METHOD = "hellinger")

library(ggsci)

mock_commu_vrlts_pcoa_plot0 <- 
  ggplot() +
  geom_point(
    data = mock_commu_vrlts_pcoa$main_df, 
    aes(PCoA1, PCoA2, colour = group),
    size = 3.2) +
  scale_colour_d3() +
  labs(x = str_c("PCoA1 (", mock_commu_vrlts_pcoa$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", mock_commu_vrlts_pcoa$axis_eig[2], "%)")) +
  # scale_x_continuous(limits = c(-0.2, 0.3)) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    )
  )
mock_commu_vrlts_pcoa_plot0




mock_commu_vrlts_pcoa_plot <- 
  ggplot() +
  geom_point(
    data = mock_commu_vrlts_pcoa$main_df, 
    aes(PCoA1, PCoA2, colour = strain_list, shape = group),
    size = 3.2) +
  scale_colour_manual(values = all_grp_pal) +
  labs(x = str_c("PCoA1 (", mock_commu_vrlts_pcoa$axis_eig[1], "%)"),
       y = str_c("PCoA2 (", mock_commu_vrlts_pcoa$axis_eig[2], "%)")) +
  # scale_x_continuous(limits = c(-0.2, 0.3)) +
  theme_bw() +
  theme(
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 13),
    axis.title.y = element_text(face = "bold",
                                size = 13),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.title = element_text(
      colour = "darkgreen", size = 18, face = "bold.italic", hjust = 0.5, vjust = 0.5
    )
  )
mock_commu_vrlts_pcoa_plot


all_grp_pal <- c(pal_d3(palette = "category20b")(10),
                 pal_d3(palette = "category20c")(10),
                 pal_d3(palette = "category20")(10))



##### save tmp RData #####
FRRN_cla_tab
FRRN_phy_tab
FRRN_rlt_taxa_spl_grp
save(FRRN_cla_tab, FRRN_phy_tab, FRRN_rlt_taxa_spl_grp, EPICON_fungal_otutab1, EPICON_fungal_otutab1_1, EPICON_fungi_taxa1, file = "./6.RData/FRRN_fig4_tmp.RData")



