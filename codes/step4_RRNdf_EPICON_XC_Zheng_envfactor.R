


##### rDNA copy number of fungi in EPICON #####
# 
# by Qiushi-Li, IM-CAS, 2026-08-20

library(tidyverse)
library(vegan)
library(broom)
library(readxl)
library(patchwork)
library(gghalves)
library(agricolae)
library(ggsci)

# load("FRRN_20241211.Rdata")

# EPICON
load("./2.database/EPICON.data.preparation.RC.bNTI.ted.2019.04.19.Rdata")


# Annotated
EPICON_fungi_taxa <- read.csv("./2.database/All_otu_blast.csv", header = T)
EPICON_fungi_taxa %>% colnames()

# view(EPICON_fungi_taxa)


EPICON_fungi_taxa1 <- EPICON_fungi_taxa %>% 
  mutate(Class = str_sub(Class, start = 4),
         Order = str_sub(Order, start = 4),
         Family = str_sub(Family, start = 4),
         Genus = str_sub(Genus, start = 4),
         Subphylum = str_sub(Subphylum, start = 4)) %>% 
  filter(Kingdom == "k__Fungi")


# rDNA table
FRRN_taxa_FG

# explore --------------
dim(fung0)
fung0$Fungi %>% table()

EPICON_fung1 <- fung0 %>% filter(Fungi == "Fungi")
dim(EPICON_fung1)
view(EPICON_fung1)

EPICON_fungal_otutab <- EPICON_fung1 %>% select(1:1251)
dim(EPICON_fungal_otutab)

EPICON_fungal_otuExtra <- EPICON_fung1 %>% select(1252:1302)

EPICON_fungal_otuExtra1 <- 
  EPICON_fungal_otuExtra %>% select(OTU, ID)
EPICON_fungal_otuExtra1

# spc tab
FRRN_spc_tab <- FRRN_taxa_FG %>% mutate(
  spc = str_c(str_split_i(Name, pattern = " ", 1), str_split_i(Name, pattern = " ", 2), sep = " ")
  ) %>% group_by(spc) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_spc_tab

# gen tab
FRRN_gen_tab <- FRRN_taxa_FG %>% select(gen, Both_ITS_LSU) %>%
  group_by(gen) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_gen_tab

# fam tab
FRRN_fam_tab <- FRRN_taxa_FG %>% select(fam, Both_ITS_LSU) %>%
  group_by(fam) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_fam_tab

# ord tab
FRRN_ord_tab <- FRRN_taxa_FG %>% select(ord, Both_ITS_LSU) %>%
  group_by(ord) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_ord_tab

# cla tab
FRRN_cla_tab <- FRRN_taxa_FG %>% select(cla, Both_ITS_LSU) %>%
  group_by(cla) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_cla_tab

# phy tab
FRRN_phy_tab <- FRRN_taxa_FG %>% select(phy, Both_ITS_LSU) %>%
  group_by(phy) %>% summarise(rDNA_m = mean(Both_ITS_LSU))
FRRN_phy_tab


################### ---------

# FGanno -----------------------
# EPICON_fungi_PL %>% select(Order) %>% view()

FunGuild_database1
FG_index_df

# saveRDS(FG_index_df, "FG_index_df.RData")

FGanno_gen <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[1]) %>% select(-Taxon_Level)
colnames(FGanno_gen)[1] <- FG_index_df$Tax_Ind[1]

FGanno_EPICON_tmp <- EPICON_fungi_taxa1 %>% 
  left_join(FGanno_gen, by = c("Genus" = "gen"))

FGanno_EPICON_na <- FGanno_EPICON_tmp %>% filter(is.na(Trophic_Mode))
FGanno_EPICON_na
FGanno_EPICON_gen <- FGanno_EPICON_tmp %>% filter(!is.na(Trophic_Mode))
FGanno_EPICON_gen

FGanno_fam <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[2]) %>% select(-Taxon_Level)
colnames(FGanno_fam)[1] <- FG_index_df$Tax_Ind[2]

FGanno_EPICON_fam <- FGanno_EPICON_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FGanno_fam, by = c("Family" = "fam"))

FGanno_EPICON <- rbind(FGanno_EPICON_gen, FGanno_EPICON_fam)

FGanno_EPICON_FFF <- 
  FGanno_EPICON %>% 
  mutate(Guild1 = str_extract(Guild, "\\|([^\\|]+)\\|"), .before = Growth_Morphology) %>% 
  mutate(Guild1 = if_else(!is.na(Guild1), str_sub(Guild1, start = 2, end = -2), Guild)) %>% 
  mutate(Guild1 = if_else(is.na(Guild1), "Unknow", Guild1))
FGanno_EPICON_FFF

# view(EPICON_fungi_taxa1)

# frDNA --------------------------
EPICON_frDNA_gen <- EPICON_fungi_taxa1 %>% 
  left_join(FRRN_gen_tab, by = c("Genus" = "gen")) %>% 
  select(ID, Genus, rDNA_m) %>% filter(!is.na(rDNA_m)) %>%
  mutate(tax_lev = "Gen", .after = ID)
#EPICON_frDNA_gen
colnames(EPICON_frDNA_gen)[3] <- "tax"

EPICON_frDNA_anno <- EPICON_frDNA_gen

EPICON_frDNA_fam <- EPICON_fungi_taxa1 %>% 
  filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(FRRN_fam_tab, by = c("Family" = "fam")) %>% 
  select(ID, Family, rDNA_m) %>% filter(!is.na(rDNA_m)) %>%
  mutate(tax_lev = "Fam", .after = ID)
#EPICON_frDNA_fam
colnames(EPICON_frDNA_fam)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_fam)

EPICON_frDNA_ord <- EPICON_fungi_taxa1 %>% 
  filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(FRRN_ord_tab, by = c("Order" = "ord")) %>% 
  select(ID, Order, rDNA_m) %>% filter(!is.na(rDNA_m)) %>%
  mutate(tax_lev = "Ord", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_ord)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_ord)


EPICON_frDNA_cla <- EPICON_fungi_taxa1 %>% 
  filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(FRRN_cla_tab, by = c("Class" = "cla")) %>% 
  select(ID, Class, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% 
  mutate(tax_lev = "Cla", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_cla)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_cla)


EPICON_frDNA_phy <- EPICON_fungi_taxa1 %>%
  filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(FRRN_phy_tab, by = c("Subphylum" = "phy")) %>% 
  select(ID, Subphylum, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% 
  mutate(tax_lev = "Phy", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_phy)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_phy)
EPICON_frDNA_anno

# done ...

###### EPICON supp ######

# otutab?
EPICON_fungal_otutab
colnames(EPICON_fungal_otutab)

# env
env

# rDNA anno
EPICON_frDNA_anno

# ok~ start ~
EPICON_fungal_otutab1a <- EPICON_fungal_otutab %>% select(env$aa)
EPICON_fungal_otutab1a


EPICON_table_sum <- 
  EPICON_fungal_otutab1a %>% apply(MARGIN = 2, sum) %>% 
  as.data.frame()

EPICON_table_sum <- 
  EPICON_table_sum %>% mutate(
    sample_id = rownames(EPICON_table_sum)
  )

colnames(EPICON_table_sum)[1] <- "reads_num"

filterd_EPICON_sample <- EPICON_table_sum %>% 
  filter(
    reads_num > 5000
  ) %>% select(
    sample_id
  ) %>% pull()


EPICON_fungal_otutab1a1 <- EPICON_fungal_otutab1a %>% 
  select(all_of(filterd_EPICON_sample))


# EPICON_fungal_otutab1a1 %>% apply(MARGIN = 1, sum) %>% as.data.frame() %>% view()


EPICON_fungal_otutab1a1 %>% apply(MARGIN = 2, sum) %>% min()

# min abd is 5145
EPICON_fungal_otutab1a1_r <- rrarefy(EPICON_fungal_otutab1a1 %>% t(), 5145) %>% t()


apply(EPICON_fungal_otutab1a1_r, sum, MARGIN = 2)



rel_per <- function(val) {
  
  tmp_col <- val
  tmp_col_sum <- sum(tmp_col)
  
  re_per_col <- tmp_col / tmp_col_sum
  
  return(re_per_col)
  
}


EPICON_fungal_otutab1a1_per <- apply(EPICON_fungal_otutab1a1_r, rel_per, MARGIN = 2) %>% as.data.frame()
colnames(EPICON_fungal_otutab1a1_per)[1]

EPICON_fungal_otutab1a1_per1 <- 
  EPICON_fungal_otutab1a1_per %>% mutate(OTU_ID = rownames(EPICON_fungal_otutab1a1_per), .before = TP01L01) %>% 
  mutate(OTU_ID = str_split_i(OTU_ID, pattern = "_", 1)) %>% left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))


# EPICON_fungal_otutab1a1_per1 %>% view()

sample_id <- colnames(EPICON_fungal_otutab1a1_per)


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab1a1_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

EPICON_rDNAm_1a1 <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)

env_timepoint <- 
  data.frame(
    Timepiont = c("TP00", "TP01", "TP02", "TP03", "TP04", "TP05",
                  "TP06", "TP07", "TP08", "TP09", "TP10", "TP11",
                  "TP12", "TP13", "TP14", "TP15", "TP16", "TP17"),
    timepoint = c(0:17)
  )

env_timepoint

EPICON_rDNAm_1a1_env <- 
  EPICON_rDNAm_1a1 %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")

EPICON_rDNAm_1a1_env$Timepiont %>% table()
EPICON_rDNAm_1a1_env$Habitat %>% table()
EPICON_rDNAm_1a1_env$Treatment %>% table()




###### supp #######
# EPICON_rDNAm_env %>% view()
# EPICON_rDNAm_env_root_rhizo <- EPICON_rDNAm_env %>% filter(Habitat == "Rhizosphere" | Habitat == "Root")
# 
# EPICON_rDNAm_env_root_rhizo$sample_id
# 
# EPICON_otutab_sub <- EPICON_fungal_otutab1_per1 %>% select(OTU_ID, all_of(EPICON_rDNAm_env_root_rhizo$sample_id), rDNA_m)
# 
# sample_id_sub <- EPICON_rDNAm_env_root_rhizo$sample_id
# 
# EPICON_otutab_sub %>% view()
# 
# EPICON_otutab_sub_df <- 
#   apply(EPICON_otutab_sub %>% select(-OTU_ID, -rDNA_m), MARGIN = 2,
#         function(x) {x*EPICON_otutab_sub$rDNA_m}) %>% as.data.frame() %>% 
#   mutate(OTU_ID = EPICON_otutab_sub$OTU_ID, .before = TP01R1) %>% left_join(
#     EPICON_fungi_PL, by = c("OTU_ID" = "ID")
#   ) %>% drop_na()
# 
# EPICON_otutab_sub_df$TP01R1
# EPICON_otutab_sub_df
# 
# EPICON_otutab_sub_df_long <- EPICON_otutab_sub_df %>% 
#   pivot_longer(cols = sample_id_sub, values_to = "rDNA_each") %>% 
#   left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("name" = "aa"))
# 
# library(ggsci)
# ggplot(EPICON_otutab_sub_df_long, aes(Timepiont, rDNA_each, fill = Treatment
#                                       )) +
#   geom_col() +
#   scale_fill_d3(palette = "category20") +
#   facet_grid(Habitat ~ Subphylum, scales = "free")



# figS8, FRRN diff compartment -------

EPICON_com_1a1_kru <- 
  kruskal(EPICON_rDNAm_1a1_env$rDNAm, EPICON_rDNAm_1a1_env$Habitat, p.adj = "fdr")
EPICON_com_1a1_kru

EPICON_1a1_plot_subdata <- 
  EPICON_com_1a1_kru$groups %>% as.data.frame() %>% 
  mutate(Habitat = rownames(.), .before = groups) %>%
  select(Habitat, groups)
EPICON_1a1_plot_subdata

EPICON_com_1a1_maxrDNA <- EPICON_rDNAm_1a1_env %>% group_by(Habitat) %>% summarise(max_rDNA = max(rDNAm))

EPICON_1a1_plot_subdata1 <- EPICON_1a1_plot_subdata %>% left_join(EPICON_com_1a1_maxrDNA, by = "Habitat")

kruskal.test(rDNAm ~ Habitat, data = EPICON_rDNAm_1a1_env) %>% tidy()
# Kruskal-Wallis chi-squared = 428.13, df = 3, p-value < 2.2e-16

figS8_EPICON_com_1a1 <- 
  ggplot(EPICON_rDNAm_1a1_env, aes(Habitat, rDNAm)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Habitat), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.2, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Habitat),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = EPICON_plot_subdata1,
            aes(x = Habitat, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "blue",
            size = 7) +
  annotate(geom = "text", x = 2.5, y = 125, label = "chi-square = 428.13, df = 3, p = 1.79e-92", size = 5) +
  scale_y_continuous(limits = c(33, 130)) +
  scale_colour_manual(values = c("Leaf" = "darkgreen",
                                 "Rhizosphere" = "#ff00ff",
                                 "Root" = "navy",
                                 "Soil" = "brown"
  )) +
  scale_fill_manual(values = c("Leaf" = "darkgreen",
                               "Rhizosphere" = "#ff00ff",
                               "Root" = "navy",
                               "Soil" = "brown"
  )) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15),
    axis.title.y = element_text(face = "bold",
                                size = 15),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 1
  )
figS8_EPICON_com_1a1

tm <- now() %>% str_split_i(pattern = " ", 1)
figS8_pdf <- str_c("figS8_", "EPICON_com_", tm, ".pdf", sep = "")
figS8_jpg <- str_c("figS8_", "EPICON_com_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S8_pdf <- str_c(fig_path, figS8_pdf)
fig_fullpath_S8_jpg <- str_c(fig_path, figS8_jpg)

ggsave(fig_fullpath_S8_pdf, figS8_EPICON_com_1a1, width = 5.62, height = 5.42)
ggsave(fig_fullpath_S8_jpg, figS8_EPICON_com_1a1, width = 5.62, height = 5.42)


# envs
EPICON_rDNAm_1a1_env1 <- 
  EPICON_rDNAm_1a1_env %>% 
  mutate(grp = str_c(Habitat, Treatment, sep = "_"))
# EPICON_rDNAm_env1$Timepiont

EPICON_rDNAm_1a1_env1$Habitat %>% table()
EPICON_rDNAm_1a1_env1$Treatment %>% table()
EPICON_rDNAm_1a1_env1$grp %>% table()

# EPICON_rDNAm_plot <- 
#   ggplot(EPICON_rDNAm_env1) +
#   geom_point(
#              aes(Timepiont, rDNAm, colour = grp),
#              position = position_dodge(width = 0.7), size = 1) +
#   geom_boxplot(
#                aes(Timepiont, rDNAm, colour = grp),
#                position = position_dodge(width = 0.7), width = 0.4) +
#   geom_smooth(aes(timepoint + 1, rDNAm, colour = grp), se = F) +
#   geom_vline(xintercept = seq(1.5, 17.5, 1), linewidth = 0.05, linetype = "dashed") +
#   #annotate(geom = "text", x = 9.5, y = 96, label = "Habitat(H): df = 3, p = 2.52e-89") +
#   #annotate(geom = "text", x = 9.5, y = 94, label = "Treatment(T): df = 2, p = 1.11e-8") +
#   #annotate(geom = "text", x = 9.5, y = 92, label = "H * T: df = 6, p = 4.61e-42") +
#   scale_colour_manual(values = c("lightgreen", "green", "darkgreen",
#                                  "pink", "red", "darkred",
#                                  "lightblue", "blue", "darkblue",
#                                  "grey80", "grey40", "black")) +
#   guides(colour = guide_legend(
#     title = NULL,
#     ncol = 4)) +
#   labs(x = "Time", y = "rDNA copy number") + 
#   theme_classic() +
#   theme(legend.position = "inside",
#         legend.position.inside = c(0.5, 0.9),
#         axis.text = element_text(size = 16, colour = "black"),
#         axis.title = element_text(size = 20, colour = "black", face = "bold"),
#         legend.title = element_text(face = "bold", size = 18),
#         legend.key.size = unit(1, "cm"),
#         legend.text = element_text(size = 12),
#         legend.box.background = element_rect(colour = "black", linewidth = 1.5))
# EPICON_rDNAm_plot 

# ggsave("EPICON_rDNAm_plot_20241211_new.jpg", EPICON_rDNAm_plot)

# facet ------------
EPICON_rDNAm_1a1_env1 <- 
  EPICON_rDNAm_1a1_env1 %>% mutate(
    Treatment = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Post_flowering_drought" ~ "Post-flowering drought",
      Treatment == "Pre_flowering_drought" ~ "Pre-flowering drought"
    )
  )

EPICON_rDNAm_1a1_env1$Treatment <- factor(EPICON_rDNAm_1a1_env1$Treatment, levels = c("Control", "Pre-flowering drought", "Post-flowering drought"))

EPICON_rDNAm_1a1_env1$Timepiont

EPICON_rDNAm_1a1_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_rDNAm_1a1_env1)
EPICON_rDNAm_1a1_aov_sum <- summary(EPICON_rDNAm_1a1_aov)

EPICON_rDNAm_1a1_aov %>% tidy()

# EPICON_rDNAm_plot_facet_Treatment <- 
#   ggplot(EPICON_rDNAm_env1) +
#   geom_point(
#     aes(Timepiont, rDNAm, colour = Treatment),
#     position = position_dodge(width = 0.8), size = 2) +
#   geom_boxplot(
#     aes(Timepiont, rDNAm, colour = Treatment),
#     position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
#   geom_smooth(data = EPICON_rDNAm_env1 %>% filter(Timepiont != "TP00"),
#               aes(timepoint + 1, rDNAm, colour = Treatment), se = F) +
#   facet_wrap(~ Habitat, ncol = 1) +
#   scale_x_discrete(labels = seq(0, 17, 1)) +
#   scale_colour_manual(values = c("black",
#                                  "blue",
#                                  "red")) +
#   guides(colour = guide_legend(
#     title = "Treatment",
#     ncol = 1)) +
#   labs(x = "Week",
#        y = "Community-weighted rDNA copy number",
#        subtitle = "Compartment (C): df = 3, p = 1.64e-208; Treatment (T): df = 2, p = 4.08e-19; Week (W): df = 17, p = 7.04e-53\n
#        C * T: df = 6, p = 8.21e-67; C * W: df = 48, p = 8.94e-85; T * W: df = 25, p = 0.0181\n
#        C * T * W: df = 69, p = 2.76e-7") + 
#   theme_bw() +
#   theme(legend.position = "inside",
#         legend.position.inside = c(0.8, 0.91),
#         axis.text = element_text(size = 10, colour = "black"),
#         axis.title = element_text(size = 14, colour = "black", face = "bold"),
#         legend.title = element_text(face = "bold", size = 12),
#         strip.text = element_text(size = 12),
#         plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5))
# EPICON_rDNAm_plot_facet_Treatment
# 
# 
EPICON_rDNAm_1a1_env1 <-
  EPICON_rDNAm_1a1_env1 %>% mutate(
    show_type = if_else(timepoint == 0, 1, 0)
  )



# fig3c, EPICON ------------------------------
EPICON_rDNAm_1a1_aov %>% tidy()

EPICON_rDNAm_1a1_env1$timepoint

fig3c_EPICON_rDNAm_1a1_plot_Treatment_Compartment <- 
  ggplot(EPICON_rDNAm_1a1_env1) +
  geom_point(
    aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_1a1_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_rDNAm_1a1_env1 %>% filter(Timepiont != "TP00"),
              aes(timepoint + 1, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 9, y = 85, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "5.08e-258;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "4.62e-98;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "3.89e-16")) +
  annotate(geom = "text", x = 9, y = 81, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "1.20e-122;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "9.56e-101;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "4.88e-22")) +
  annotate(geom = "text", x = 9, y = 77, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "3.50e-12")) +
  scale_x_discrete(labels = seq(0, 17, 1)) +
  scale_colour_manual(values = c("darkgreen",
                                 "#ff00ff",
                                 "navy",
                                 "brown"
  )) +
  scale_y_continuous(limits = c(35, 87)) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  #scale_alpha_discrete(values = c(0, 1)) +
  guides(
    colour = guide_legend(
      title = "Compartment",
      ncol = 1,
      order = 1,
      direction = "vertical",
      override.aes = list(size = 3.5)),
    linetype = guide_legend(
      title = "Treatment",
      ncol = 1,
      order = 2,
      direction = "vertical")
  ) +
  labs(x = "Week",
       y = "Community-weighted rDNA copy number") + 
  theme_bw() +
  theme(legend.position = "right",
        #legend.position.inside = c(0.5, 0.05),
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 12),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        strip.text = element_text(size = 12),
        plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5),
        aspect.ratio = 0.618)
fig3c_EPICON_rDNAm_1a1_plot_Treatment_Compartment

tm <- now() %>% str_split_i(pattern = " ", 1)
fig3c_pdf <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".pdf", sep = "")
fig3c_jpg <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_3c_pdf <- str_c(fig_path, fig3c_pdf)
fig_fullpath_3c_jpg <- str_c(fig_path, fig3c_jpg)

ggsave(fig_fullpath_3c_pdf, fig3c_EPICON_rDNAm_1a1_plot_Treatment_Compartment, width = 11.5, height = 5.8)
ggsave(fig_fullpath_3c_jpg, fig3c_EPICON_rDNAm_1a1_plot_Treatment_Compartment, width = 11.5, height = 5.8)



fig3c_1_EPICON_rDNAm_1a1_plot_Treatment_Compartment <- 
  ggplot() +
  #geom_point(
  #  aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_rDNAm_1a1_env1 %>% filter(Timepiont != "TP00"),
              aes(timepoint, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 9, y = 85, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "5.08e-258;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "4.62e-98;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "3.89e-16")) +
  annotate(geom = "text", x = 9, y = 81, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "1.20e-122;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "9.56e-101;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "4.88e-22")) +
  annotate(geom = "text", x = 9, y = 77, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "3.50e-12")) +
  scale_x_continuous(labels = seq(1, 17, 1),
                     breaks = seq(1, 17, 1)) +
  scale_colour_manual(values = c("darkgreen",
                                 "#ff00ff",
                                 "navy",
                                 "brown"
  )) +
  scale_y_continuous(limits = c(35, 87)) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  #scale_alpha_discrete(values = c(0, 1)) +
  guides(
    colour = guide_legend(
      title = "Compartment",
      ncol = 1,
      order = 1,
      direction = "vertical",
      override.aes = list(size = 3.5)),
    linetype = guide_legend(
      title = "Treatment",
      ncol = 1,
      order = 2,
      direction = "vertical")
  ) +
  labs(x = "Week",
       y = "Community-weighted rDNA copy number") + 
  theme_bw() +
  theme(legend.position = "right",
        #legend.position.inside = c(0.5, 0.05),
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 12),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 15, colour = "black", face = "bold"),
        strip.text = element_text(size = 12),
        plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5),
        aspect.ratio = 0.618)
fig3c_1_EPICON_rDNAm_1a1_plot_Treatment_Compartment

tm <- now() %>% str_split_i(pattern = " ", 1)
fig3c_pdf <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".pdf", sep = "")
fig3c_jpg <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_3c_pdf <- str_c(fig_path, fig3c_pdf)
fig_fullpath_3c_jpg <- str_c(fig_path, fig3c_jpg)

ggsave(fig_fullpath_3c_pdf, fig3c_1_EPICON_rDNAm_1a1_plot_Treatment_Compartment, width = 11.5, height = 5.8)
ggsave(fig_fullpath_3c_jpg, fig3c_1_EPICON_rDNAm_1a1_plot_Treatment_Compartment, width = 11.5, height = 5.8)



# OTU_rDNAm_new <- EPICON_fungal_otutab1_per1 %>% select(OTU_ID, rDNA_m)
# SAM_rDNAm_new <- EPICON_rDNAm_env1 %>% select(sample_id, rDNAm)
# 
# write_xlsx(OTU_rDNAm_new, "OTU_rDNAm_new.xlsx")
# write_xlsx(SAM_rDNAm_new, "SAM_rDNAm_new.xlsx")


# different trait group in EPICON
FGanno_EPICON_FFF1 <- 
  FGanno_EPICON_FFF %>% mutate(
    Guild2 = case_when(
      str_detect(Guild1, "Saprotroph") ~ "Saprotroph fungi",
      str_detect(Guild1, "Plant Pathogen") ~ "Plant pathogen fungi",
      str_detect(Guild1, "Arbuscular Mycorrhizal") ~ "Arbuscular mycorrhizal fungi",
      .default = "Others"
    )
  )


EPICON_fungal_otutab1a1_1 <- EPICON_fungal_otutab1a1 %>% mutate(
  OTU_ID = str_split_i(rownames(EPICON_fungal_otutab1a1), i = 1, pattern = "_"), .before = TP01L01
)

FGanno_EPICON_FFF1$Guild2 %>% table()

# EPICON PP otutab -------------
EPICON_fungal_otutab_1a1_PP <- EPICON_fungal_otutab1a1_1 %>% 
  left_join(FGanno_EPICON_FFF1 %>% select(ID, Guild2), by = c("OTU_ID" = "ID")) %>%
  filter(Guild2 == "Plant pathogen fungi") %>% select(-Guild2)

EPICON_fungal_otutab_1a1_PP1 <- EPICON_fungal_otutab_1a1_PP %>% 
  mutate(sum_abd = apply(EPICON_fungal_otutab_1a1_PP %>% select(-OTU_ID), sum, MARGIN = 1)) %>% 
  filter(sum_abd != 0) %>% select(-sum_abd)

EPICON_fungal_otutab_1a1_PP1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% hist()
EPICON_fungal_otutab_1a1_PP1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% min()


EPICON_fungal_otutab_1a1_PP_r <- 
  rrarefy(EPICON_fungal_otutab_1a1_PP1 %>% 
            select(-OTU_ID) %>% t() , 776) %>% t() %>% 
  as.data.frame()

EPICON_fungal_otutab_1a1_PP_r %>% apply(sum, MARGIN = 2)

EPICON_fungal_otutab_1a1_PP_per <- apply(EPICON_fungal_otutab_1a1_PP_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = EPICON_fungal_otutab_1a1_PP1$OTU_ID, .before = TP01L01)


EPICON_fungal_otutab_1a1_PP_per1 <- 
  EPICON_fungal_otutab_1a1_PP_per %>% 
  left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))

sample_id <- colnames(EPICON_fungal_otutab_1a1_PP_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab_1a1_PP_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

EPICON_fungal_otutab_1a1_PP_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
EPICON_fungal_otutab_1a1_PP_rDNAm


EPICON_1a1_PP_rDNAm <- 
  EPICON_fungal_otutab_1a1_PP_rDNAm %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")


EPICON_1a1_PP_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_1a1_PP_rDNAm)
EPICON_1a1_PP_rDNAm_aov_sum <- summary(EPICON_1a1_PP_rDNAm_aov)

EPICON_1a1_PP_rDNAm_aov %>% tidy()

# figS10, plant pathogen --------------------------------
figS9a_EPICON_1a1_PP <- 
  ggplot(EPICON_1a1_PP_rDNAm) +
  #geom_point(
  #  aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_PP_rDNAm %>% filter(Timepiont != "TP00"),
              aes(timepoint, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 9, y = 90, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "1.46e-209;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "3.21e-50;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "3.46e-17")) +
  annotate(geom = "text", x = 9, y = 87, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "2.39e-52;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "3.41e-58;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "5.00e-10")) +
  annotate(geom = "text", x = 9, y = 84, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "6.72e-8")) +
  #annotate(geom = "text", x = 9, y = 100,
  #         label = "Compartment (C): df = 3, p = 1.46e-164; Treatment (T): df = 2, p = 3.19e-19; Week (W): df = 17, p = 1.01e-77\n
  #         C * T: df = 6, p = 8.23e-90; C * W: df = 48, p = 9.80e-95; T * W: df = 25, p = 2.72e-6\n
  #         C * T * W: df = 69, p = 7.71e-7", hjust = 0.5, vjust = 0.5,
  #         fontface = "bold",
  #         size = 4) +
  # facet_wrap(~ Treatment, ncol = 1) +
  scale_x_continuous(labels = seq(1, 17, 1), breaks = seq(1, 17, 1), limits = c(1, 17)) +
  #scale_y_continuous(limits = c(50, 85)) +
  scale_colour_manual(values = c("darkgreen",
                                 "#ff00ff",
                                 "navy",
                                 "brown"
  )) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  #scale_alpha_discrete(values = c(0, 1)) +
  guides(
    colour = guide_legend(
      title = "Compartment",
      ncol = 1,
      order = 1,
      direction = "vertical",
      override.aes = list(size = 3.5)),
    linetype = guide_legend(
      title = "Treatment",
      ncol = 1,
      order = 2,
      direction = "vertical")
  ) +
  labs(x = "Week",
       y = "Community-weighted rDNA copy number of plant pathogenic fungi") + 
  theme_bw() +
  theme(legend.position = "right",
        #legend.position.inside = c(0.5, 0.05),
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 12),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        strip.text = element_text(size = 12),
        plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5),
        aspect.ratio = 0.618)
figS9a_EPICON_1a1_PP

tm <- now() %>% str_split_i(pattern = " ", 1)
figS9a_pdf <- str_c("figS9a_", "EPICON_PP_", tm, ".pdf", sep = "")
figS9a_jpg <- str_c("figS9a_", "EPICON_PP_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S9a_pdf <- str_c(fig_path, figS9a_pdf)
fig_fullpath_S9a_jpg <- str_c(fig_path, figS9a_jpg)

ggsave(fig_fullpath_S9a_pdf, figS9a_EPICON_1a1_PP, width = 11.5, height = 5.8)
ggsave(fig_fullpath_S9a_jpg, figS9a_EPICON_1a1_PP, width = 11.5, height = 5.8)


# EPICON Sap otutab ------
EPICON_fungal_otutab_1a1_Sap <- EPICON_fungal_otutab1a1_1 %>% 
  left_join(FGanno_EPICON_FFF1 %>% select(ID, Guild2), by = c("OTU_ID" = "ID")) %>%
  filter(Guild2 == "Saprotroph fungi") %>% select(-Guild2)

EPICON_fungal_otutab_1a1_Sap1 <- EPICON_fungal_otutab_1a1_Sap %>% 
  mutate(sum_abd = apply(EPICON_fungal_otutab_1a1_Sap %>% select(-OTU_ID), sum, MARGIN = 1)) %>% 
  filter(sum_abd != 0) %>% select(-sum_abd)

# EPICON_fungal_otutab_1a1_Sap1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% min()


EPICON_fungal_otutab_1a1_Sap_r <- 
  rrarefy(EPICON_fungal_otutab_1a1_Sap1 %>% select(-OTU_ID) %>% t() , 338) %>% t() %>% as.data.frame()

EPICON_fungal_otutab_1a1_Sap_r %>% apply(sum, MARGIN = 2)


EPICON_fungal_otutab_1a1_Sap_per <- apply(EPICON_fungal_otutab_1a1_Sap_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = EPICON_fungal_otutab_1a1_Sap1$OTU_ID, .before = TP01L01)

EPICON_fungal_otutab_1a1_Sap_per1 <- 
  EPICON_fungal_otutab_1a1_Sap_per %>% 
  left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))

sample_id <- colnames(EPICON_fungal_otutab_1a1_Sap_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab_1a1_Sap_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

EPICON_fungal_otutab_1a1_Sap_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
EPICON_fungal_otutab_1a1_Sap_rDNAm

EPICON_1a1_Sap_rDNAm <- 
  EPICON_fungal_otutab_1a1_Sap_rDNAm %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")


EPICON_1a1_Sap_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_1a1_Sap_rDNAm)
EPICON_1a1_Sap_rDNAm_aov_sum <- summary(EPICON_1a1_Sap_rDNAm_aov)

EPICON_1a1_Sap_rDNAm_aov %>% tidy()

# figS9, saprotroph -------------------------------
figS9b_EPICON_1a1_Sap <- 
  ggplot(EPICON_1a1_Sap_rDNAm) +
  #geom_point(
  #  aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_Sap_rDNAm %>% filter(Timepiont != "TP00"),
              aes(timepoint, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 9, y = 83, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "1.36e-141;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "2.97e-32;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "9.81e-14")) +
  annotate(geom = "text", x = 9, y = 79, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "4.50e-39;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "8.97e-44;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "7.90e-6")) +
  annotate(geom = "text", x = 9, y = 75, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "5.85e-8")) +
  #annotate(geom = "text", x = 9, y = 100,
  #         label = "Compartment (C): df = 3, p = 1.46e-164; Treatment (T): df = 2, p = 3.19e-19; Week (W): df = 17, p = 1.01e-77\n
  #         C * T: df = 6, p = 8.23e-90; C * W: df = 48, p = 9.80e-95; T * W: df = 25, p = 2.72e-6\n
  #         C * T * W: df = 69, p = 7.71e-7", hjust = 0.5, vjust = 0.5,
  #         fontface = "bold",
  #         size = 4) +
  # facet_wrap(~ Treatment, ncol = 1) +
  scale_x_continuous(labels = seq(1, 17, 1), breaks = seq(1, 17, 1), limits = c(1, 17)) +
  # scale_y_continuous(limits = c(30, 82)) +
  scale_colour_manual(values = c("darkgreen",
                                 "#ff00ff",
                                 "navy",
                                 "brown"
  )) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  #scale_alpha_discrete(values = c(0, 1)) +
  guides(
    colour = guide_legend(
      title = "Compartment",
      ncol = 1,
      order = 1,
      direction = "vertical",
      override.aes = list(size = 3.5)),
    linetype = guide_legend(
      title = "Treatment",
      ncol = 1,
      order = 2,
      direction = "vertical")
  ) +
  labs(x = "Week",
       y = "Community-weighted rDNA copy number of saprotrophic fungi") + 
  theme_bw() +
  theme(legend.position = "right",
        #legend.position.inside = c(0.5, 0.05),
        legend.box = "vertical",
        legend.title = element_text(face = "bold", size = 12),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 12, colour = "black", face = "bold"),
        strip.text = element_text(size = 12),
        plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5),
        aspect.ratio = 0.618)
figS9b_EPICON_1a1_Sap


tm <- now() %>% str_split_i(pattern = " ", 1)
figS9b_pdf <- str_c("figS9b_", "EPICON_Sap_", tm, ".pdf", sep = "")
figS9b_jpg <- str_c("figS9b_", "EPICON_Sap_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S9b_pdf <- str_c(fig_path, figS9b_pdf)
fig_fullpath_S9b_jpg <- str_c(fig_path, figS9b_jpg)

ggsave(fig_fullpath_S9b_pdf, figS9b_EPICON_1a1_Sap, width = 11.5, height = 5.8)
ggsave(fig_fullpath_S9b_jpg, figS9b_EPICON_1a1_Sap, width = 11.5, height = 5.8)

# figS9
figS9 <- figS9a_EPICON_1a1_PP + figS9b_EPICON_1a1_Sap + plot_layout(guides = "collect")
figS9

tm <- now() %>% str_split_i(pattern = " ", 1)
figS9_pdf <- str_c("figS9_", "EPICON_PP_Sap_", tm, ".pdf", sep = "")
figS9_jpg <- str_c("figS9_", "EPICON_PP_Sap_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S9_pdf <- str_c(fig_path, figS9_pdf)
fig_fullpath_S9_jpg <- str_c(fig_path, figS9_jpg)

ggsave(fig_fullpath_S9_pdf, figS9, width = 20, height = 5.8)
ggsave(fig_fullpath_S9_jpg, figS9, width = 20, height = 5.8)

##### EPICON lm #####
# new function #
# EPICON_leaf_lm_data <- EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Leaf")
# colnames(EPICON_leaf_lm_data)[5] <- "group"
# 
# 
# EPICON_1a1_leaf_lm <-
#   subplot_data_corr(DF = EPICON_leaf_lm_data, Y_val = "rDNAm", X_val = "timepoint",
#                     GROUP_list = c("Control", "Pre-flowering drought", "Post-flowering drought"),
#                     METHOD = "spearman", adj_METHOD = "fdr")
# EPICON_1a1_leaf_lm
# should be checked

# EPICON lm result ----------------
# old version #
lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Treatment == subgrp)
    
  } else {
    
    df_tmp <- df
    
  }
  
  df_tmp <- df_tmp %>% drop_na()
  
  sub_lmR <- lm(df_tmp[[yval]] ~ df_tmp[[xval]])
  lmR_sumy <- summary(sub_lmR)
  
  if(lmR_sumy$coefficients[2] > 0) {
    lmR_r <- round(sqrt(lmR_sumy$r.squared), 3) %>% signif(., 3)
  } else {
    lmR_r <- round(sqrt(lmR_sumy$r.squared), 3) %>% signif(., 3) * -1
  }
  
  lmR_P_tmpf <- lmR_sumy$fstatistic
  lmR_P_tmp <- pf(lmR_P_tmpf[1], lmR_P_tmpf[2], lmR_P_tmpf[3], lower.tail = F)
  
  lmR_p <- round(lmR_P_tmp, 3) %>% signif(., 3)
  
  lmR_slope <- lmR_sumy$coefficients[2]
  
  if(is.na(lmR_p)) {
    
    Rp_sig <- "NA"
    
  } else if(lmR_p < 0.001) {
    
    Rp_sig <- "***"
    
  } else if(lmR_p <= 0.01) {
    
    Rp_sig <- "**"
    
  } else if(lmR_p <= 0.05){
    
    Rp_sig <- "*"
    
  } else if(lmR_p <= 1){
    
    Rp_sig <- "NS"
    
  }
  
  if(trans != "no_trans") {
    
    df_rlt <- data.frame(
      Treatment = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_P_tmp,
      slope = lmR_slope
    )
    
  } else {
    
    df_rlt <- data.frame(
      Treatment = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_P_tmp,
      slope = lmR_slope
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


EPICON_leaf_1a1_lm <-
  map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
          .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Leaf"),
                             yval = "rDNAm", xval = "timepoint",
                             subgrp = ., trans = "no_trans"))
EPICON_leaf_1a1_lm

EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Leaf") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_m = mean(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_soil_1a1_lm <-
  map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
          .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Soil"),
                             yval = "rDNAm", xval = "timepoint",
                             subgrp = ., trans = "no_trans"))
EPICON_soil_1a1_lm

EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Soil") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_m = mean(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Root") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_m = mean(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_rDNAm_1a1_env1 %>% filter(Habitat == "Rhizosphere") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_m = mean(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()

# done ...

###### XC-Li, leat and root fungi ######

# meta data
XC_METADATA <- read.csv("./2.database/fun.env_510.csv", header = T)

XC_METADATA$Sample.No. %>% length()

# XC ROOT
XC_ROOT <- read_csv("./2.database/bfp.fungi.matrix.csv")

XC_ROOT_TAX <- XC_ROOT %>% select(1:11)
XC_ROOT_TABLE <- XC_ROOT %>% select(-c(2:11)) 


# filtered TABLE by TAX
# TAX
XC_ROOT_TAX1 <- 
  XC_ROOT_TAX %>% filter(
    Phylum != "."
  ) %>% select(
    -blast.unite.O1
  ) %>% select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species, EMF) %>% 
  mutate(Species = str_sub(Species, start = 4)) %>% 
  mutate(Species = str_replace(Species, pattern = "\\.", " "))
XC_ROOT_TAX1


# filtered
XC_ROOT_TABLE1 <- XC_ROOT_TABLE %>% filter(
  OTU %in% XC_ROOT_TAX1$OTU
)


XC_ROOT_abd <- 
  XC_ROOT_TABLE1 %>% select(-1) %>% apply(MARGIN = 2, sum) %>% 
  as.data.frame()

XC_ROOT_abd1 <- 
  XC_ROOT_abd %>% mutate(
    sample_id = rownames(XC_ROOT_abd)
  )

colnames(XC_ROOT_abd1)[1] <- "reads_num"

filterd_XC_ROOT_sample <- XC_ROOT_abd1 %>% 
  filter(
    reads_num > 5000
  ) %>% select(
    sample_id
  ) %>% pull()


XC_ROOT_TABLE2 <- XC_ROOT_TABLE1 %>% 
  select(1, all_of(filterd_XC_ROOT_sample))

XC_ROOT_TABLE2


###### new corr function #####
subplot_data_corr <- function(DF, X_val, Y_val, GROUP_list = c("all"), METHOD, adj_METHOD) {
  
  
  # main function
  corr_rlts <- function(df, yval, xval, subgrp = c("all"), corr_method) {
    
    
    # data
    if(subgrp != "all") {
      
      df_tmp <- df %>% filter(group == subgrp)
      
    } else {
      
      df_tmp <- df
      
    }
    
    
    # cor test
    lmR_spearman <- cor.test(df_tmp[[yval]], df_tmp[[xval]], method = corr_method)
    lmR_spe_rho <- lmR_spearman$estimate %>% round(3)
    lmR_spe_pval <- lmR_spearman$p.value  
    
    
    # results
    df_rlt <- data.frame(
      group = subgrp,
      anno_x1 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.8 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(rho) == ", lmR_spe_rho),
      #psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_spe_pval
    )
    
    return(df_rlt)
    
  }
  
  
  
  # corr temp results
  corr_tmp_results <- 
    map_dfr(GROUP_list, ~ corr_rlts(df = DF,
                                    yval = Y_val,
                                    xval = X_val,
                                    corr_method = METHOD,
                                    subgrp = .))
  
  # final output
  corr_results <- corr_tmp_results %>% 
    mutate(
      pval_sig = case_when(
        pval <= 0.001 ~ "***",
        pval > 0.001 & pval <= 0.01 ~ "**",
        pval > 0.01 & pval <= 0.05 ~ "*",
        .default = "NS"
      ),
      .after = pval
    ) %>% 
    mutate(
      pval_sig_tmp = str_c("italic(p) == ", pval %>% formatC(format = "e", digits = 3))
    ) %>% 
    mutate(
      rlt_sig = str_c(rsig, " ~~ ", pval_sig_tmp)
    ) %>% 
    mutate(
      pval_adj = p.adjust(pval, method = adj_METHOD)
    ) %>% 
    mutate(
      pval_adj_sig = case_when(
        pval_adj <= 0.001 ~ "***",
        pval_adj > 0.001 & pval <= 0.01 ~ "**",
        pval_adj > 0.01 & pval <= 0.05 ~ "*",
        .default = "NS"
      ),
      .after = pval_adj
    ) %>% 
    mutate(
      pval_adj_sig_tmp = str_c("italic(p) == ", pval_adj %>% formatC(format = "e", digits = 3))
    ) %>% 
    mutate(
      rlt_adj_sig = str_c(rsig, " ~~ ", pval_adj_sig_tmp)
    ) %>% select(-pval_sig_tmp, -pval_adj_sig_tmp)
  
  
  return(corr_results)
  
  
}


# root -----------
XC_METADATA
XC_ROOT_TABLE2
XC_ROOT_TAX1


# FRRN tab
FRRN_phy_tab
FRRN_cla_tab
FRRN_ord_tab
FRRN_fam_tab
FRRN_gen_tab
FRRN_spc_tab

# spc
FRRN_XC_ROOT_spc <- XC_ROOT_TAX1 %>% left_join(FRRN_spc_tab, by = c("Species" = "spc")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_spc

FRRN_XC_ROOT_spc_na <- XC_ROOT_TAX1 %>% left_join(FRRN_spc_tab, by = c("Species" = "spc")) %>% 
  filter(is.na(rDNA_m))

# gen
FRRN_XC_ROOT_gen <- FRRN_XC_ROOT_spc_na %>% select(-rDNA_m) %>% left_join(FRRN_gen_tab, by = c("Genus" = "gen")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_gen

FRRN_XC_ROOT_gen_na <- FRRN_XC_ROOT_spc_na %>% select(-rDNA_m) %>% left_join(FRRN_gen_tab, by = c("Genus" = "gen")) %>% 
  filter(is.na(rDNA_m))

# fam
FRRN_XC_ROOT_fam <- FRRN_XC_ROOT_gen_na %>% select(-rDNA_m) %>% left_join(FRRN_fam_tab, by = c("Family" = "fam")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_fam

FRRN_XC_ROOT_fam_na <- FRRN_XC_ROOT_gen_na %>% select(-rDNA_m) %>% left_join(FRRN_fam_tab, by = c("Family" = "fam")) %>% 
  filter(is.na(rDNA_m))


# ord
FRRN_XC_ROOT_ord <- FRRN_XC_ROOT_fam_na %>% select(-rDNA_m) %>% left_join(FRRN_ord_tab, by = c("Order" = "ord")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_ord

FRRN_XC_ROOT_ord_na <- FRRN_XC_ROOT_fam_na %>% select(-rDNA_m) %>% left_join(FRRN_ord_tab, by = c("Order" = "ord")) %>% 
  filter(is.na(rDNA_m))


# cla
FRRN_XC_ROOT_cla <- FRRN_XC_ROOT_ord_na %>% select(-rDNA_m) %>% left_join(FRRN_cla_tab, by = c("Class" = "cla")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_cla

FRRN_XC_ROOT_cla_na <- FRRN_XC_ROOT_ord_na %>% select(-rDNA_m) %>% left_join(FRRN_cla_tab, by = c("Class" = "cla")) %>% 
  filter(is.na(rDNA_m))


# phy
FRRN_XC_ROOT_phy <- FRRN_XC_ROOT_cla_na %>% select(-rDNA_m) %>% left_join(FRRN_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(!is.na(rDNA_m))
FRRN_XC_ROOT_phy

FRRN_XC_ROOT_phy_na <- FRRN_XC_ROOT_cla_na %>% select(-rDNA_m) %>% left_join(FRRN_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(is.na(rDNA_m))
# FRRN_XCroot_phy_na %>% dim()
# FRRN_XCroot_phy_na %>% view()



# FRRN XCroot anno --------------------------
FRRN_XC_ROOT_anno <- rbind(FRRN_XC_ROOT_phy, FRRN_XC_ROOT_cla, FRRN_XC_ROOT_ord,
                           FRRN_XC_ROOT_fam, FRRN_XC_ROOT_gen, FRRN_XC_ROOT_spc)
# dim(FRRN_XC_ROOT_anno)
# done

# colnames(XC_ROOT_TABLE)

# rrarefy
XC_ROOT_TABLE2 %>% select(-1) %>% apply(MARGIN = 2, sum) %>% min()

# min abd is 5201
XC_ROOT_TABLE2_r <- rrarefy(XC_ROOT_TABLE2 %>% select(-OTU) %>% t(), 5201) %>% t()

# apply(XC_ROOT_TABLE2_r, sum, MARGIN = 2)

XC_ROOT_TABLE_per <- apply(XC_ROOT_TABLE2_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU = XC_ROOT_TABLE2$OTU, .before = bet001)

# apply(XC_ROOT_TABLE_per %>% select(-OTU), sum, MARGIN = 2)

# view(XC_ROOT_TABLE_per)

XC_ROOT_TABLE_per1 <- 
  XC_ROOT_TABLE_per %>% 
  left_join(FRRN_XC_ROOT_anno %>% select(OTU, rDNA_m), by = "OTU")

sample_id <- colnames(XC_ROOT_TABLE_per1 %>% select(-OTU, -rDNA_m))
sample_id

rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- XC_ROOT_TABLE_per1 %>% select(OTU, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

XC_ROOT_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
XC_ROOT_rDNAm

XC_ROOT_rDNAm1 <- XC_ROOT_rDNAm %>% left_join(XC_METADATA, by = c("sample_id" = "Sample.No."))

#view(XCroot_rDNAm1)
#XCroot_rDNAm1 %>% select(group)


# ROOT MAT
XC_ROOT_group_list <- XC_ROOT_rDNAm1$group %>% unique()

XC_ROOT_MAT_cor_subdata <- 
  subplot_data_corr(DF = XC_ROOT_rDNAm1, Y_val = "rDNAm", X_val = "MAT",
                    GROUP_list = XC_ROOT_group_list, METHOD = "spearman", adj_METHOD = "fdr")
XC_ROOT_MAT_cor_subdata


# help("geom_hex")
XC_ROOT_rDNAm_MAT <- 
  ggplot(XC_ROOT_rDNAm1, aes(MAT, rDNAm)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 0.4) +
  # geom_bin2d() +
  # geom_hex() +
  # geom_smooth(colour = "red", method = "lm") +
  geom_smooth(method = "loess", span = 1, colour = "red") +
  geom_text(data = XC_ROOT_MAT_cor_subdata,
            aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
            parse = T,
            colour = "blue",
            size = 4.5,
            vjust = 3) +
  # scale_fill_gradient2(low = "grey90", high = "brown") +
  # facet_wrap(~ group, scales = "free") +
  labs(x = "Mean annual temperature", y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
XC_ROOT_rDNAm_MAT

tm <- now() %>% str_split_i(pattern = " ", 1)
figS5a_pdf <- str_c("figS5a_", "XC_ROOT_MAT_", tm, ".pdf", sep = "")
figS5a_jpg <- str_c("figS5a_", "XC_ROOT_MAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S5a_pdf <- str_c(fig_path, figS5a_pdf)
fig_fullpath_S5a_jpg <- str_c(fig_path, figS5a_jpg)

ggsave(fig_fullpath_S5a_pdf, XC_ROOT_rDNAm_MAT, width = 5.05, height = 4.97)
ggsave(fig_fullpath_S5a_jpg, XC_ROOT_rDNAm_MAT, width = 5.05, height = 4.97)


# "#3B4992FF" "#EE0000FF" "#008B45FF"


# ROOT MAP
XC_ROOT_MAP_cor_subdata <- 
  subplot_data_corr(DF = XC_ROOT_rDNAm1, Y_val = "rDNAm", X_val = "MAP",
                    GROUP_list = XC_ROOT_group_list, METHOD = "spearman", adj_METHOD = "fdr")
XC_ROOT_MAP_cor_subdata

XC_ROOT_rDNAm_MAP <- 
  ggplot(XC_ROOT_rDNAm1, aes(MAP, rDNAm)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 50) +
  # geom_smooth(colour = "red", method = "lm") +
  # geom_hex() +
  geom_smooth(method = "loess", span = 1, colour = "red") +
  geom_text(data = XC_ROOT_MAP_cor_subdata,
            aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
            parse = T,
            colour = "blue",
            size = 4.5,
            vjust = 3) +
  scale_x_continuous(
    breaks = seq(500, 2500, 1000),
    labels = seq(500, 2500, 1000)
  ) +
  # scale_fill_gradient2(low = "grey90", high = "brown") +
  # facet_wrap(~ group, scales = "free") +
  labs(x = "Mean annual precipitation", y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
XC_ROOT_rDNAm_MAP

tm <- now() %>% str_split_i(pattern = " ", 1)
figS5b_pdf <- str_c("figS5b_", "XC_ROOT_MAP_", tm, ".pdf", sep = "")
figS5b_jpg <- str_c("figS5b_", "XC_ROOT_MAP_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S5b_pdf <- str_c(fig_path, figS5b_pdf)
fig_fullpath_S5b_jpg <- str_c(fig_path, figS5b_jpg)

ggsave(fig_fullpath_S5b_pdf, XC_ROOT_rDNAm_MAP, width = 5.05, height = 4.97)
ggsave(fig_fullpath_S5b_jpg, XC_ROOT_rDNAm_MAP, width = 5.05, height = 4.97)

# XC ROOT MAT and MAP 
XC_ROOT_MATMAP <- XC_ROOT_rDNAm_MAT + XC_ROOT_rDNAm_MAP
ggsave("figS5_XC_ROOT_MATMAP.pdf", XC_ROOT_MATMAP, width = 9.47, height = 4.7)
ggsave("figS5_XC_ROOT_MATMAP.jpg", XC_ROOT_MATMAP, width = 9.47, height = 4.7)



# XC ROOT EMF ----------

XC_ROOT_EMFind <- XC_ROOT_TAX %>% filter(EMF == "emf") %>% select(OTU) %>% pull()
XC_ROOT_TABLE2_EMF <- XC_ROOT_TABLE2 %>% filter(OTU %in% XC_ROOT_EMFind)

EMF_sample_sub <- apply(XC_ROOT_TABLE2_EMF %>% select(-1), sum, MARGIN = 2) %>%
  as.data.frame()

colnames(EMF_sample_sub) <- "Freq"
# EMF_sample_sub$Freq %>% hist()


EMF_sample_sub1 <- EMF_sample_sub %>% mutate(sample_id = rownames(EMF_sample_sub), .before = Freq) %>% 
  filter(Freq > 5000)


XC_ROOT_TABLE2_EMF1 <- XC_ROOT_TABLE2_EMF %>% select(OTU, all_of(EMF_sample_sub1$sample_id))

# apply(XC_ROOT_TABLE2_EMF1 %>% select(-OTU), sum, MARGIN = 2) %>% min()

XC_ROOT_TABLE2_EMF1_r <- rrarefy(XC_ROOT_TABLE2_EMF1 %>% select(-OTU) %>% t(), 5021) %>% t() %>% as.data.frame()

XC_ROOT_TABLE2_EMF1_r <- XC_ROOT_TABLE2_EMF1_r %>% mutate(OTU = XC_ROOT_TABLE2_EMF1$OTU, .before = bet001)


XC_ROOT_TABLE2_EMF1_per <- apply(XC_ROOT_TABLE2_EMF1_r %>% select(-OTU), rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU = XC_ROOT_TABLE2_EMF1_r$OTU, .before = bet001)

XC_ROOT_TABLE2_EMF1_per1 <- 
  XC_ROOT_TABLE2_EMF1_per %>% 
  left_join(FRRN_XC_ROOT_anno %>% select(OTU, rDNA_m), by = "OTU")

sample_id <- colnames(XC_ROOT_TABLE2_EMF1_per1 %>% select(-OTU, -rDNA_m))
sample_id

rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- XC_ROOT_TABLE2_EMF1_per1 %>% select(OTU, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

XC_ROOT_EMF_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
XC_ROOT_EMF_rDNAm

XC_ROOT_EMF_rDNAm1 <- XC_ROOT_EMF_rDNAm %>% left_join(XC_METADATA, by = c("sample_id" = "Sample.No."))

# XC_ROOT_EMF_rDNAm1 %>% view()

# XC ROOT EMF LAT
XC_ROOT_EMF_rDNAm2 <- XC_ROOT_EMF_rDNAm1 %>% select(-group)
colnames(XC_ROOT_EMF_rDNAm2)[10] <- "group"

XC_ROOT_EMF_group_list <- XC_ROOT_EMF_rDNAm2$group %>% unique()

XC_ROOT_EMF_LAT_cor_subdata <- 
  subplot_data_corr(DF = XC_ROOT_EMF_rDNAm2, Y_val = "rDNAm", X_val = "Latitude",
                    GROUP_list = XC_ROOT_EMF_group_list, METHOD = "spearman", adj_METHOD = "fdr")
XC_ROOT_EMF_LAT_cor_subdata


# XC_ROOT_EMF_rDNAm_LAT <- 
#   ggplot(XC_ROOT_EMF_rDNAm2, aes(Latitude, rDNAm, colour = group)) +
#   geom_jitter(size = 2.8, alpha = 0.4, width = 0.4) +
#   # geom_hex() +
#   geom_text(data = XC_ROOT_EMF_LAT_cor_subdata,
#             aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
#             parse = T,
#             colour = "blue",
#             size = 4.5,
#             vjust = 2) + 
#   geom_smooth(colour = "red") +
#   scale_colour_nejm() +
#   facet_wrap(~ group, scales = "free") +
#   guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
#                                nrow = 1,
#                                direction = "horizontal")) +
#   labs(y = "Community-weighted rDNA copy number\nof ectomycorrhizal fungi") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12, face = "bold.italic"),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "none",
#         aspect.ratio = 1)
# XC_ROOT_EMF_rDNAm_LAT
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS16_pdf <- str_c("figS16_", "XCroot_EMF_LAT_", tm, ".pdf", sep = "")
# figS16_jpg <- str_c("figS16_", "XCroot_EMF_LAT_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S16_pdf <- str_c(fig_path, figS16_pdf)
# fig_fullpath_S16_jpg <- str_c(fig_path, figS16_jpg)
# 
# ggsave(fig_fullpath_S16_pdf, XCroot_EMF_rDNAm_LAT, width = 11.6, height = 4.45)
# ggsave(fig_fullpath_S16_jpg, XCroot_EMF_rDNAm_LAT, width = 11.6, height = 4.45)

# XC ROOT full ECM vs LAT
XC_ROOT_EMF_rDNAm3 <- XC_ROOT_EMF_rDNAm1 %>% mutate(
  group = "ROOT_ECM"
)

# view(XC_ROOT_EMF_rDNAm2)

XC_ROOT_EMF_group_list <- XC_ROOT_EMF_rDNAm3$group %>% unique()

XC_ROOT_EMF_LAT_cor_subdata <- 
  subplot_data_corr(DF = XC_ROOT_EMF_rDNAm3, Y_val = "rDNAm", X_val = "Latitude",
                    GROUP_list = XC_ROOT_EMF_group_list, METHOD = "spearman", adj_METHOD = "fdr")
XC_ROOT_EMF_LAT_cor_subdata


XC_ROOT_EMF_rDNAm_LAT_full <- 
  ggplot(XC_ROOT_EMF_rDNAm3, aes(Latitude, rDNAm, colour = Family)) +
  geom_jitter(size = 2.5, alpha = 0.5, width = 0.4) +
  # geom_hex() +
  geom_text(data = XC_ROOT_EMF_LAT_cor_subdata,
            aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
            parse = T,
            colour = "blue",
            size = 4.5,
            vjust = 4.5) + 
  geom_smooth(colour = "red", method = "lm") +
  scale_colour_nejm() +
  # facet_wrap(~ group, scales = "free") +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal")) +
  labs(y = "Community-weighted rDNA copy number\nof ectomycorrhizal fungi") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.94))
XC_ROOT_EMF_rDNAm_LAT_full

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7a_pdf <- str_c("figS7a_", "XCroot_EMF_LAT_", tm, ".pdf", sep = "")
figS7a_jpg <- str_c("figS7a_", "XCroot_EMF_LAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7a_pdf <- str_c(fig_path, figS7a_pdf)
fig_fullpath_S7a_jpg <- str_c(fig_path, figS7a_jpg)

ggsave(fig_fullpath_S7a_pdf, XC_ROOT_EMF_rDNAm_LAT_full, width = 5.11, height = 4.67)
ggsave(fig_fullpath_S7a_jpg, XC_ROOT_EMF_rDNAm_LAT_full, width = 5.11, height = 4.67)




# Zheng data another version ----------
Z_soil_meta_cfb <- read.csv("./2.database/Zheng_CFB.env.2022.03.28.csv", header = T)
Z_soil_taxa_cfb <- read.csv("./2.database/Zheng_CFB.fung.ID.2022.03.28.csv", header = T)
Z_soil_otutab_cfb <- read.csv("./2.database/Zheng_CFB.fung.2022.03.28.csv", header = T)


# Z another anno
FRRN_phy_tab
FRRN_cla_tab
FRRN_ord_tab
FRRN_fam_tab
FRRN_gen_tab
FRRN_spc_tab


# taxa tab
Z_soil_taxa_cfb1 <- Z_soil_taxa_cfb %>% select(OTU, Phylum, Class, Order, Family, Genus, Species, Lifestyle2) %>% 
  mutate(Species = str_replace(Species, pattern = "_", replacement = " "))
#view(Z_soil_taxa_cfb1)


# spc
FRRN_Zsoil_cfb_spc <- Z_soil_taxa_cfb1 %>% left_join(FRRN_spc_tab, by = c("Species" = "spc")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_spc

FRRN_Zsoil_cfb_spc_na <- Z_soil_taxa_cfb1 %>% left_join(FRRN_spc_tab, by = c("Species" = "spc")) %>% 
  filter(is.na(rDNA_m))


# gen
FRRN_Zsoil_cfb_gen <- FRRN_Zsoil_cfb_spc_na %>% select(-rDNA_m) %>% left_join(FRRN_gen_tab, by = c("Genus" = "gen")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_gen

FRRN_Zsoil_cfb_gen_na <- FRRN_Zsoil_cfb_spc_na %>% select(-rDNA_m) %>% left_join(FRRN_gen_tab, by = c("Genus" = "gen")) %>% 
  filter(is.na(rDNA_m))

# fam
FRRN_Zsoil_cfb_fam <- FRRN_Zsoil_cfb_gen_na %>% select(-rDNA_m) %>% left_join(FRRN_fam_tab, by = c("Family" = "fam")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_fam

FRRN_Zsoil_cfb_fam_na <- FRRN_Zsoil_cfb_gen_na %>% select(-rDNA_m) %>% left_join(FRRN_fam_tab, by = c("Family" = "fam")) %>% 
  filter(is.na(rDNA_m))


# ord
FRRN_Zsoil_cfb_ord <- FRRN_Zsoil_cfb_fam_na %>% select(-rDNA_m) %>% left_join(FRRN_ord_tab, by = c("Order" = "ord")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_ord

FRRN_Zsoil_cfb_ord_na <- FRRN_Zsoil_cfb_fam_na %>% select(-rDNA_m) %>% left_join(FRRN_ord_tab, by = c("Order" = "ord")) %>% 
  filter(is.na(rDNA_m))


# cla
FRRN_Zsoil_cfb_cla <- FRRN_Zsoil_cfb_ord_na %>% select(-rDNA_m) %>% left_join(FRRN_cla_tab, by = c("Class" = "cla")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_cla

FRRN_Zsoil_cfb_cla_na <- FRRN_Zsoil_cfb_ord_na %>% select(-rDNA_m) %>% left_join(FRRN_cla_tab, by = c("Class" = "cla")) %>% 
  filter(is.na(rDNA_m))


# phy
FRRN_Zsoil_cfb_phy <- FRRN_Zsoil_cfb_cla_na %>% select(-rDNA_m) %>% left_join(FRRN_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(!is.na(rDNA_m))
FRRN_Zsoil_cfb_phy

FRRN_Zsoil_cfb_phy_na <- FRRN_Zsoil_cfb_cla_na %>% select(-rDNA_m) %>% left_join(FRRN_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(is.na(rDNA_m))
#FRRN_Zsoil_cfb_phy_na %>% dim()
#FRRN_Zsoil_cfb_phy_na %>% view()

# FRRN Zsoil anno --------------------------
FRRN_Zsoil_cfb_anno <- rbind(FRRN_Zsoil_cfb_phy, FRRN_Zsoil_cfb_cla, FRRN_Zsoil_cfb_ord, 
                             FRRN_Zsoil_cfb_fam, FRRN_Zsoil_cfb_gen, FRRN_Zsoil_cfb_spc)
dim(FRRN_Zsoil_cfb_anno)
# done


# calcu ---------- 
colnames(FRRN_Zsoil_cfb_anno)

Z_soil_otutab_cfb %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min()
# min 5276

Z_soil_otutab_cfb_r <- rrarefy(Z_soil_otutab_cfb %>% select(-X) %>% t() , 5276) %>% t() %>% as.data.frame()
Z_soil_otutab_cfb_r %>% apply(sum, MARGIN = 2)


colnames(Z_soil_otutab_cfb_r) 

Zsoil_otutab_cfb_per <- apply(Z_soil_otutab_cfb_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb$X, .before = X301)

Zsoil_otutab_cfb_per1 <- 
  Zsoil_otutab_cfb_per %>% 
  left_join(FRRN_Zsoil_cfb_anno %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(Zsoil_otutab_cfb_per1 %>% select(-OTU_ID, -rDNA_m, -Lifestyle2))
sample_id



rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- Zsoil_otutab_cfb_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

Zsoil_cfb_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
Zsoil_cfb_rDNAm

Zsoil_cfb_rDNAm1 <- Zsoil_cfb_rDNAm %>% left_join(Z_soil_meta_cfb, by = c("sample_id" = "otu_sample"))


Zsoil_cfb_rDNAm2 <- Zsoil_cfb_rDNAm1 %>% select(sample_id, rDNAm, Latitude1, Longitude1, Altitude1,
                                                TN, TC, ACa, AFe, AK, AMg, TP, C.N, C.P, N.P, soil_D) %>% 
  pivot_longer(names_to = "Factors", values_to = "Value", cols = 3:16)





#view(Zsoil_cfb_rDNAm2)
colnames(Zsoil_cfb_rDNAm2)[3] <- "group"

Zsoil_cbf_envFactor <- Zsoil_cfb_rDNAm2$group %>% unique()


Zsoil_env_cor_subdata <- 
  subplot_data_corr(DF = Zsoil_cfb_rDNAm2, Y_val = "rDNAm", X_val = "Value",
                    GROUP_list = Zsoil_cbf_envFactor, METHOD = "spearman", adj_METHOD = "fdr")



# removed figs
# Zsoil_cfb_rDNAm_env <- 
#   ggplot(Zsoil_cfb_rDNAm2, aes(Value, rDNAm)) +
#   # geom_jitter(size = 1.4, alpha = 0.4, width = 0.2) +
#   geom_hex() +
#   geom_text(data = Zsoil_env_cor_subdata,
#             aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
#             parse = T,
#             colour = "blue",
#             size = 4.5,
#             vjust = 2) + 
#   geom_smooth(colour = "red") +
#   # scale_colour_nejm() +
#   facet_wrap(~ group, scales = "free") +
#   scale_fill_gradient2(low = "#3B4992FF", high = "#008B45FF", mid = "white", midpoint = 5) +
#   guides(fill = guide_colorbar(override.aes = list(size = 5, alpha = 1),
#                                nrow = 1,
#                                direction = "horizontal")) +
#   labs(x = "Environmental factors", y = "Community-weighted rDNA copy number") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         # legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12, face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "inside",
#         legend.position.inside = c(0.75, 0.12),
#         aspect.ratio = 1)
# Zsoil_cfb_rDNAm_env
# 
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS19_pdf <- str_c("figS19_", "Zsoil_cfb_env_", tm, ".pdf", sep = "")
# figS19_jpg <- str_c("figS19_", "Zsoil_cfb_env_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S19_pdf <- str_c(fig_path, figS19_pdf)
# fig_fullpath_S19_jpg <- str_c(fig_path, figS19_jpg)
# 
# ggsave(fig_fullpath_S19_pdf, Zsoil_cfb_rDNAm_env, width = 11.5, height = 11.8)
# ggsave(fig_fullpath_S19_jpg, Zsoil_cfb_rDNAm_env, width = 11.5, height = 11.8)


# sub otu table
FRRN_Zsoil_cfb_anno_ECM <- FRRN_Zsoil_cfb_anno %>% filter(Lifestyle2 == "EcM")
Z_soil_otutab_cfb_ECM <- Z_soil_otutab_cfb %>% filter(X %in% FRRN_Zsoil_cfb_anno_ECM$OTU)

FRRN_Zsoil_cfb_anno_PP <- FRRN_Zsoil_cfb_anno %>% filter(Lifestyle2 == "Plant pathogen")
Z_soil_otutab_cfb_PP <- Z_soil_otutab_cfb %>% filter(X %in% FRRN_Zsoil_cfb_anno_PP$OTU)

FRRN_Zsoil_cfb_anno_Sap <- FRRN_Zsoil_cfb_anno %>% filter(Lifestyle2 == "Saprotroph")
Z_soil_otutab_cfb_Sap <- Z_soil_otutab_cfb %>% filter(X %in% FRRN_Zsoil_cfb_anno_Sap$OTU)


Z_soil_otutab_cfb_ECM %>% select(-X) %>% apply(sum, MARGIN = 2) %>% view()
Z_soil_otutab_cfb_ECM$X %>% length()

# ECM calcu -----------------
Z_soil_otutab_cfb_ECM_r <- rrarefy(Z_soil_otutab_cfb_ECM %>% select(-X) %>% t() , 257) %>% t() %>% as.data.frame()
Z_soil_otutab_cfb_ECM_r %>% apply(sum, MARGIN = 2)

colnames(Z_soil_otutab_cfb_ECM_r) 

Z_soil_otutab_cfb_ECM_per <- apply(Z_soil_otutab_cfb_ECM_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb_ECM$X, .before = X301)

Z_soil_otutab_cfb_ECM_per1 <- 
  Z_soil_otutab_cfb_ECM_per %>% 
  left_join(FRRN_Zsoil_cfb_anno_ECM %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(Z_soil_otutab_cfb_ECM_per1 %>% select(-OTU_ID, -rDNA_m, -Lifestyle2))
sample_id


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- Z_soil_otutab_cfb_ECM_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

Zsoil_cfb_ECM_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
Zsoil_cfb_ECM_rDNAm

Zsoil_cfb_ECM_rDNAm1 <- Zsoil_cfb_ECM_rDNAm %>% left_join(Z_soil_meta_cfb, by = c("sample_id" = "otu_sample"))

#Zsoil_cfb_ECM_rDNAm1 %>% select(Latitude1, Latitude) %>% view()

Zsoil_cfb_ECM_rDNAm2 <- Zsoil_cfb_ECM_rDNAm1 %>% select(sample_id, rDNAm, Latitude1, Longitude1, Altitude1,
                                                TN, TC, ACa, AFe, AK, AMg, TP, C.N, C.P, N.P, soil_D) %>% 
  pivot_longer(names_to = "Factors", values_to = "Value", cols = 3:16)



# supp
Zsoil_cfb_ECM_rDNAm2

# data process
Zsoil_cfb_ECM_rDNAm3 <- Zsoil_cfb_ECM_rDNAm2
colnames(Zsoil_cfb_ECM_rDNAm3)[3] <- "group"

Zsoil_cfb_ECM_group_list <- Zsoil_cfb_ECM_rDNAm3$group %>% unique()

Zsoil_cfb_ECM_cor_subdata <- 
  subplot_data_corr(DF = Zsoil_cfb_ECM_rDNAm3, Y_val = "rDNAm", X_val = "Value",
                    GROUP_list = Zsoil_cfb_ECM_group_list, METHOD = "spearman", adj_METHOD = "bonferroni")
Zsoil_cfb_ECM_cor_subdata

# removed figs
# Zsoil_cfb_ECM_rDNAm_env <- 
#   ggplot(Zsoil_cfb_ECM_rDNAm3, aes(Value, rDNAm)) +
#   # geom_jitter(size = 1.4, alpha = 0.4, width = 0.2) +
#   geom_hex() +
#   geom_text(data = Zsoil_cfb_ECM_cor_subdata,
#             aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
#             parse = T,
#             colour = "blue",
#             size = 4.5,
#             vjust = 2) + 
#   geom_smooth(colour = "red") +
#   # scale_colour_nejm() +
#   facet_wrap(~ group, scales = "free") +
#   scale_fill_gradient2(low = "#3B4992FF", high = "#008B45FF", mid = "white", midpoint = 6) +
#   guides(fill = guide_colorbar(override.aes = list(size = 5, alpha = 1),
#                                nrow = 1,
#                                direction = "horizontal")) +
#   labs(x = "Environmental factors", y = "Community-weighted rDNA copy number of ectomycorrhizal fungi") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         # legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12, face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "inside",
#         legend.position.inside = c(0.75, 0.12),
#         aspect.ratio = 1)
# Zsoil_cfb_ECM_rDNAm_env
# 
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS20_pdf <- str_c("figS20_", "Zsoil_ECM_env_", tm, ".pdf", sep = "")
# figS20_jpg <- str_c("figS20_", "Zsoil_ECM_env_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S20_pdf <- str_c(fig_path, figS20_pdf)
# fig_fullpath_S20_jpg <- str_c(fig_path, figS20_jpg)
# 
# ggsave(fig_fullpath_S20_pdf, Zsoil_cfb_ECM_rDNAm_env, width = 11.5, height = 11.8)
# ggsave(fig_fullpath_S20_jpg, Zsoil_cfb_ECM_rDNAm_env, width = 11.5, height = 11.8)


# plant pathogen
PP_ind <- Z_soil_otutab_cfb_PP %>% select(-X) %>% apply(sum, MARGIN = 2) %>% as.data.frame()

colnames(PP_ind) <- "Freq"


PP_ind <- PP_ind %>% mutate(sample_id = rownames(PP_ind)) %>% filter(Freq > 300) %>%
  select(sample_id) %>% pull()

Z_soil_otutab_cfb_PP1 <- Z_soil_otutab_cfb_PP %>% select(X, all_of(PP_ind))
Z_soil_otutab_cfb_PP1 %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min()

# PP calcu -----------------
Z_soil_otutab_cfb_PP_r <- rrarefy(Z_soil_otutab_cfb_PP1 %>% select(-X) %>% t() , 308) %>% t() %>% as.data.frame()
Z_soil_otutab_cfb_PP_r %>% apply(sum, MARGIN = 2)

colnames(Z_soil_otutab_cfb_PP_r) 

Z_soil_otutab_cfb_PP_per <- apply(Z_soil_otutab_cfb_PP_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb_PP1$X, .before = X302)

Z_soil_otutab_cfb_PP_per1 <- 
  Z_soil_otutab_cfb_PP_per %>% 
  left_join(FRRN_Zsoil_cfb_anno_PP %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(Z_soil_otutab_cfb_PP_per1 %>% select(-OTU_ID, -rDNA_m, -Lifestyle2))
sample_id



rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- Z_soil_otutab_cfb_PP_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

Zsoil_cfb_PP_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
Zsoil_cfb_PP_rDNAm

Zsoil_cfb_PP_rDNAm1 <- Zsoil_cfb_PP_rDNAm %>% left_join(Z_soil_meta_cfb, by = c("sample_id" = "otu_sample"))

# Zsoil_cfb_PP_rDNAm1 %>% select(Latitude1, Latitude) %>% view()

Zsoil_cfb_PP_rDNAm2 <- Zsoil_cfb_PP_rDNAm1 %>% select(sample_id, rDNAm, Latitude1, Longitude1, Altitude1,
                                                        TN, TC, ACa, AFe, AK, AMg, TP, C.N, C.P, N.P, soil_D) %>% 
  pivot_longer(names_to = "Factors", values_to = "Value", cols = 3:16)


# supp
Zsoil_cfb_PP_rDNAm2

# data process
Zsoil_cfb_PP_rDNAm3 <- Zsoil_cfb_PP_rDNAm2
colnames(Zsoil_cfb_PP_rDNAm3)[3] <- "group"

Zsoil_cfb_PP_group_list <- Zsoil_cfb_PP_rDNAm3$group %>% unique()

Zsoil_cfb_PP_cor_subdata <- 
  subplot_data_corr(DF = Zsoil_cfb_PP_rDNAm3, Y_val = "rDNAm", X_val = "Value",
                    GROUP_list = Zsoil_cfb_PP_group_list, METHOD = "spearman", adj_METHOD = "bonferroni")
Zsoil_cfb_PP_cor_subdata

# removed figs
# Zsoil_cfb_PP_rDNAm_env <- 
#   ggplot(Zsoil_cfb_PP_rDNAm3, aes(Value, rDNAm)) +
#   # geom_jitter(size = 1.4, alpha = 0.4, width = 0.2) +
#   geom_hex() +
#   geom_text(data = Zsoil_cfb_PP_cor_subdata,
#             aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
#             parse = T,
#             colour = "blue",
#             size = 4.5,
#             vjust = 2) + 
#   geom_smooth(colour = "red") +
#   # scale_colour_nejm() +
#   facet_wrap(~ group, scales = "free") +
#   scale_fill_gradient2(low = "#3B4992FF", high = "#008B45FF", mid = "white", midpoint = 5) +
#   guides(fill = guide_colorbar(override.aes = list(size = 5, alpha = 1),
#                                nrow = 1,
#                                direction = "horizontal")) +
#   labs(x = "Environmental factors", y = "Community-weighted rDNA copy number of plant pathogen fungi") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         # legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12, face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "inside",
#         legend.position.inside = c(0.75, 0.12),
#         aspect.ratio = 1)
# Zsoil_cfb_PP_rDNAm_env
# 
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS21_pdf <- str_c("figS21_", "Zsoil_PP_env_", tm, ".pdf", sep = "")
# figS21_jpg <- str_c("figS21_", "Zsoil_PP_env_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S21_pdf <- str_c(fig_path, figS21_pdf)
# fig_fullpath_S21_jpg <- str_c(fig_path, figS21_jpg)
# 
# ggsave(fig_fullpath_S21_pdf, Zsoil_cfb_PP_rDNAm_env, width = 11.5, height = 11.8)
# ggsave(fig_fullpath_S21_jpg, Zsoil_cfb_PP_rDNAm_env, width = 11.5, height = 11.8)



# Sap calcu -----------------
Z_soil_otutab_cfb_Sap %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min()

Z_soil_otutab_cfb_Sap_r <- rrarefy(Z_soil_otutab_cfb_Sap %>% select(-X) %>% t() , 2656) %>% t() %>% as.data.frame()

Z_soil_otutab_cfb_Sap_r %>% apply(sum, MARGIN = 2)


Z_soil_otutab_cfb_Sap_per <- apply(Z_soil_otutab_cfb_Sap_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb_Sap$X, .before = X301)

Z_soil_otutab_cfb_Sap_per1 <- 
  Z_soil_otutab_cfb_Sap_per %>% 
  left_join(FRRN_Zsoil_cfb_anno_Sap %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(Z_soil_otutab_cfb_Sap_per1 %>% select(-OTU_ID, -rDNA_m, -Lifestyle2))
sample_id



rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- Z_soil_otutab_cfb_Sap_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
  tmp_val <- tmp_df[, 2] * tmp_df[, 3]
  
  tmp_df1 <- tmp_df %>% 
    mutate(tmp_val = tmp_val)
  
  na_rel_per <- tmp_df1 %>% filter(is.na(rDNA_m)) %>% select(2) %>% pull() %>% sum()
  rDNA_rel_per <- tmp_df1 %>% filter(!is.na(rDNA_m)) %>% select(4) %>% pull() %>% sum()
  
  rDNA_rel_per_rlt <- rDNA_rel_per / (1 - na_rel_per)
  
  rlt_df <- data.frame(
    sample_id = sample_id[1],
    rDNAm = rDNA_rel_per_rlt
  )
  
  return(rlt_df)
  
}

Zsoil_cfb_Sap_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
Zsoil_cfb_Sap_rDNAm

Zsoil_cfb_Sap_rDNAm1 <- Zsoil_cfb_Sap_rDNAm %>% left_join(Z_soil_meta_cfb, by = c("sample_id" = "otu_sample"))

#Zsoil_cfb_Sap_rDNAm1 %>% select(Latitude1, Latitude) %>% view()

Zsoil_cfb_Sap_rDNAm2 <- Zsoil_cfb_Sap_rDNAm1 %>% select(sample_id, rDNAm, Latitude1, Longitude1, Altitude1,
                                                        TN, TC, ACa, AFe, AK, AMg, TP, C.N, C.P, N.P, soil_D) %>% 
  pivot_longer(names_to = "Factors", values_to = "Value", cols = 3:16)


# supp
Zsoil_cfb_Sap_rDNAm2

# data process
Zsoil_cfb_Sap_rDNAm3 <- Zsoil_cfb_Sap_rDNAm2
colnames(Zsoil_cfb_Sap_rDNAm3)[3] <- "group"

Zsoil_cfb_Sap_group_list <- Zsoil_cfb_Sap_rDNAm3$group %>% unique()

Zsoil_cfb_Sap_cor_subdata <- 
  subplot_data_corr(DF = Zsoil_cfb_Sap_rDNAm3, Y_val = "rDNAm", X_val = "Value",
                    GROUP_list = Zsoil_cfb_Sap_group_list, METHOD = "spearman", adj_METHOD = "bonferroni")
Zsoil_cfb_Sap_cor_subdata

# removed figs
# Zsoil_cfb_Sap_rDNAm_env <- 
#   ggplot(Zsoil_cfb_Sap_rDNAm3, aes(Value, rDNAm)) +
#   # geom_jitter(size = 1.4, alpha = 0.4, width = 0.2) +
#   geom_hex() +
#   geom_text(data = Zsoil_cfb_Sap_cor_subdata,
#             aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
#             parse = T,
#             colour = "blue",
#             size = 4.5,
#             vjust = 2) + 
#   geom_smooth(colour = "red") +
#   # scale_colour_nejm() +
#   facet_wrap(~ group, scales = "free") +
#   scale_fill_gradient2(low = "#3B4992FF", high = "#008B45FF", mid = "white", midpoint = 6) +
#   guides(fill = guide_colorbar(override.aes = list(size = 5, alpha = 1),
#                                nrow = 1,
#                                direction = "horizontal")) +
#   labs(x = "Environmental factors", y = "Community-weighted rDNA copy number of Saprotroph fungi") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         # legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12, face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "inside",
#         legend.position.inside = c(0.75, 0.12),
#         aspect.ratio = 1)
# Zsoil_cfb_Sap_rDNAm_env
# 
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS22_pdf <- str_c("figS22_", "Zsoil_Sap_env_", tm, ".pdf", sep = "")
# figS22_jpg <- str_c("figS22_", "Zsoil_Sap_env_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S22_pdf <- str_c(fig_path, figS22_pdf)
# fig_fullpath_S22_jpg <- str_c(fig_path, figS22_jpg)
# 
# ggsave(fig_fullpath_S22_pdf, Zsoil_cfb_Sap_rDNAm_env, width = 11.5, height = 11.8)
# ggsave(fig_fullpath_S22_jpg, Zsoil_cfb_Sap_rDNAm_env, width = 11.5, height = 11.8)



# supp
Zsoil_cfb_EPS_rDNAm2 <- rbind(Zsoil_cfb_ECM_rDNAm3 %>% mutate(Guild = "Ectomycorrhizal fungi"), 
                              Zsoil_cfb_PP_rDNAm3 %>% mutate(Guild = "Plant pathogen"), 
                              Zsoil_cfb_Sap_rDNAm3 %>% mutate(Guild = "Saprotroph fungi"))
EPS_x_midpoint <- 
  Zsoil_cfb_EPS_rDNAm2 %>% group_by(group) %>% summarise(
    x_midpoint = (range(Value)[2] - range(Value)[1]) * 0.05 + range(Value)[1]
  )




Zsoil_cfb_EPS_cor_subdata <- rbind(Zsoil_cfb_ECM_cor_subdata %>% mutate(Guild = "Ectomycorrhizal fungi"),
                                   Zsoil_cfb_PP_cor_subdata %>% mutate(Guild = "Plant pathogen"),
                                   Zsoil_cfb_Sap_cor_subdata %>% mutate(Guild = "Saprotroph fungi")) %>% 
  left_join(EPS_x_midpoint, by = "group")


Zsoil_cfb_EPS_rDNAm_LAT <- 
  ggplot(Zsoil_cfb_EPS_rDNAm2 %>% filter(group == "Latitude1"), aes(Value, rDNAm, colour = Guild)) +
  geom_jitter(size = 2, alpha = 0.4, width = 0.4) +
  geom_text(data = Zsoil_cfb_EPS_cor_subdata %>% 
              filter(Guild == "Ectomycorrhizal fungi") %>% 
              filter(group == "Latitude1"),
            aes(x = x_midpoint, y = Inf, label = rlt_sig),
            parse = T,
            colour = "#3B4992FF",
            size = 4,
            vjust = 3, hjust = 0) + 
  geom_text(data = Zsoil_cfb_EPS_cor_subdata %>% 
              filter(Guild == "Plant pathogen") %>% 
              filter(group == "Latitude1"),
            aes(x = x_midpoint, y = Inf, label = rlt_sig),
            parse = T,
            colour = "#EE0000FF",
            size = 4,
            vjust = 4.5, hjust = 0) + 
  geom_text(data = Zsoil_cfb_EPS_cor_subdata %>% 
              filter(Guild == "Saprotroph fungi") %>% 
              filter(group == "Latitude1"),
            aes(x = x_midpoint, y = Inf, label = rlt_sig),
            parse = T,
            colour = "#008B45FF",
            size = 4,
            vjust = 6, hjust = 0) + 
  geom_smooth(alpha = 0.2, show.legend = F, method = "lm") +
  scale_colour_aaas() +
  # facet_wrap(~ group, scales = "free") +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
                               ncol = 1,
                               direction = "vertical",
                               title = NULL)) +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.87))
Zsoil_cfb_EPS_rDNAm_LAT

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7b_pdf <- str_c("figS7b_", "Zsoil_EPS_LAT_", tm, ".pdf", sep = "")
figS7b_jpg <- str_c("figS7b_", "Zsoil_EPS_LAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7b_pdf <- str_c(fig_path, figS7b_pdf)
fig_fullpath_S7b_jpg <- str_c(fig_path, figS7b_jpg)

ggsave(fig_fullpath_S7b_pdf, Zsoil_cfb_EPS_rDNAm_LAT, width = 5.36, height = 5.05)
ggsave(fig_fullpath_S7b_jpg, Zsoil_cfb_EPS_rDNAm_LAT, width = 5.36, height = 5.05)

# figS7
figS7 <- XC_ROOT_EMF_rDNAm_LAT_full + Zsoil_cfb_EPS_rDNAm_LAT
figS7

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7_pdf <- str_c("figS7_", "XC_ROOT_Zsoil_LAT_", tm, ".pdf", sep = "")
figS7_jpg <- str_c("figS7_", "XC_ROOT_Zsoil_LAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7_pdf <- str_c(fig_path, figS7_pdf)
fig_fullpath_S7_jpg <- str_c(fig_path, figS7_jpg)

ggsave(fig_fullpath_S7_pdf, figS7, width = 10.4, height = 5.06)
ggsave(fig_fullpath_S7_jpg, figS7, width = 10.4, height = 5.06)
# Saving 10.4 x 5.06 in image


Zsoil_cfb_EPS_rDNAm2_Lat <- 
  Zsoil_cfb_EPS_rDNAm2 %>% filter(group == "Latitude1") %>% 
  mutate(group = str_sub(group, end = -2))


lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm2 %>% filter(Guild == "Ectomycorrhizal fungi") %>% filter(group == "Latitude1")) %>% 
  summary() %>% tidy()

lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm2 %>% filter(Guild == "Plant pathogen") %>% filter(group == "Latitude1")) %>% 
  summary() %>% tidy()

lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm2 %>% filter(Guild == "Saprotroph fungi") %>% filter(group == "Latitude1")) %>% 
  summary() %>% tidy()


# combine XC and Zheng, root and soil ------------------
# XCroot_rDNAm_lat
# XCroot_rDNAm2

##### may have bugs #####
XC_ROOT_c_tmp <- XC_ROOT_rDNAm1 %>% select(sample_id, rDNAm, Latitude) %>% 
  mutate(Compartment = "Root")

colnames(XC_ROOT_c_tmp)


Zsoil_cfb_rDNAm1_comb_tmp <- Zsoil_cfb_rDNAm1 %>% select(sample_id, rDNAm, Latitude1) %>% 
  mutate(Compartment = "Soil")

colnames(Zsoil_cfb_rDNAm1_comb_tmp)[3] <- "Latitude"

RS_rDNAm <- rbind(XC_ROOT_c_tmp, Zsoil_cfb_rDNAm1_comb_tmp)
RS_rDNAm



colnames(RS_rDNAm)[4] <- "group"

RS_cor_subdata <- 
  subplot_data_corr(DF = RS_rDNAm, Y_val = "rDNAm", X_val = "Latitude",
                    GROUP_list = c("Root", "Soil"), METHOD = "spearman", adj_METHOD = "fdr")
RS_cor_subdata


# lm(rDNAm ~ Latitude, data = XCroot_rDNAm1) %>% tidy()
# # italic(r) == 0.198 ~~ italic(p) == 1.60e-14
# 
# lm(rDNAm ~ Latitude1, data = Zsoil_cfb_rDNAm1) %>% tidy()
# # italic(r) == 0.135 ~~ italic(p) == 0.037

# fix(RS_cor_subdata)

# fig3a, RS rDNAm ~ Lat --------------------------
fig3a_RS_rDNAm_lat <- 
  ggplot(RS_rDNAm, aes(Latitude, rDNAm)) +
  geom_jitter(aes(colour = group), size = 1.8, alpha = 0.4, width = 0.4) +
  geom_smooth(aes(colour = group), method = "lm", alpha = 0.1) +
  geom_text(data = RS_cor_subdata,
            aes(x = anno_x1, y = anno_y1,
                label = rlt_sig, 
                colour = group),
            parse = T,
            show.legend = F,
            size = 4.5,
            vjust = 0,
            hjust = 0) +
  guides(colour = guide_legend(
    title = NULL,
    nrow = 1
  )) +
  scale_colour_manual(values = c("blue",
                                 "red")) +
  scale_y_continuous(limits = c(20, 180)) +
  labs(y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", face = "bold", size = 15),
    legend.title = element_text(face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.27, 0.93),
    aspect.ratio = 1
  )
fig3a_RS_rDNAm_lat

tm <- now() %>% str_split_i(pattern = " ", 1)
fig3a_pdf <- str_c("fig3a_", "RS_rDNAm_lat_", tm, ".pdf", sep = "")
fig3a_jpg <- str_c("fig3a_", "RS_rDNAm_lat_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_3a_pdf <- str_c(fig_path, fig3a_pdf)
fig_fullpath_3a_jpg <- str_c(fig_path, fig3a_jpg)


ggsave(fig_fullpath_3a_pdf, fig3a_RS_rDNAm_lat, width = 6.1, height = 5.8)
ggsave(fig_fullpath_3a_jpg, fig3a_RS_rDNAm_lat, width = 6.1, height = 5.8)


# EPICON lm result ----------------

lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Treatment == subgrp)
    
  } else {
    
    df_tmp <- df
    
  }
  
  df_tmp <- df_tmp %>% drop_na()
  
  sub_lmR <- lm(df_tmp[[yval]] ~ df_tmp[[xval]])
  lmR_sumy <- summary(sub_lmR)
  
  if(lmR_sumy$coefficients[2] > 0) {
    lmR_r <- round(sqrt(lmR_sumy$r.squared), 3) %>% signif(., 3)
  } else {
    lmR_r <- round(sqrt(lmR_sumy$r.squared), 3) %>% signif(., 3) * -1
  }
  
  lmR_P_tmpf <- lmR_sumy$fstatistic
  lmR_P_tmp <- pf(lmR_P_tmpf[1], lmR_P_tmpf[2], lmR_P_tmpf[3], lower.tail = F)
  
  lmR_p <- round(lmR_P_tmp, 3) %>% signif(., 3)
  
  lmR_slope <- lmR_sumy$coefficients[2]
  
  if(is.na(lmR_p)) {
    
    Rp_sig <- "NA"
    
  } else if(lmR_p < 0.001) {
    
    Rp_sig <- "***"
    
  } else if(lmR_p <= 0.01) {
    
    Rp_sig <- "**"
    
  } else if(lmR_p <= 0.05){
    
    Rp_sig <- "*"
    
  } else if(lmR_p <= 1){
    
    Rp_sig <- "NS"
    
  }
  
  if(trans != "no_trans") {
    
    df_rlt <- data.frame(
      Treatment = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_P_tmp,
      slope = lmR_slope
    )
    
  } else {
    
    df_rlt <- data.frame(
      Treatment = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_P_tmp,
      slope = lmR_slope
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


# EPICON_leaf_lm <- 
#   map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
#           .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_env1 %>% filter(Habitat == "Leaf"),
#                              yval = "rDNAm", xval = "timepoint",
#                              subgrp = ., trans = "no_trans"))
# EPICON_leaf_lm
# 
# EPICON_rDNAm_env1 %>% filter(Habitat == "Leaf") %>% group_by(Treatment, Timepiont) %>% summarise(
#   rDNAm_med = median(rDNAm),
#   rDNAm_sd = sd(rDNAm)
# ) %>% view()
# 
# 
# EPICON_soil_lm <- 
#   map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
#           .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_env1 %>% filter(Habitat == "Soil"),
#                              yval = "rDNAm", xval = "timepoint",
#                              subgrp = ., trans = "no_trans"))
# EPICON_soil_lm
# 
# EPICON_rDNAm_env1 %>% filter(Habitat == "Soil") %>% group_by(Treatment, Timepiont) %>% summarise(
#   rDNAm_med = median(rDNAm),
#   rDNAm_sd = sd(rDNAm)
# ) %>% view()
# 
# 
# EPICON_rDNAm_env1 %>% filter(Habitat == "Root") %>% group_by(Treatment, Timepiont) %>% summarise(
#   rDNAm_med = median(rDNAm),
#   rDNAm_sd = sd(rDNAm)
# ) %>% view()
# 
# 
# EPICON_rDNAm_env1 %>% filter(Habitat == "Rhizosphere") %>% group_by(Treatment, Timepiont) %>% summarise(
#   rDNAm_med = median(rDNAm),
#   rDNAm_sd = sd(rDNAm)
# ) %>% view()
# 
# FRRN_phy_tab
# EPICON_rDNAm_env1 %>% dim()

# done.
