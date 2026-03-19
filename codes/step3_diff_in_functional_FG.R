
# step3: differences among fungal traits in FunGuild Database.
# fig2a, fig2b, fig2c, figS3, figS4
# Qiushi-Li, IM-CAS
# 2024.11.30


# packages
library(tidyverse)
library(readxl)

# install.packages("rcompanion")
library(rcompanion)
library(agricolae)
# devtools::install_github('erocoar/gghalves')
library(gghalves)
library(ggsci)
library(patchwork)
library(broom)

library(writexl)

# load RData ----------
#load("rCNV_version_1.0_20250219.RData")

# FunGuild database --------------------------
FunGuild_database <- read_csv("./database/FunGuild_v2024.csv")
colnames(FunGuild_database)

Tax_Lev <- c(13, 9, 7, 3)
Tax_Ind <- c("gen", "fam", "ord", "phy")
FG_index_df <- data.frame(
  Tax_Lev,
  Tax_Ind
)
FG_index_df

# rCNV total tmp table
rCNV_rlt_taxa_spl

rCNV_taxa_FG <- rCNV_rlt_taxa_spl %>%
  select(project, Both_ITS_LSU, Name, phy, cla, ord, fam, Gen)
colnames(rCNV_taxa_FG)[8] <- "gen"

rCNV_taxa_FG <- rCNV_taxa_FG %>% mutate(
  phy = str_sub(phy, start = 3),
  cla = str_sub(cla, start = 3),
  ord = str_sub(ord, start = 3),
  fam = str_sub(fam, start = 3),
  gen = str_sub(gen, start = 3)
)


# Full align ----
FunGuild_database1 <- FunGuild_database %>% select(-Notes, -`Citation/Source`)
colnames(FunGuild_database1)

# view(FunGuild_database1)

FG_full_tmp_gen <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[1]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_gen)[1] <- FG_index_df$Tax_Ind[1]

rCNV_taxa_FG_full_gen <- rCNV_taxa_FG %>% 
  left_join(FG_full_tmp_gen)

rCNV_taxa_FG_full_gen_na <- rCNV_taxa_FG_full_gen %>% filter(is.na(Trophic_Mode))
rCNV_taxa_FG_full_gen_na

rCNV_taxa_FG_full_gen_anno <- rCNV_taxa_FG_full_gen %>% filter(!is.na(Trophic_Mode))

FG_full_tmp_fam <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[2]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_fam)[1] <- FG_index_df$Tax_Ind[2]

rCNV_taxa_FG_full_fam <- rCNV_taxa_FG_full_gen_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_fam)

rCNV_taxa_FG_full_fam_na <- rCNV_taxa_FG_full_fam %>% filter(is.na(Trophic_Mode))
rCNV_taxa_FG_full_fam_na

rCNV_taxa_FG_full_fam_anno <- rCNV_taxa_FG_full_fam %>% filter(!is.na(Trophic_Mode))

FG_full_tmp_ord <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[3]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_ord)[1] <- FG_index_df$Tax_Ind[3]

rCNV_taxa_FG_full_ord <- rCNV_taxa_FG_full_fam_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_ord)

rCNV_taxa_FG_full_ord_na <- rCNV_taxa_FG_full_ord %>% filter(is.na(Trophic_Mode))
rCNV_taxa_FG_full_ord_na

rCNV_taxa_FG_full_ord_anno <- rCNV_taxa_FG_full_ord %>% filter(!is.na(Trophic_Mode))
rCNV_taxa_FG_full_ord_anno


FG_full_tmp_phy <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[4]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_phy)[1] <- FG_index_df$Tax_Ind[4]

rCNV_taxa_FG_full_phy <- rCNV_taxa_FG_full_ord_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_phy)

rCNV_taxa_FG_full_phy_na <- rCNV_taxa_FG_full_phy %>% filter(is.na(Trophic_Mode))
rCNV_taxa_FG_full_phy_na

rCNV_taxa_FG_full_phy_anno <- rCNV_taxa_FG_full_phy %>% filter(!is.na(Trophic_Mode))
rCNV_taxa_FG_full_phy_anno


rCNV_taxa_FGanno_full <- rbind(rCNV_taxa_FG_full_gen_anno, rCNV_taxa_FG_full_fam_anno, rCNV_taxa_FG_full_ord_anno, rCNV_taxa_FG_full_phy_anno, 
                               rCNV_taxa_FG_full_phy_na)
rCNV_taxa_FGanno_full

rCNV_taxa_FGanno_full %>% filter(is.na(Trophic_Mode))
# 9 unannotated


rCNV_taxa_FGanno_full <- rCNV_taxa_FGanno_full %>% filter(!is.na(Trophic_Mode))
dim(rCNV_taxa_FGanno_full)
# write_xlsx(rCNV_taxa_FGanno_full, "rCNV_taxa_FunGuild_20250129.xlsx")


# may needed adjusted
rCNV_taxa_FGanno_full1 <- 
  rCNV_taxa_FGanno_full %>% 
  mutate(Guild1 = str_extract(rCNV_taxa_FGanno_full$Guild, "\\|([^\\|]+)\\|"), .before = Growth_Morphology) %>% 
  mutate(Guild1 = if_else(!is.na(Guild1), str_sub(Guild1, start = 2, end = -2), Guild))

# write_xlsx(rCNV_taxa_FGanno_full1, "rCNV_taxa_FG_20250129.xlsx")
rCNV_taxa_FGanno_full1$Guild1 %>% table()

rCNV_FG_subgrp_ind <- c("Plant Pathogen",
                        "Dung Saprotroph", "Pollen Saprotroph", "Wood Saprotroph", "Plant Saprotroph", "Undefined Saprotroph",
                        "Endophyte", "Epiphyte", "Arbuscular Mycorrhizal", "Ectomycorrhizal")

rCNV_taxa_FGanno_subgrp <- 
  rCNV_taxa_FGanno_full1 %>% 
  filter(Guild1 %in% rCNV_FG_subgrp_ind) %>% 
  mutate(Guild2 = if_else(str_detect(Guild1, pattern = "Saprotroph"), "Saprotroph Fungi", Guild1), .before = Growth_Morphology)

rCNV_taxa_FGanno_subgrp1 <- 
  rCNV_taxa_FGanno_subgrp %>% 
  filter(Guild2 == "Plant Pathogen" | Guild2 == "Saprotroph Fungi" | Guild2 == "Ectomycorrhizal") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota") %>% 
  mutate(Guild2 = if_else(str_detect(Guild2, pattern = "Plant Pathogen"), "Plant Pathogen Fungi", Guild2)) %>% 
  mutate(Guild2 = if_else(str_detect(Guild2, pattern = "Ectomycorrhizal"), "Ectomycorrhizal Fungi", Guild2))

# fig2a, Asc Bas combine ------------------------
rCNV_taxa_FGanno_subgrp2 <- 
  rCNV_taxa_FGanno_subgrp1 %>% mutate(Guild3 = str_c(phy, Guild2, sep = "-"), .after = Guild2)

rCNV_taxa_FGanno_subgrp2$Guild3 <- factor(rCNV_taxa_FGanno_subgrp2$Guild3, levels = c(
  "Ascomycota-Ectomycorrhizal Fungi", "Ascomycota-Saprotroph Fungi", "Ascomycota-Plant Pathogen Fungi",
  "Basidiomycota-Ectomycorrhizal Fungi", "Basidiomycota-Saprotroph Fungi", "Basidiomycota-Plant Pathogen Fungi"
))


# kru_rlt_FG_subgrp2_subdata1


rCNV_FG_subgrp2_kru <- 
  kruskal(rCNV_taxa_FGanno_subgrp2$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp2$Guild3, p.adj = "fdr")
rCNV_FG_subgrp2_kru 

# kru pvalue --------

rCNV_taxa_FGanno_subgrp2$Guild3 %>% table()
rCNV_taxa_FGanno_subgrp2$Guild2 %>% table()


kruskal.test(Both_ITS_LSU ~ Guild3, data = rCNV_taxa_FGanno_subgrp2 %>% filter(phy == "Ascomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 24.055, df = 1, p-value = 9.363e-07

kruskal.test(Both_ITS_LSU ~ Guild3, data = rCNV_taxa_FGanno_subgrp2 %>% filter(phy == "Basidiomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 3.6203, df = 1, p-value = 0.05708

kruskal.test(Both_ITS_LSU ~ phy, data = rCNV_taxa_FGanno_subgrp2 %>% filter(Guild2 == "Plant Pathogen Fungi"))


rCNV_FG_subgrp2_plot_subdata1 <- 
  rCNV_FG_subgrp2_kru$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups)
rCNV_FG_subgrp2_plot_subdata1

rCNV_FG_subgrp2_maxrCNV <- rCNV_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_subgrp2_plot_subdata1 <- rCNV_FG_subgrp2_plot_subdata1 %>% left_join(rCNV_FG_subgrp2_maxrCNV, by = "Guild3")

rCNV_FG_subgrp2_MeSd <- rCNV_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                                        sd_rDNA = sd(Both_ITS_LSU))
rCNV_FG_subgrp2_plot_subdata1 <- 
  rCNV_FG_subgrp2_plot_subdata1 %>% left_join(rCNV_FG_subgrp2_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

kru_rlt_FG_subgrp2_subdata2 <- rCNV_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% count(Guild3) %>% 
  mutate(lab = str_c("n =", n, sep = " "))


rCNV_AscBas_SRH <- scheirerRayHare(Both_ITS_LSU ~ phy * Guild2, data = rCNV_taxa_FGanno_subgrp2)
rCNV_AscBas_SRH

#fix(rCNV_FG_subgrp2_plot_subdata1)


rCNV_taxa_FGanno_subgrp2 %>% filter(Both_ITS_LSU > 700)

fig2a_rCNV_FG_AscBas <- 
  ggplot(rCNV_taxa_FGanno_subgrp2,
         aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_subgrp2_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 25, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = rCNV_FG_subgrp2_plot_subdata1,
            aes(x = Guild3, y = max_rDNA, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = kru_rlt_FG_subgrp2_subdata2,
            aes(x = Guild3, y = -25, label = lab), colour = "black", hjust = 0.5) +
  geom_segment(x = 3.5, xend = 3.5, y = -20, yend = 580, linetype = 2) +
  annotate(geom = "text", x = 2, y = 510, label = "Ascomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 5, y = 510, label = "Basidiomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 3.5, y = 630, label = expression(Phylum ~~ (P): H == "126.657;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~ 
                                                                 Guild ~~ (G): H == "23.006;" ~~ df == "2;" ~~ italic(p) == "1.010e-05;" ~~
                                                                 "P x G:" ~~ H == "30.811;" ~~ df == "2;" ~~ italic(p) == "2.039e-07"),
           size = 3.8, colour = "black", fontface = "bold") +
  scale_fill_d3(palette = "category20") +
  scale_colour_d3(palette = "category20") +
  scale_x_discrete(labels = c("Ectomycorrhizal\nfungi", "Saprotrophic\nfungi", "Plant pathogenic\nfungi",
                              "Ectomycorrhizal\nfungi", "Saprotrophic\nfungi", "Plant pathogenic\nfungi"),
                   expand = c(0, 0)) +
  scale_y_continuous(limits = c(-50, 680),
                     expand = c(0, 0)) +
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
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5))
fig2a_rCNV_FG_AscBas
# > 700 un shown

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2a_pdf <- str_c("fig2a_", "rCNV_FG_AscBas_", tm, ".pdf", sep = "")
fig2a_jpg <- str_c("fig2a_", "rCNV_FG_AscBas_", tm, ".jpg", sep = "")

ggsave(fig2a_pdf, fig2a_rCNV_FG_AscBas, width = 10.6, height = 4.44)
ggsave(fig2a_jpg, fig2a_rCNV_FG_AscBas, width = 10.6, height = 4.44)

# fig2a, IQR -----------------
rCNV_taxa_FGanno_subgrp2

IQR_outlier <- function(val) {
  
  Q1 <- quantile(val, 0.25)
  Q3 <- quantile(val, 0.75)
  IQR <- Q3 - Q1
  
  # filter
  outliers_ind <- (val < (Q1 - 1.5 * IQR) | val > (Q3 + 1.5 * IQR))
  
  return(outliers_ind)
  
}


rCNV_taxa_FGanno_subgrp2_IQR <- 
  rCNV_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()


rCNV_FG_subgrp2_IQR_kru <- 
  kruskal(rCNV_taxa_FGanno_subgrp2_IQR$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp2_IQR$Guild3, p.adj = "fdr")
rCNV_FG_subgrp2_IQR_kru 

# kru pvalue --------

kruskal.test(Both_ITS_LSU ~ Guild3, data = rCNV_taxa_FGanno_subgrp2_IQR %>% filter(phy == "Ascomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 23.882, df = 1, p-value = 1.024e-06

kruskal.test(Both_ITS_LSU ~ Guild3, data = rCNV_taxa_FGanno_subgrp2_IQR %>% filter(phy == "Basidiomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 2.8651, df = 1, p-value = 0.09052

rCNV_FG_subgrp2_IQR_plot_subdata1 <- 
  rCNV_FG_subgrp2_IQR_kru$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups)
rCNV_FG_subgrp2_IQR_plot_subdata1

rCNV_FG_subgrp2_IQR_maxrCNV <- rCNV_taxa_FGanno_subgrp2_IQR %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_subgrp2_IQR_plot_subdata1 <- rCNV_FG_subgrp2_IQR_plot_subdata1 %>% left_join(rCNV_FG_subgrp2_IQR_maxrCNV, by = "Guild3")

rCNV_FG_subgrp2_IQR_MeSd <- rCNV_taxa_FGanno_subgrp2_IQR %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                                        sd_rDNA = sd(Both_ITS_LSU))
rCNV_FG_subgrp2_IQR_plot_subdata1 <- 
  rCNV_FG_subgrp2_IQR_plot_subdata1 %>% left_join(rCNV_FG_subgrp2_IQR_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

kru_rlt_FG_subgrp2_IQR_subdata2 <- rCNV_taxa_FGanno_subgrp2_IQR %>% group_by(Guild3) %>% count(Guild3) %>% 
  mutate(lab = str_c("n =", n, sep = " "))


rCNV_AscBas_IQR_SRH <- scheirerRayHare(Both_ITS_LSU ~ phy * Guild2, data = rCNV_taxa_FGanno_subgrp2_IQR)
rCNV_AscBas_IQR_SRH
rCNV_AscBas_IQR_SRH$p.value


fig2a_rCNV_FG_AscBas_IQR <- 
  ggplot(rCNV_taxa_FGanno_subgrp2_IQR,
         aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_subgrp2_IQR_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 25, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = rCNV_FG_subgrp2_IQR_plot_subdata1,
            aes(x = Guild3, y = max_rDNA, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = kru_rlt_FG_subgrp2_IQR_subdata2,
            aes(x = Guild3, y = 0, label = lab), vjust = 2, colour = "black", hjust = 0.5) +
  geom_segment(x = 3.5, xend = 3.5, y = -20, yend = 560, linetype = 2) +
  annotate(geom = "text", x = 2, y = 510, label = "Ascomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 5, y = 510, label = "Basidiomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 3.5, y = 585, label = expression(Phylum ~~ (P): H == "133.091;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~ 
                                                                 Guild ~~ (G): H == "16.288;" ~~ df == "2;" ~~ italic(p) == "2.91e-04;" ~~
                                                                 "P x G:" ~~ H == "27.845;" ~~ df == "2;" ~~ italic(p) == "8.99e-07"),
           size = 4, colour = "black", fontface = "bold") +
  #annotate(geom = "text", x = 3.5, y = 1230, label = expression(Guild ~~ (G): H == "16.190;" ~~ df == "2;" ~~ italic(p) == "3.05e-04"),
  #         size = 4, colour = "black", fontface = "bold") +
  #annotate(geom = "text", x = 3.5, y = 1140, label = expression("P x G:" ~~ H == "29.394;" ~~ df == "2;" ~~ italic(p) == "4.14e-07"),
  #         size = 4, colour = "black", fontface = "bold") +
  scale_fill_d3(palette = "category20") +
  scale_colour_d3(palette = "category20") +
  scale_x_discrete(labels = c("Ectomycorrhizal\nfungi", "Saprotrophic\nfungi", "Plant pathogenic\nfungi",
                              "Ectomycorrhizal\nfungi", "Saprotrophic\nfungi", "Plant pathogenic\nfungi")) +
  scale_y_continuous(limits = c(-50, 625),
                     expand = c(0, 0)) +
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
    aspect.ratio = 0.5)
fig2a_rCNV_FG_AscBas_IQR

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2a_IQR_pdf <- str_c("fig2a_", "rCNV_FG_AscBas_IQR_", tm, ".pdf", sep = "")
fig2a_IQR_jpg <- str_c("fig2a_", "rCNV_FG_AscBas_IQR_", tm, ".jpg", sep = "")

ggsave(fig2a_IQR_pdf, fig2a_rCNV_FG_AscBas_IQR, width = 10.9, height = 6.17)
ggsave(fig2a_IQR_jpg, fig2a_rCNV_FG_AscBas_IQR, width = 10.9, height = 6.17)


# fig2b, AMF vs ECM -------------------------
rCNV_taxa_FGanno_subgrp$Guild1 %>% table()

rCNV_taxa_FGanno_subgrp_AMFECM <- 
  rCNV_taxa_FGanno_subgrp %>% filter(Guild1 == "Arbuscular Mycorrhizal" | Guild1 == "Ectomycorrhizal") %>% 
  mutate(Guild2 = str_c(Guild1, "Fungi", sep = " "))

rCNV_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% range()
rCNV_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% pull() %>% mean()
rCNV_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% pull() %>% sd()

kruskal.test(Both_ITS_LSU ~ Guild2, data = rCNV_taxa_FGanno_subgrp_AMFECM)
# Kruskal-Wallis chi-squared = 91.32, df = 1, p-value < 2.2e-16

kruskal.test(Both_ITS_LSU ~ Guild2, data = rCNV_taxa_FGanno_subgrp_AMFECM) %>% tidy()

rCNV_FG_AMFECM_plot_subdata1 <- rCNV_taxa_FGanno_subgrp_AMFECM %>% 
  group_by(Guild2) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                 sd_rDNA = sd(Both_ITS_LSU),
                                 max_rDNA = max(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))


rCNV_FG_AMFECM_plot_subdata2 <- rCNV_taxa_FGanno_subgrp_AMFECM %>% 
  group_by(Guild2) %>% count(Guild2) %>% mutate(lab = str_c("n =", n, sep = " ")) %>% 
  left_join(rCNV_FG_AMFECM_plot_subdata1, by = "Guild2") %>% 
  mutate(lab = str_c(lab.x, lab.y, sep = " "))


rCNV_taxa_FGanno_subgrp_AMFECM %>% filter(Both_ITS_LSU > 300)

fig2b_rCNV_FG_AMFECM <- 
  ggplot(rCNV_taxa_FGanno_subgrp_AMFECM, aes(Guild2, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  #geom_text(data = rCNV_FG_AMFECM_plot_subdata1,
  #          aes(x = Guild2, y = max_rDNA, label = lab), vjust = -1, colour = "black") +
  geom_text(data = rCNV_FG_AMFECM_plot_subdata2,
            aes(x = Guild2, y = 0, label = lab), vjust = 2, colour = "black") + 
  annotate(geom = "text", x = 1.5, y = 270, label = expression("chi-square" == "91.320;" ~~ df == "1;" ~~ italic(p) == "1.22e-21"), size = 5) +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Arbuscular mycorrhizal\nfungi", "Ectomycorrhizal\nfungi"),
                   expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 270),
                     expand = c(0.1, 0)) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 1
    )
fig2b_rCNV_FG_AMFECM

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2b_pdf <- str_c("fig2b_", "rCNV_FG_AMFECM_", tm, ".pdf", sep = "")
fig2b_jpg <- str_c("fig2b_", "rCNV_FG_AMFECM_", tm, ".jpg", sep = "")

ggsave(fig2b_pdf, fig2b_rCNV_FG_AMFECM, width = 4.95, height = 4.74)
ggsave(fig2b_jpg, fig2b_rCNV_FG_AMFECM, width = 4.95, height = 4.74)


# IQR -------------
rCNV_taxa_FGanno_subgrp$Guild1 %>% table()


rCNV_taxa_FGanno_subgrp_AMFECM_IQR <- 
  rCNV_taxa_FGanno_subgrp_AMFECM %>% group_by(Guild2) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()

rCNV_FG_AMFECM_IQR_plot_subdata <- kruskal.test(Both_ITS_LSU ~ Guild2, data = rCNV_taxa_FGanno_subgrp_AMFECM_IQR) %>% tidy()

rCNV_FG_AMFECM_IQR_plot_subdata1 <- rCNV_taxa_FGanno_subgrp_AMFECM_IQR %>% 
  group_by(Guild2) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                 sd_rDNA = sd(Both_ITS_LSU),
                                 max_rDNA = max(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))


rCNV_FG_AMFECM_IQR_plot_subdata2 <- rCNV_taxa_FGanno_subgrp_AMFECM_IQR %>% 
  group_by(Guild2) %>% count(Guild2) %>% mutate(lab = str_c("n =", n, sep = " ")) %>% 
  left_join(rCNV_FG_AMFECM_IQR_plot_subdata1, by = "Guild2") %>% 
  mutate(lab = str_c(lab.x, lab.y, sep = " "))


fig2b_rCNV_FG_AMFECM_IQR <- 
  ggplot(rCNV_taxa_FGanno_subgrp_AMFECM_IQR, aes(Guild2, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  #geom_text(data = rCNV_FG_AMFECM_plot_subdata1,
  #          aes(x = Guild2, y = max_rDNA, label = lab), vjust = -1, colour = "black") +
  geom_text(data = rCNV_FG_AMFECM_IQR_plot_subdata2,
            aes(x = Guild2, y = 0, label = lab), vjust = 2, colour = "black") + 
  annotate(geom = "text", x = 1.5, y = 270, label = expression("chi-square" == "85.055;" ~~ df == "1;" ~~ italic(p) == "2.90e-20"), size = 5) +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Arbuscular mycorrhizal\nfungi", "Ectomycorrhizal\nfungi")) +
  scale_y_continuous(limits = c(0, 270),
                     expand = c(0.1, 0)) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 1)
fig2b_rCNV_FG_AMFECM_IQR

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2b_IQR_pdf <- str_c("fig2b_", "rCNV_FG_AMFECM_IQR_", tm, ".pdf", sep = "")
fig2b_IQR_jpg <- str_c("fig2b_", "rCNV_FG_AMFECM_IQR_", tm, ".jpg", sep = "")

ggsave(fig2b_IQR_pdf, fig2b_rCNV_FG_AMFECM_IQR, width = 5.65, height = 5.33)
ggsave(fig2b_IQR_jpg, fig2b_rCNV_FG_AMFECM_IQR, width = 5.65, height = 5.33)

# FigS3, Saprotroph ---------------
rCNV_taxa_FGanno_subgrp_sap <- 
  rCNV_taxa_FGanno_subgrp %>% filter(Guild2 == "Saprotroph Fungi") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota")
rCNV_taxa_FGanno_subgrp_sap

kru_sap_rlt <- kruskal(rCNV_taxa_FGanno_subgrp_sap$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp_sap$Guild1, p.adj = "fdr")
kru_sap_rlt

rCNV_FG_sap_maxrDNA <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild1) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_sap_plot_subdata1 <- 
  kru_sap_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild1 = rownames(.), .before = groups) %>%
  select(Guild1, groups) %>% left_join(rCNV_FG_sap_maxrDNA, by = "Guild1")

rCNV_FG_sap_plot_subdata2 <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild1) %>% count(Guild1)


rCNV_taxa_FGanno_subgrp_sap$Guild1 <- factor(rCNV_taxa_FGanno_subgrp_sap$Guild1,
                                             levels = c("Dung Saprotroph", "Plant Saprotroph", "Pollen Saprotroph",
                                                        "Wood Saprotroph", "Undefined Saprotroph"))

rCNV_FG_sap_plot <- 
  ggplot(rCNV_taxa_FGanno_subgrp_sap, aes(Guild1, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild1), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild1),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_sap_plot_subdata1,
            aes(x = Guild1, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "red",
            size = 5) +
  geom_text(data = rCNV_FG_sap_plot_subdata2,
            aes(x = Guild1, y = 0, label = n), vjust = 2, colour = "red") +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Plant\nsaprotroph\nfungi", "Wood\nsaprotroph\nfungi",
                              "Dung\nsaprotroph\nfungi", "Undefined\nsaprotroph\nfungi")) +
  scale_y_continuous(limits = c(0, 750),
                     expand = c(0.1, 0.1)) +
  labs(y = "rDNA copy number", x = "Saprotroph type") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm")
  )
rCNV_FG_sap_plot


# FigS3, Saprotroph AscBas split -------------

rCNV_taxa_FGanno_subgrp_sap <- 
  rCNV_taxa_FGanno_subgrp_sap %>% mutate(Guild3 = str_c(phy, Guild1, sep = "-"))

kru_sap_AscBas_rlt <- kruskal(rCNV_taxa_FGanno_subgrp_sap$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp_sap$Guild3, p.adj = "fdr")
kru_sap_AscBas_rlt

rCNV_FG_sap_AscBas_maxrDNA <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_sap_AscBas_plot_subdata1 <- 
  kru_sap_AscBas_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups) %>% left_join(rCNV_FG_sap_AscBas_maxrDNA, by = "Guild3")

rCNV_FG_sap_AscBas_plot_subdata1

rCNV_taxa_FGanno_subgrp_sap_MeSd <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                 sd_rDNA = sd(Both_ITS_LSU))

rCNV_FG_sap_AscBas_plot_subdata1 <- rCNV_FG_sap_AscBas_plot_subdata1 %>% left_join(rCNV_taxa_FGanno_subgrp_sap_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

rCNV_FG_sap_AscBas_plot_subdata2 <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% count(Guild3) %>% mutate(
    lab = str_c("n =", n, sep = " ")
  )
rCNV_FG_sap_AscBas_plot_subdata2

rCNV_taxa_FGanno_subgrp_sap_tail <- 
  rCNV_taxa_FGanno_subgrp_sap %>% filter(Guild3 == "Basidiomycota-Dung Saprotroph")

rCNV_taxa_FGanno_subgrp_sap$Guild3 <- factor(rCNV_taxa_FGanno_subgrp_sap$Guild3,
                                             levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                        "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                        "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                        "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


rCNV_taxa_FGanno_subgrp_sap1 <- rbind(rCNV_taxa_FGanno_subgrp_sap, rCNV_taxa_FGanno_subgrp_sap_tail)

rCNV_taxa_FGanno_subgrp_sap1$Guild3 <- factor(rCNV_taxa_FGanno_subgrp_sap1$Guild3,
                                              levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                         "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                         "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                         "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


rCNV_taxa_FGanno_subgrp_sap_SRH <- rCNV_taxa_FGanno_subgrp_sap %>% filter(Guild3 != "Basidiomycota-Dung Saprotroph")
rCNV_taxa_FGanno_subgrp_sap_SRH

scheirerRayHare(Both_ITS_LSU ~ phy * Guild1, data = rCNV_taxa_FGanno_subgrp_sap_SRH)

rCNV_FG_sap_AscBas_plot_subdata1 <- rCNV_FG_sap_AscBas_plot_subdata1 %>% 
  mutate(groups = if_else(groups == "bc", "---", groups))

#fix(rCNV_FG_sap_AscBas_plot_subdata1)

figS3_rCNV_FG_sap_AscBas <- 
  ggplot(rCNV_taxa_FGanno_subgrp_sap, aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(data = rCNV_taxa_FGanno_subgrp_sap1,
                   side = "r", colour = NA, alpha = 0.8, aes(fill = Guild1)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild1),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_sap_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black") +
  geom_text(data = rCNV_FG_sap_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 30, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = rCNV_FG_sap_AscBas_plot_subdata2,
            aes(x = Guild3, y = 0, label = lab), vjust = 2, colour = "black") +
  geom_segment(x = 4.5, xend = 4.5, y = -20, yend = 510, linetype = 2) +
  annotate(geom = "text", x = 2.5, y = 500, label = "Ascomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 6.5, y = 500, label = "Basidiomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 4.5, y = 600, label = expression(Phylum ~~ (P): H == "120.283;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~
                                                                 Guild ~~ (G): H == "2.787;" ~~ df == "3;" ~~ italic(p) == "0.425;" ~~ 
                                                                 "P x G:" ~~ H == "7.330;" ~~ df == "2;" ~~ italic(p) == "0.0256"),
           size = 3.5, colour = "black", fontface = "bold") +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Plant\nsaprotrophic\nfungi", "Wood\nsaprotrophic\nfungi", "Dung\nsaprotrophic\nfungi", "Undefined\nsaprotrophic\nfungi",
                              "Plant\nsaprotrophic\nfungi", "Wood\nsaprotrophic\nfungi", "Dung\nsaprotrophic\nfungi", "Undefined\nsaprotrophic\nfungi")) +
  scale_y_continuous(limits = c(-50, 620)) +
  labs(y = "rDNA copy number", x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 0.5
  )
figS3_rCNV_FG_sap_AscBas

tm <- now() %>% str_split_i(pattern = " ", 1)
figS3_pdf <- str_c("figS3_", "rCNV_FG_sap_AscBas_", tm, ".pdf", sep = "")
figS3_jpg <- str_c("figS3_", "rCNV_FG_sap_AscBas_", tm, ".jpg", sep = "")

ggsave(figS3_pdf, figS3_rCNV_FG_sap_AscBas, width = 10.6, height = 6.02)
ggsave(figS3_jpg, figS3_rCNV_FG_sap_AscBas, width = 10.6, height = 6.02)


# IQR sap -----------------
rCNV_taxa_FGanno_subgrp_sap_IQR <- 
  rCNV_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()

kru_sap_AscBas_IQR_rlt <- kruskal(rCNV_taxa_FGanno_subgrp_sap_IQR$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp_sap_IQR$Guild3, p.adj = "fdr")
kru_sap_AscBas_IQR_rlt

rCNV_FG_sap_AscBas_IQR_maxrDNA <- 
  rCNV_taxa_FGanno_subgrp_sap_IQR %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_sap_AscBas_IQR_plot_subdata1 <- 
  kru_sap_AscBas_IQR_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups) %>% left_join(rCNV_FG_sap_AscBas_IQR_maxrDNA, by = "Guild3")

rCNV_FG_sap_AscBas_IQR_plot_subdata1

rCNV_taxa_FGanno_subgrp_sap_IQR_MeSd <- 
  rCNV_taxa_FGanno_subgrp_sap_IQR %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                     sd_rDNA = sd(Both_ITS_LSU))

rCNV_FG_sap_AscBas_IQR_plot_subdata1 <- rCNV_FG_sap_AscBas_IQR_plot_subdata1 %>% left_join(rCNV_taxa_FGanno_subgrp_sap_IQR_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

rCNV_FG_sap_AscBas_IQR_plot_subdata2 <- 
  rCNV_taxa_FGanno_subgrp_sap_IQR %>% group_by(Guild3) %>% count(Guild3) %>% mutate(
    lab = str_c("n =", n, sep = " ")
  )
rCNV_FG_sap_AscBas_IQR_plot_subdata2

rCNV_taxa_FGanno_subgrp_sap_IQR_tail <- 
  rCNV_taxa_FGanno_subgrp_sap_IQR %>% filter(Guild3 == "Basidiomycota-Dung Saprotroph")

rCNV_taxa_FGanno_subgrp_sap_IQR$Guild3 <- factor(rCNV_taxa_FGanno_subgrp_sap_IQR$Guild3,
                                                 levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                            "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                            "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                            "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


rCNV_taxa_FGanno_subgrp_IQR_sap1 <- rbind(rCNV_taxa_FGanno_subgrp_sap_IQR, rCNV_taxa_FGanno_subgrp_sap_IQR_tail)

rCNV_taxa_FGanno_subgrp_IQR_sap1$Guild3 <- factor(rCNV_taxa_FGanno_subgrp_IQR_sap1$Guild3,
                                              levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                         "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                         "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                         "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


rCNV_taxa_FGanno_subgrp_IQR_sap_SRH <- rCNV_taxa_FGanno_subgrp_sap_IQR %>% filter(Guild3 != "Basidiomycota-Dung Saprotroph")
rCNV_taxa_FGanno_subgrp_IQR_sap_SRH

scheirerRayHare(Both_ITS_LSU ~ phy * Guild1, data = rCNV_taxa_FGanno_subgrp_IQR_sap_SRH)

rCNV_FG_sap_AscBas_IQR_plot_subdata1 <- rCNV_FG_sap_AscBas_IQR_plot_subdata1 %>% mutate(groups = if_else(groups == "bcd", "---", groups))

figS3_rCNV_FG_sap_AscBas_IQR <- 
  ggplot(rCNV_taxa_FGanno_subgrp_sap_IQR, aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(data = rCNV_taxa_FGanno_subgrp_IQR_sap1,
                   side = "r", colour = NA, alpha = 0.8, aes(fill = Guild1)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild1),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_sap_AscBas_IQR_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black") +
  geom_text(data = rCNV_FG_sap_AscBas_IQR_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 20, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = rCNV_FG_sap_AscBas_IQR_plot_subdata2,
            aes(x = Guild3, y = 0, label = lab), vjust = 2, colour = "black") +
  geom_segment(x = 4.5, xend = 4.5, y = -20, yend = 370, linetype = 2) +
  annotate(geom = "text", x = 2.5, y = 350, label = "Ascomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 6.5, y = 350, label = "Basidiomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 4.5, y = 400, label = expression(Phylum ~~ (P): H == "117.674;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~
                                                                 Guild ~~ (G): H == "3.495;" ~~ df == "3;" ~~ italic(p) == "0.321;" ~~ 
                                                                 "P x G:" ~~ H == "6.525;" ~~ df == "2;" ~~ italic(p) == "0.0383"),
           size = 4, colour = "black", fontface = "bold") +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Plant\nsaprotrophic\nfungi", "Wood\nsaprotrophic\nfungi", "Dung\nsaprotrophic\nfungi", "Undefined\nsaprotrophic\nfungi",
                              "Plant\nsaprotrophic\nfungi", "Wood\nsaprotrophic\nfungi", "Dung\nsaprotrophic\nfungi", "Undefined\nsaprotrophic\nfungi")) +
  scale_y_continuous(limits = c(0, 400),
                     expand = c(0.1, 0.1)) +
  labs(y = "rDNA copy number", x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 0.5
  )
figS3_rCNV_FG_sap_AscBas_IQR

tm <- now() %>% str_split_i(pattern = " ", 1)
figS3_IQR_pdf <- str_c("figS3_", "rCNV_FG_sap_AscBas_IQR_", tm, ".pdf", sep = "")
figS3_IQR_jpg <- str_c("figS3_", "rCNV_FG_sap_AscBas_IQR_", tm, ".jpg", sep = "")

ggsave(figS3_IQR_pdf, figS3_rCNV_FG_sap_AscBas_IQR, width = 10.6, height = 6.02)
ggsave(figS3_IQR_jpg, figS3_rCNV_FG_sap_AscBas_IQR, width = 10.6, height = 6.02)

# PP split---------------
rCNV_taxa_FGanno_subgrp_pp <- 
  rCNV_taxa_FGanno_subgrp %>% filter(Guild2 == "Plant Pathogen") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota")
rCNV_taxa_FGanno_subgrp_pp


# Asc Bas split -------------
rCNV_taxa_FGanno_subgrp_pp <- 
  rCNV_taxa_FGanno_subgrp_pp %>% mutate(Guild3 = str_c(phy, Guild1, sep = "-"))


kru_pp_AscBas_rlt <- kruskal(rCNV_taxa_FGanno_subgrp_pp$Both_ITS_LSU, rCNV_taxa_FGanno_subgrp_pp$Guild3, p.adj = "fdr")
kru_pp_AscBas_rlt

rCNV_FG_pp_AscBas_maxrDNA <- 
  rCNV_taxa_FGanno_subgrp_pp %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_FG_pp_AscBas_plot_subdata1 <- 
  kru_pp_AscBas_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups) %>% left_join(rCNV_FG_pp_AscBas_maxrDNA, by = "Guild3")

rCNV_FG_pp_AscBas_plot_subdata1

rCNV_FG_pp_AscBas_plot_subdata2 <- 
  rCNV_taxa_FGanno_subgrp_pp %>% group_by(Guild3) %>% count(Guild3)
rCNV_FG_pp_AscBas_plot_subdata2


kruskal.test(Both_ITS_LSU ~ phy, data = rCNV_taxa_FGanno_subgrp_pp)

rCNV_FG_pp_AscBas <- 
  ggplot(rCNV_taxa_FGanno_subgrp_pp, aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(data = rCNV_taxa_FGanno_subgrp_pp,
                   side = "r", colour = NA, alpha = 0.8, aes(fill = phy)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = phy),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_FG_pp_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "red",
            size = 5) +
  geom_text(data = rCNV_FG_pp_AscBas_plot_subdata2,
            aes(x = Guild3, y = 0, label = n), vjust = 2, colour = "red") +
  annotate(geom = "text", x = 1.5, y = 630, label = expression("chi-square" == "0.681;" ~~ df == "1;" ~~ italic(p) == "0.409"), size = 5) +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Ascomycota plant pathogen\nfungi", "Basidiomycota plant pathogen\nfungi")) +
  scale_y_continuous(limits = c(0, 680),
                     expand = c(0.1, 0.1)) +
  labs(y = "rDNA copy number", x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black",
                             hjust = 0.5, vjust = 0.5),
    axis.title.x = element_text(face = "bold",
                                size = 15,
                                vjust = -1.5),
    axis.title.y = element_text(face = "bold",
                                size = 15,
                                vjust = 2),
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    aspect.ratio = 1
  )
rCNV_FG_pp_AscBas



# Fungal traits database ------------------------
Fungaltrait_database <- read_excel("./database/Fungaltrait_table.xlsx", sheet = 1)
FT_Ee <- Fungaltrait_database %>% select(GENUS, Ectomycorrhiza_exploration_type_template...9) %>% drop_na()

rCNV_taxa_FGanno_subgrp_AMFECM
rCNV_taxa_FGanno_subgrp_ECM <- rCNV_taxa_FGanno_subgrp_AMFECM %>% 
  filter(Guild1 == "Ectomycorrhizal")

colnames(FT_Ee)[1] <- "gen"
colnames(FT_Ee)[2] <- "Ectomycorrhiza_exploration_type_template"

rCNV_ECM_FG_FT_Ee <- 
  rCNV_taxa_FGanno_subgrp_ECM %>% left_join(FT_Ee, by = "gen") %>% drop_na()

rCNV_ECM_FG_FT_Ee$Ectomycorrhiza_exploration_type_template %>% table()

rCNV_ECM_FG_FT_Ee_filterd <- 
  rCNV_ECM_FG_FT_Ee %>% 
  filter(Ectomycorrhiza_exploration_type_template != "mat") %>% 
  filter(Ectomycorrhiza_exploration_type_template != "unknown")

rCNV_ECM_FG_FT_Ee_filterd <- 
  rCNV_ECM_FG_FT_Ee_filterd %>% 
  mutate(Ect_exp_type = if_else(str_detect(Ectomycorrhiza_exploration_type_template, "short-distance"), "short-distance", 
                                Ectomycorrhiza_exploration_type_template))

rCNV_ECM_FG_FT_Ee_filterd$Ect_exp_type %>% table()

rCNV_ECM_FG_FT_Ee_filterd$Ect_exp_type <- 
  factor(rCNV_ECM_FG_FT_Ee_filterd$Ect_exp_type,
         levels = c("contact",
                    "short-distance",
                    "medium-distance_fringe","medium-distance_smooth",
                    "long-distance"))

# fig2c, ECM forage type -----------
kru_ECM <- kruskal(rCNV_ECM_FG_FT_Ee_filterd$Both_ITS_LSU, rCNV_ECM_FG_FT_Ee_filterd$Ect_exp_type, p.adj = "fdr")
kru_ECM$groups

rCNV_ECM_FG_FT_Ee_filterd

rCNV_ECM_maxrDNA <- 
  rCNV_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% summarise(max_rDNA = max(Both_ITS_LSU))


rCNV_ECM_FG_FT_Ee_filterd_subdata1 <- 
  kru_ECM$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type = rownames(.), .before = groups) %>% 
  select(Ect_exp_type, groups) %>% left_join(rCNV_ECM_maxrDNA, by = "Ect_exp_type")

rCNV_ECM_FG_FT_Ee_filterd_MeSd <- rCNV_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% 
  summarise(m_rDNA = mean(Both_ITS_LSU),
            sd_rDNA = sd(Both_ITS_LSU))

rCNV_ECM_FG_FT_Ee_filterd_subdata1 <- 
  rCNV_ECM_FG_FT_Ee_filterd_subdata1 %>% left_join(rCNV_ECM_FG_FT_Ee_filterd_MeSd, by = "Ect_exp_type") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

rCNV_ECM_FG_FT_Ee_filterd_subdata2 <- 
  rCNV_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% count(Ect_exp_type) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

kruskal.test(Both_ITS_LSU ~ Ect_exp_type, data = rCNV_ECM_FG_FT_Ee_filterd)
# Kruskal-Wallis chi-squared = 31.114, df = 4, p-value = 2.902e-06

fig2c_rCNV_ECM_Ee <- 
  ggplot(rCNV_ECM_FG_FT_Ee_filterd, aes(Ect_exp_type, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 20, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_subdata2,
            aes(x = Ect_exp_type, y = 0, label = lab), vjust = 2, colour = "black") +
  annotate(geom = "text", x = 3, y = 270, label = expression("chi-square" == "31.114;" ~~ df == "4;" ~~ italic(p) == "2.902e-06"), size = 5) +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Contact",
                              "Short\ndistance",
                              "Medium\ndistance fringe", "Medium\ndistance smooth",
                              "Long\ndistance"),
                   expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20, 280)) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black"),
    axis.title.x = element_text(face = "bold",
                                size = 15),
    axis.title.y = element_text(face = "bold",
                                size = 15)
  )
fig2c_rCNV_ECM_Ee

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2c_pdf <- str_c("fig2c_", "rCNV_ECM_Ee_", tm, ".pdf", sep = "")
fig2c_jpg <- str_c("fig2c_", "rCNV_ECM_Ee_", tm, ".jpg", sep = "")

ggsave(fig2c_pdf, fig2c_rCNV_ECM_Ee, width = 7.9, height = 4.74)
ggsave(fig2c_jpg, fig2c_rCNV_ECM_Ee, width = 7.9, height = 4.74)

# fig1e, IQR ------------
rCNV_ECM_FG_FT_Ee_filterd_IQR <- 
  rCNV_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()

kru_ECM_IQR <- kruskal(rCNV_ECM_FG_FT_Ee_filterd_IQR$Both_ITS_LSU, rCNV_ECM_FG_FT_Ee_filterd_IQR$Ect_exp_type, p.adj = "fdr")
kru_ECM_IQR$groups

rCNV_ECM_FG_FT_Ee_filterd_IQR

rCNV_ECM_IQR_maxrDNA <- 
  rCNV_ECM_FG_FT_Ee_filterd_IQR %>% group_by(Ect_exp_type) %>% summarise(max_rDNA = max(Both_ITS_LSU))


rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata1 <- 
  kru_ECM_IQR$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type = rownames(.), .before = groups) %>% 
  select(Ect_exp_type, groups) %>% left_join(rCNV_ECM_IQR_maxrDNA, by = "Ect_exp_type")

rCNV_ECM_FG_FT_Ee_IQR_filterd_MeSd <- rCNV_ECM_FG_FT_Ee_filterd_IQR %>% group_by(Ect_exp_type) %>% 
  summarise(m_rDNA = mean(Both_ITS_LSU),
            sd_rDNA = sd(Both_ITS_LSU))

rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata1 <- 
  rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata1 %>% left_join(rCNV_ECM_FG_FT_Ee_IQR_filterd_MeSd, by = "Ect_exp_type") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata2 <- 
  rCNV_ECM_FG_FT_Ee_filterd_IQR %>% group_by(Ect_exp_type) %>% count(Ect_exp_type) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

kruskal.test(Both_ITS_LSU ~ Ect_exp_type, data = rCNV_ECM_FG_FT_Ee_filterd_IQR)

fig1e_rCNV_ECM_Ee_IQR <- 
  ggplot(rCNV_ECM_FG_FT_Ee_filterd_IQR, aes(Ect_exp_type, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 20, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_IQR_filterd_subdata2,
            aes(x = Ect_exp_type, y = 0, label = lab), vjust = 2, colour = "black") +
  annotate(geom = "text", x = 3, y = 280, label = expression("chi-square" == "43.406;" ~~ df == "4;" ~~ italic(p) == "8.524e-09"), size = 5.5) +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Contact",
                              "Short\ndistance",
                              "Medium\ndistance fringe", "Medium\ndistance smooth",
                              "Long\ndistance")) +
  scale_y_continuous(limits = c(0, 290),
                     expand = c(0.1, 0)) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black"),
    axis.title.x = element_text(face = "bold",
                                size = 15),
    axis.title.y = element_text(face = "bold",
                                size = 15)
  )
fig1e_rCNV_ECM_Ee_IQR

ggsave("fig1e_rCNV_ECM_Ee_IQR_20250210.pdf", fig1e_rCNV_ECM_Ee_IQR, width = 8.53, height = 4.46)
ggsave("fig1e_rCNV_ECM_Ee_IQR_20250210.jpg", fig1e_rCNV_ECM_Ee_IQR, width = 8.53, height = 4.46)

# figS4, forage type between Asc and Bas ---------------------
rCNV_ECM_FG_FT_Ee_filterd_phy <-
  rCNV_ECM_FG_FT_Ee_filterd %>% 
  filter(phy != "Mucoromycota") %>% 
  filter(Ect_exp_type == "contact" | Ect_exp_type == "short-distance") %>% 
  mutate(Ect_exp_type1 = str_c(phy, Ect_exp_type, sep = "-"))

kru_ECM_phy <- kruskal(rCNV_ECM_FG_FT_Ee_filterd_phy$Both_ITS_LSU,
                       rCNV_ECM_FG_FT_Ee_filterd_phy$Ect_exp_type1, p.adj = "fdr")
kru_ECM_phy

rCNV_ECM_phy_maxrDNA <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_ECM_FG_FT_Ee_filterd_phy_subdata1 <- 
  kru_ECM_phy$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type1 = rownames(.), .before = groups) %>% 
  select(Ect_exp_type1, groups) %>% left_join(rCNV_ECM_phy_maxrDNA, by = "Ect_exp_type1")

rCNV_ECM_FG_FT_Ee_filterd_phy_subdata2 <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% count(Ect_exp_type1)

rCNV_ECM_FG_FT_Ee_filterd_phy_MeSd <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                              sd_rDNA = sd(Both_ITS_LSU))
rCNV_ECM_FG_FT_Ee_filterd_phy_subdata2 <-
  rCNV_ECM_FG_FT_Ee_filterd_phy_subdata2 %>% left_join(rCNV_ECM_FG_FT_Ee_filterd_phy_MeSd, by = "Ect_exp_type1") %>% 
  mutate(lab = str_c("n = ", n, " ", "(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

rCNV_ECM_FG_FT_Ee_filterd_phy$phy %>% table()
rCNV_ECM_FG_FT_Ee_filterd_phy$Ect_exp_type1 %>% table()

scheirerRayHare(Both_ITS_LSU ~ phy + Ect_exp_type, data = rCNV_ECM_FG_FT_Ee_filterd_phy)

figS4_rCNV_ECM_Ee_AscBas <- 
  ggplot(rCNV_ECM_FG_FT_Ee_filterd_phy, aes(Ect_exp_type1, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_phy_subdata1,
            aes(x = Ect_exp_type1, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_phy_subdata2,
            aes(x = Ect_exp_type1, y = 0, label = lab), vjust = 2, hjust = 0.5, colour = "black") +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Contact", "Short distance",
                              "Contact", "Short distance")) +
  scale_y_continuous(limits = c(0, 300),
                     expand = c(0.1, 0.1)) +
  geom_segment(x = 2.5, xend = 2.5, y = -20, yend = 270, linetype = 2) +
  annotate(geom = "text", x = 1.5, y = 250, label = "Ascomycota", size = 6, fontface = "bold.italic") +
  annotate(geom = "text", x = 3.5, y = 250, label = "Basidiomycota", size = 6, fontface = "bold.italic") +
  annotate(geom = "text", x = 2.5, y = 295, label = expression(Phylum ~~ (P): H == "4.798;" ~~ df == "1;" ~~ italic(p) == "0.0285" ~~
                                                                  Guild ~~ (G): H == "10.954;" ~~ df == "1;" ~~ italic(p) == "0.000933" ~~
                                                                  "P x G:" ~~ H == "10.757;" ~~ df == "1;" ~~ italic(p) == "0.00104"),
           size = 4, colour = "black", fontface = "bold") +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black"),
    axis.title.x = element_text(face = "bold",
                                size = 15),
    axis.title.y = element_text(face = "bold",
                                size = 15),
    aspect.ratio = 0.5
  )
figS4_rCNV_ECM_Ee_AscBas

tm <- now() %>% str_split_i(pattern = " ", 1)
figS4_pdf <- str_c("figS4_", "rCNV_ECM_Ee_AscBas_", tm, ".pdf", sep = "")
figS4_jpg <- str_c("figS4_", "rCNV_ECM_Ee_AscBas_", tm, ".jpg", sep = "")

ggsave(figS4_pdf, figS4_rCNV_ECM_Ee_AscBas, width = 9.96, height = 5.12)
ggsave(figS4_jpg, figS4_rCNV_ECM_Ee_AscBas, width = 9.96, height = 5.12)

# figS4, IQR -----------
rCNV_ECM_FG_FT_Ee_filterd_phy_IQR <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()

kru_ECM_phy_IQR <- kruskal(rCNV_ECM_FG_FT_Ee_filterd_phy_IQR$Both_ITS_LSU,
                       rCNV_ECM_FG_FT_Ee_filterd_phy_IQR$Ect_exp_type1, p.adj = "fdr")
kru_ECM_phy_IQR

rCNV_ECM_phy_IQR_maxrDNA <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy_IQR %>% group_by(Ect_exp_type1) %>% summarise(max_rDNA = max(Both_ITS_LSU))

rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata1 <- 
  kru_ECM_phy_IQR$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type1 = rownames(.), .before = groups) %>% 
  select(Ect_exp_type1, groups) %>% left_join(rCNV_ECM_phy_IQR_maxrDNA, by = "Ect_exp_type1")

rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata2 <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy_IQR %>% group_by(Ect_exp_type1) %>% count(Ect_exp_type1)

rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_MeSd <- 
  rCNV_ECM_FG_FT_Ee_filterd_phy_IQR %>% group_by(Ect_exp_type1) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                              sd_rDNA = sd(Both_ITS_LSU))
rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata2 <-
  rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata2 %>% left_join(rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_MeSd, by = "Ect_exp_type1") %>% 
  mutate(lab = str_c("n = ", n, " ", "(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))


scheirerRayHare(Both_ITS_LSU ~ phy + Ect_exp_type, data = rCNV_ECM_FG_FT_Ee_filterd_phy_IQR)

figS4_rCNV_ECM_Ee_AscBas_IQR <- 
  ggplot(rCNV_ECM_FG_FT_Ee_filterd_phy_IQR, aes(Ect_exp_type1, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata1,
            aes(x = Ect_exp_type1, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = rCNV_ECM_FG_FT_Ee_filterd_phy_IQR_subdata2,
            aes(x = Ect_exp_type1, y = 0, label = lab), vjust = 2, hjust = 0.5, colour = "black") +
  scale_colour_d3(palette = "category20") +
  scale_fill_d3(palette = "category20") +
  scale_x_discrete(labels = c("Contact", "Short distance",
                              "Contact", "Short distance")) +
  scale_y_continuous(limits = c(0, 300),
                     expand = c(0.1, 0.1)) +
  geom_segment(x = 2.5, xend = 2.5, y = -20, yend = 270, linetype = 2) +
  annotate(geom = "text", x = 1.5, y = 250, label = "Ascomycota", size = 6, fontface = "bold.italic") +
  annotate(geom = "text", x = 3.5, y = 250, label = "Basidiomycota", size = 6, fontface = "bold.italic") +
  annotate(geom = "text", x = 2.5, y = 295, label = expression(Phylum ~~ (P): H == "4.326;" ~~ df == "1;" ~~ italic(p) == "0.0375" ~~
                                                                 Guild ~~ (G): H == "11.397;" ~~ df == "1;" ~~ italic(p) == "0.000" ~~
                                                                 "P x G:" ~~ H == "10.039;" ~~ df == "1;" ~~ italic(p) == "0.00153"),
           size = 4, colour = "black", fontface = "bold") +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold",
                             size = 12,
                             colour = "black"),
    axis.title.x = element_text(face = "bold",
                                size = 15),
    axis.title.y = element_text(face = "bold",
                                size = 15),
    aspect.ratio = 0.5
  )
figS4_rCNV_ECM_Ee_AscBas_IQR

tm <- now() %>% str_split_i(pattern = " ", 1)
figS4_IQR_pdf <- str_c("figS4_", "rCNV_ECM_Ee_AscBas_IQR_", tm, ".pdf", sep = "")
figS4_IQR_jpg <- str_c("figS4_", "rCNV_ECM_Ee_AscBas_IQR_", tm, ".jpg", sep = "")

ggsave(figS4_IQR_pdf, figS4_rCNV_ECM_Ee_AscBas_IQR, width = 9.96, height = 5.12)
ggsave(figS4_IQR_pdf, figS4_rCNV_ECM_Ee_AscBas_IQR, width = 9.96, height = 5.12)

# done...

