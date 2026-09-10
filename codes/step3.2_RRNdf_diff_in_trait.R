
# step3: differences among fungal traits in FunGuild Database.
# fig2a, fig2b, fig2c, figS3, figS4
# Qiushi-Li, IM-CAS
# 2024.11.30


# packages
library(tidyverse)
library(readxl)

# install.packages("rcompanion")
library(rcompanion) # SHR test
library(agricolae) # kruskal.test
library(PMCMRplus)    # kwAllPairsDunnTest

# devtools::install_github('erocoar/gghalves')
library(gghalves)
library(ggsci)
library(patchwork)
library(broom)

library(writexl)

# load RData ----------
#load("rCNV_version_1.0_20250219.RData")

# FunGuild database --------------------------
FunGuild_database <- read_csv("./2.database/FunGuild_v2024.csv")
colnames(FunGuild_database)

Tax_Lev <- c(13, 9, 7, 3)
Tax_Ind <- c("gen", "fam", "ord", "phy")
FG_index_df <- data.frame(
  Tax_Lev,
  Tax_Ind
)
FG_index_df

# FRRN total tmp table
FRRN_rlt_taxa_spl2
# write_xlsx(FRRN_rlt_taxa_spl2, "FRRN_taxa.xlsx")



FRRN_taxa_FG <- FRRN_rlt_taxa_spl2 %>%
  select(project, Both_ITS_LSU, Name, phy, cla, ord, fam, Gen)
colnames(FRRN_taxa_FG)[8] <- "gen"

FRRN_taxa_FG <- FRRN_taxa_FG %>% mutate(
  cla = str_sub(cla, start = 3),
  ord = str_sub(ord, start = 3),
  fam = str_sub(fam, start = 3),
  gen = str_sub(gen, start = 3)
)


# Full align ----
FunGuild_database1 <- FunGuild_database %>% select(-Notes, -`Citation/Source`)
colnames(FunGuild_database1)

# view(FunGuild_database1)

FunGuild_database1$Trophic_Mode %>% table()


FG_full_tmp_gen <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[1]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_gen)[1] <- FG_index_df$Tax_Ind[1]

FRRN_taxa_FG_full_gen <- FRRN_taxa_FG %>% 
  left_join(FG_full_tmp_gen)

FRRN_taxa_FG_full_gen_na <- FRRN_taxa_FG_full_gen %>% filter(is.na(Trophic_Mode))
FRRN_taxa_FG_full_gen_na

FRRN_taxa_FG_full_gen_anno <- FRRN_taxa_FG_full_gen %>% filter(!is.na(Trophic_Mode))

FG_full_tmp_fam <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[2]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_fam)[1] <- FG_index_df$Tax_Ind[2]

FRRN_taxa_FG_full_fam <- FRRN_taxa_FG_full_gen_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_fam)

FRRN_taxa_FG_full_fam_na <- FRRN_taxa_FG_full_fam %>% filter(is.na(Trophic_Mode))
FRRN_taxa_FG_full_fam_na

FRRN_taxa_FG_full_fam_anno <- FRRN_taxa_FG_full_fam %>% filter(!is.na(Trophic_Mode))

FG_full_tmp_ord <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[3]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_ord)[1] <- FG_index_df$Tax_Ind[3]

FRRN_taxa_FG_full_ord <- FRRN_taxa_FG_full_fam_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_ord)

FRRN_taxa_FG_full_ord_na <- FRRN_taxa_FG_full_ord %>% filter(is.na(Trophic_Mode))
FRRN_taxa_FG_full_ord_na

FRRN_taxa_FG_full_ord_anno <- FRRN_taxa_FG_full_ord %>% filter(!is.na(Trophic_Mode))
FRRN_taxa_FG_full_ord_anno


FG_full_tmp_phy <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[4]) %>% select(-Taxon_Level)
colnames(FG_full_tmp_phy)[1] <- FG_index_df$Tax_Ind[4]

FRRN_taxa_FG_full_phy <- FRRN_taxa_FG_full_ord_na %>%
  select(-(Trophic_Mode:Confidence_Ranking)) %>% 
  left_join(FG_full_tmp_phy)

FRRN_taxa_FG_full_phy_na <- FRRN_taxa_FG_full_phy %>% filter(is.na(Trophic_Mode))
FRRN_taxa_FG_full_phy_na

FRRN_taxa_FG_full_phy_anno <- FRRN_taxa_FG_full_phy %>% filter(!is.na(Trophic_Mode))
FRRN_taxa_FG_full_phy_anno


FRRN_taxa_FGanno_full <- rbind(FRRN_taxa_FG_full_gen_anno, FRRN_taxa_FG_full_fam_anno, FRRN_taxa_FG_full_ord_anno, FRRN_taxa_FG_full_phy_anno, 
                               FRRN_taxa_FG_full_phy_na)
FRRN_taxa_FGanno_full

FRRN_taxa_FGanno_full %>% filter(is.na(Trophic_Mode))
# 10 unannotated


# may needed adjusted
FRRN_taxa_FGanno_full1 <- 
  FRRN_taxa_FGanno_full %>% 
  mutate(Guild1 = str_extract(FRRN_taxa_FGanno_full$Guild, "\\|([^\\|]+)\\|"), .before = Growth_Morphology) %>% 
  mutate(Guild1 = if_else(!is.na(Guild1), str_sub(Guild1, start = 2, end = -2), Guild))

# write_xlsx(rCNV_taxa_FGanno_full1, "rCNV_taxa_FG_new.xlsx")

FRRN_taxa_FGanno_full1$Guild1 %>% table()
FRRN_taxa_FGanno_full1 %>% filter(phy == "Glomeromycota") %>% 
  filter(Guild1 != "Arbuscular Mycorrhizal") %>% view()

FRRN_FG_subgrp_ind <- c("Plant Pathogen",
                        "Dung Saprotroph", "Pollen Saprotroph", "Wood Saprotroph", "Plant Saprotroph", "Undefined Saprotroph",
                        "Endophyte", "Epiphyte", "Arbuscular Mycorrhizal", "Ectomycorrhizal")

FRRN_taxa_FGanno_subgrp$phy %>% table()

FRRN_taxa_FGanno_subgrp <- 
  FRRN_taxa_FGanno_full1 %>% 
  filter(Guild1 %in% FRRN_FG_subgrp_ind) %>% 
  mutate(Guild2 = if_else(str_detect(Guild1, pattern = "Saprotroph"), "Saprotroph Fungi", Guild1), .before = Growth_Morphology)

# view(FRRN_taxa_FGanno_subgrp)


FRRN_taxa_FGanno_subgrp1 <- 
  FRRN_taxa_FGanno_subgrp %>% 
  filter(Guild2 == "Plant Pathogen" | Guild2 == "Saprotroph Fungi" | Guild2 == "Ectomycorrhizal") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota") %>% 
  mutate(Guild2 = if_else(str_detect(Guild2, pattern = "Plant Pathogen"), "Plant Pathogen Fungi", Guild2)) %>% 
  mutate(Guild2 = if_else(str_detect(Guild2, pattern = "Ectomycorrhizal"), "Ectomycorrhizal Fungi", Guild2))

FRRN_taxa_FGanno_subgrp1$Guild2 %>% table()



# fig2a, Asc Bas combine ------------------------
FRRN_taxa_FGanno_subgrp2 <- 
  FRRN_taxa_FGanno_subgrp1 %>% mutate(Guild3 = str_c(phy, Guild2, sep = "-"), .after = Guild2)

FRRN_taxa_FGanno_subgrp2$Guild3 <- factor(FRRN_taxa_FGanno_subgrp2$Guild3, levels = c(
  "Ascomycota-Ectomycorrhizal Fungi", "Ascomycota-Saprotroph Fungi", "Ascomycota-Plant Pathogen Fungi",
  "Basidiomycota-Ectomycorrhizal Fungi", "Basidiomycota-Saprotroph Fungi", "Basidiomycota-Plant Pathogen Fungi"
))


# kru_rlt_FG_subgrp2_subdata1


FRRN_FG_subgrp2_kru_tmp <- 
  kruskal(FRRN_taxa_FGanno_subgrp2$Both_ITS_LSU, FRRN_taxa_FGanno_subgrp2$Guild3, p.adj = "fdr")
FRRN_FG_subgrp2_kru_tmp


FRRN_FG_subgrp2_plot_subdata1 <- 
  FRRN_FG_subgrp2_kru$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups)
FRRN_FG_subgrp2_plot_subdata1

FRRN_FG_subgrp2_maxrDNA <- FRRN_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

FRRN_FG_subgrp2_plot_subdata1 <- FRRN_FG_subgrp2_plot_subdata1 %>% left_join(FRRN_FG_subgrp2_maxrDNA, by = "Guild3")

FRRN_FG_subgrp2_MeSd <- FRRN_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                                        sd_rDNA = sd(Both_ITS_LSU))
FRRN_FG_subgrp2_plot_subdata1 <- 
  FRRN_FG_subgrp2_plot_subdata1 %>% left_join(FRRN_FG_subgrp2_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

kru_rlt_FG_subgrp2_subdata2 <- FRRN_taxa_FGanno_subgrp2 %>% group_by(Guild3) %>% count(Guild3) %>% 
  mutate(lab = str_c("n =", n, sep = " "))


# kru pvalue --------
FRRN_taxa_FGanno_subgrp2$Guild3 %>% table()
FRRN_taxa_FGanno_subgrp2$Guild2 %>% table()

####### should be carefully checked #######
kruskal.test(Both_ITS_LSU ~ Guild3, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Ascomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 19.288, df = 1, p-value = 1.124e-05

kruskal.test(Both_ITS_LSU ~ Guild3, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Basidiomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 3.3805, df = 1, p-value = 0.06597
# Kruskal-Wallis chi-squared = 2.331, df = 1, p-value = 0.1268

kruskal.test(Both_ITS_LSU ~ phy, data = FRRN_taxa_FGanno_subgrp2 %>% filter(Guild2 == "Plant Pathogen Fungi"))


# SRH test
FRRN_AscBas_SRH <- scheirerRayHare(Both_ITS_LSU ~ phy * Guild2, data = FRRN_taxa_FGanno_subgrp2)
FRRN_AscBas_SRH
# phy:Guild2, p = 2.2398e-06

# poc test
# Ascomycota kruskal.test
kru_Asc <- kruskal.test(Both_ITS_LSU ~ Guild2, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Ascomycota"))
kru_Asc
# Kruskal-Wallis chi-squared = 19.401, df = 2, p-value = 6.126e-05

kruskal.test(Both_ITS_LSU ~ Guild3, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Ascomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 19.288, df = 1, p-value = 1.124e-05

# Basidiomycota kruskal.test
kru_Bas <- kruskal.test(Both_ITS_LSU ~ Guild2, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Basidiomycota"))
kru_Bas
# Kruskal-Wallis chi-squared = 31.061, df = 2, p-value = 1.799e-07

kruskal.test(Both_ITS_LSU ~ Guild3, data = FRRN_taxa_FGanno_subgrp2 %>% filter(phy == "Basidiomycota") %>% 
               filter(Guild2 != "Ectomycorrhizal Fungi"))
# Kruskal-Wallis chi-squared = 2.331, df = 1, p-value = 0.1268

# p.adj 
pvalue_AscBac <- c(kru_Asc$p.value, kru_Bas$p.value)
names(pvalue_AscBac) <- c("Asc", "Bas")
p_adj <- p.adjust(pvalue_AscBac, method = "fdr")
p_adj


# DunnTest
dunn_all <- kwAllPairsDunnTest(Both_ITS_LSU ~ Guild3, data = FRRN_taxa_FGanno_subgrp2,
                               p.adjust.method = "fdr")
summary(dunn_all)

# fix(FRRN_FG_subgrp2_plot_subdata1)
# 370


FRRN_taxa_FGanno_subgrp2 %>% filter(Both_ITS_LSU > 700)
FRRN_taxa_FGanno_subgrp$phy %>% table()


fig2a_FRRN_FG_AscBas <- 
  ggplot(FRRN_taxa_FGanno_subgrp2,
         aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_FG_subgrp2_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 25, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = FRRN_FG_subgrp2_plot_subdata1,
            aes(x = Guild3, y = max_rDNA, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = kru_rlt_FG_subgrp2_subdata2,
            aes(x = Guild3, y = -25, label = lab), colour = "black", hjust = 0.5) +
  geom_segment(x = 3.5, xend = 3.5, y = -20, yend = 580, linetype = 2) +
  annotate(geom = "text", x = 2, y = 510, label = "Ascomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 5, y = 510, label = "Basidiomycota", size = 7, fontface = "bold.italic") +
  annotate(geom = "text", x = 3.5, y = 630, label = expression(Phylum ~~ (P): H == "124.783;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~ 
                                                                 Guild ~~ (G): H == "22.366;" ~~ df == "2;" ~~ italic(p) == "1.391e-05;" ~~
                                                                 "P x G:" ~~ H == "26.018;" ~~ df == "2;" ~~ italic(p) == "2.240e-06"),
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
fig2a_FRRN_FG_AscBas
# > 700 un shown

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2a_pdf <- str_c("fig2a_", "FRRN_FG_AscBas_", tm, ".pdf", sep = "")
fig2a_jpg <- str_c("fig2a_", "FRRN_FG_AscBas_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_2a_pdf <- str_c(fig_path, fig2a_pdf)
fig_fullpath_2a_jpg <- str_c(fig_path, fig2a_jpg)

ggsave(fig_fullpath_2a_pdf, fig2a_FRRN_FG_AscBas, width = 10.6, height = 4.44)
ggsave(fig_fullpath_2a_jpg, fig2a_FRRN_FG_AscBas, width = 10.6, height = 4.44)


# fig2b, AMF vs ECM -------------------------
FRRN_taxa_FGanno_subgrp$Guild1 %>% table()
FRRN_taxa_FGanno_subgrp$phy %>% table()



FRRN_taxa_FGanno_subgrp$Guild2 %>% table()

FRRN_taxa_FGanno_subgrp %>% filter(Guild2 == "Ectomycorrhizal") %>% 
  filter(!project %in% FRRN_taxa_FGanno_subgrp1$project) %>% view()

FRRN_taxa_FGanno_subgrp1$project %>% table()

FRRN_taxa_FGanno_subgrp_AMFECM <- 
  FRRN_taxa_FGanno_subgrp %>% filter(Guild1 == "Arbuscular Mycorrhizal" | Guild1 == "Ectomycorrhizal") %>% 
  mutate(Guild2 = str_c(Guild1, "Fungi", sep = " "))

FRRN_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% range()
FRRN_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% pull() %>% mean()
FRRN_taxa_FGanno_subgrp_AMFECM %>% filter(Guild1 == "Arbuscular Mycorrhizal") %>% select(Both_ITS_LSU) %>% pull() %>% sd()

kruskal.test(Both_ITS_LSU ~ Guild2, data = FRRN_taxa_FGanno_subgrp_AMFECM)
# Kruskal-Wallis chi-squared = 91.823, df = 1, p-value < 2.2e-16

kruskal.test(Both_ITS_LSU ~ Guild2, data = FRRN_taxa_FGanno_subgrp_AMFECM) %>% tidy()

FRRN_FG_AMFECM_plot_subdata1 <- FRRN_taxa_FGanno_subgrp_AMFECM %>% 
  group_by(Guild2) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                 sd_rDNA = sd(Both_ITS_LSU),
                                 max_rDNA = max(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))


FRRN_FG_AMFECM_plot_subdata2 <- FRRN_taxa_FGanno_subgrp_AMFECM %>% 
  group_by(Guild2) %>% count(Guild2) %>% mutate(lab = str_c("n =", n, sep = " ")) %>% 
  left_join(FRRN_FG_AMFECM_plot_subdata1, by = "Guild2") %>% 
  mutate(lab = str_c(lab.x, lab.y, sep = " "))


FRRN_taxa_FGanno_subgrp_AMFECM %>% filter(Both_ITS_LSU > 300)

fig2b_FRRN_FG_AMFECM <- 
  ggplot(FRRN_taxa_FGanno_subgrp_AMFECM, aes(Guild2, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild2), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild2),
                        size = 1,
                        alpha = 0.8) +
  #geom_text(data = rCNV_FG_AMFECM_plot_subdata1,
  #          aes(x = Guild2, y = max_rDNA, label = lab), vjust = -1, colour = "black") +
  geom_text(data = FRRN_FG_AMFECM_plot_subdata2,
            aes(x = Guild2, y = 0, label = lab), vjust = 2, colour = "black") + 
  annotate(geom = "text", x = 1.5, y = 270, label = expression("chi-square" == "91.823;" ~~ df == "1;" ~~ italic(p) == "9.48e-22"), size = 5) +
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
fig2b_FRRN_FG_AMFECM

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2b_pdf <- str_c("fig2b_", "FRRN_FG_AMFECM_", tm, ".pdf", sep = "")
fig2b_jpg <- str_c("fig2b_", "FRRN_FG_AMFECM_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_2b_pdf <- str_c(fig_path, fig2b_pdf)
fig_fullpath_2b_jpg <- str_c(fig_path, fig2b_jpg)

ggsave(fig_fullpath_2b_pdf, fig2b_FRRN_FG_AMFECM, width = 4.95, height = 4.74)
ggsave(fig_fullpath_2b_jpg, fig2b_FRRN_FG_AMFECM, width = 4.95, height = 4.74)


# FigS3, Saprotroph ---------------
FRRN_taxa_FGanno_subgrp_sap <- 
  FRRN_taxa_FGanno_subgrp %>% filter(Guild2 == "Saprotroph Fungi") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota")
FRRN_taxa_FGanno_subgrp_sap

kru_sap_rlt <- kruskal(FRRN_taxa_FGanno_subgrp_sap$Both_ITS_LSU, FRRN_taxa_FGanno_subgrp_sap$Guild1, p.adj = "fdr")
kru_sap_rlt

FRRN_FG_sap_maxrDNA <- 
  FRRN_taxa_FGanno_subgrp_sap %>% group_by(Guild1) %>% summarise(max_rDNA = max(Both_ITS_LSU))

FRRN_FG_sap_plot_subdata1 <- 
  kru_sap_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild1 = rownames(.), .before = groups) %>%
  select(Guild1, groups) %>% left_join(FRRN_FG_sap_maxrDNA, by = "Guild1")

FRRN_FG_sap_plot_subdata2 <- 
  FRRN_taxa_FGanno_subgrp_sap %>% group_by(Guild1) %>% count(Guild1)


FRRN_taxa_FGanno_subgrp_sap$Guild1 <- factor(FRRN_taxa_FGanno_subgrp_sap$Guild1,
                                             levels = c("Dung Saprotroph", "Plant Saprotroph", "Pollen Saprotroph",
                                                        "Wood Saprotroph", "Undefined Saprotroph"))

FRRN_FG_sap_plot <- 
  ggplot(FRRN_taxa_FGanno_subgrp_sap, aes(Guild1, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, aes(fill = Guild1), alpha = 0.8) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild1),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_FG_sap_plot_subdata1,
            aes(x = Guild1, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "red",
            size = 5) +
  geom_text(data = FRRN_FG_sap_plot_subdata2,
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
FRRN_FG_sap_plot


# FigS3, Saprotroph AscBas split -------------

FRRN_taxa_FGanno_subgrp_sap <- 
  FRRN_taxa_FGanno_subgrp_sap %>% mutate(Guild3 = str_c(phy, Guild1, sep = "-"))

kru_sap_AscBas_rlt <- kruskal(FRRN_taxa_FGanno_subgrp_sap$Both_ITS_LSU,
                              FRRN_taxa_FGanno_subgrp_sap$Guild3, p.adj = "fdr")
kru_sap_AscBas_rlt

FRRN_FG_sap_AscBas_maxrDNA <- 
  FRRN_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

FRRN_FG_sap_AscBas_plot_subdata1 <- 
  kru_sap_AscBas_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups) %>% left_join(FRRN_FG_sap_AscBas_maxrDNA, by = "Guild3")

FRRN_FG_sap_AscBas_plot_subdata1

FRRN_taxa_FGanno_subgrp_sap_MeSd <- 
  FRRN_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                 sd_rDNA = sd(Both_ITS_LSU))

FRRN_FG_sap_AscBas_plot_subdata1 <- FRRN_FG_sap_AscBas_plot_subdata1 %>% 
  left_join(FRRN_taxa_FGanno_subgrp_sap_MeSd, by = "Guild3") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

FRRN_FG_sap_AscBas_plot_subdata2 <- 
  FRRN_taxa_FGanno_subgrp_sap %>% group_by(Guild3) %>% count(Guild3) %>% mutate(
    lab = str_c("n =", n, sep = " ")
  )
FRRN_FG_sap_AscBas_plot_subdata2

FRRN_taxa_FGanno_subgrp_sap_tail <- 
  FRRN_taxa_FGanno_subgrp_sap %>% filter(Guild3 == "Basidiomycota-Dung Saprotroph")

FRRN_taxa_FGanno_subgrp_sap$Guild3 <- factor(FRRN_taxa_FGanno_subgrp_sap$Guild3,
                                             levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                        "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                        "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                        "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


FRRN_taxa_FGanno_subgrp_sap1 <- rbind(FRRN_taxa_FGanno_subgrp_sap, FRRN_taxa_FGanno_subgrp_sap_tail)

FRRN_taxa_FGanno_subgrp_sap1$Guild3 <- factor(FRRN_taxa_FGanno_subgrp_sap1$Guild3,
                                              levels = c("Ascomycota-Plant Saprotroph", "Ascomycota-Wood Saprotroph",
                                                         "Ascomycota-Dung Saprotroph", "Ascomycota-Undefined Saprotroph",
                                                         "Basidiomycota-Plant Saprotroph", "Basidiomycota-Wood Saprotroph",
                                                         "Basidiomycota-Dung Saprotroph", "Basidiomycota-Undefined Saprotroph"))


FRRN_taxa_FGanno_subgrp_sap_SRH <- 
  FRRN_taxa_FGanno_subgrp_sap %>% filter(Guild3 != "Basidiomycota-Dung Saprotroph")
FRRN_taxa_FGanno_subgrp_sap_SRH

scheirerRayHare(Both_ITS_LSU ~ phy * Guild1, data = FRRN_taxa_FGanno_subgrp_sap_SRH)

FRRN_FG_sap_AscBas_plot_subdata1 <- FRRN_FG_sap_AscBas_plot_subdata1 %>% 
  mutate(groups = if_else(groups == "bc", "---", groups))

# fix(FRRN_FG_sap_AscBas_plot_subdata1)

figS3_FRRN_FG_sap_AscBas <- 
  ggplot(FRRN_taxa_FGanno_subgrp_sap, aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(data = FRRN_taxa_FGanno_subgrp_sap1,
                   side = "r", colour = NA, alpha = 0.8, aes(fill = Guild1)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Guild1),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_FG_sap_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black") +
  geom_text(data = FRRN_FG_sap_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA + 30, label = groups, vjust = -0.5),
            colour = "blue",
            size = 5) +
  geom_text(data = FRRN_FG_sap_AscBas_plot_subdata2,
            aes(x = Guild3, y = 0, label = lab), vjust = 2, colour = "black") +
  geom_segment(x = 4.5, xend = 4.5, y = -20, yend = 510, linetype = 2) +
  annotate(geom = "text", x = 2.5, y = 500, label = "Ascomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 6.5, y = 500, label = "Basidiomycota", size = 8, fontface = "bold.italic") +
  annotate(geom = "text", x = 4.5, y = 600, label = expression(Phylum ~~ (P): H == "118.163;" ~~ df == "1;" ~~ italic(p) == "0.000;" ~~
                                                                 Saprotrophic ~~ group ~~ (S): H == "3.205;" ~~ df == "3;" ~~ italic(p) == "0.361;" ~~ 
                                                                 "P x G:" ~~ H == "7.901;" ~~ df == "2;" ~~ italic(p) == "0.0192"),
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
figS3_FRRN_FG_sap_AscBas


tm <- now() %>% str_split_i(pattern = " ", 1)
figS3_pdf <- str_c("figS3_", "FRRN_FG_sap_AscBas_", tm, ".pdf", sep = "")
figS3_jpg <- str_c("figS3_", "FRRN_FG_sap_AscBas_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S3_pdf <- str_c(fig_path, figS3_pdf)
fig_fullpath_S3_jpg <- str_c(fig_path, figS3_jpg)


ggsave(fig_fullpath_S3_pdf, figS3_FRRN_FG_sap_AscBas, width = 10.6, height = 6.02)
ggsave(fig_fullpath_S3_jpg, figS3_FRRN_FG_sap_AscBas, width = 10.6, height = 6.02)



# PP split---------------
FRRN_taxa_FGanno_subgrp_pp <- 
  FRRN_taxa_FGanno_subgrp %>% filter(Guild2 == "Plant Pathogen") %>% 
  filter(phy == "Ascomycota" | phy == "Basidiomycota")
FRRN_taxa_FGanno_subgrp_pp

# Asc Bas split -------------
FRRN_taxa_FGanno_subgrp_pp <- 
  FRRN_taxa_FGanno_subgrp_pp %>% mutate(Guild3 = str_c(phy, Guild1, sep = "-"))


kru_pp_AscBas_rlt <- kruskal(FRRN_taxa_FGanno_subgrp_pp$Both_ITS_LSU,
                             FRRN_taxa_FGanno_subgrp_pp$Guild3, p.adj = "fdr")
kru_pp_AscBas_rlt

FRRN_FG_pp_AscBas_maxrDNA <- 
  FRRN_taxa_FGanno_subgrp_pp %>% group_by(Guild3) %>% summarise(max_rDNA = max(Both_ITS_LSU))

FRRN_FG_pp_AscBas_plot_subdata1 <- 
  kru_pp_AscBas_rlt$groups %>% as.data.frame() %>% 
  mutate(Guild3 = rownames(.), .before = groups) %>%
  select(Guild3, groups) %>% left_join(FRRN_FG_pp_AscBas_maxrDNA, by = "Guild3")

FRRN_FG_pp_AscBas_plot_subdata1

FRRN_FG_pp_AscBas_plot_subdata2 <- 
  FRRN_taxa_FGanno_subgrp_pp %>% group_by(Guild3) %>% count(Guild3)
FRRN_FG_pp_AscBas_plot_subdata2

kruskal.test(Both_ITS_LSU ~ phy, data = FRRN_taxa_FGanno_subgrp_pp)

FRRN_FG_pp_AscBas <- 
  ggplot(FRRN_taxa_FGanno_subgrp_pp, aes(Guild3, Both_ITS_LSU)) +
  geom_half_violin(data = FRRN_taxa_FGanno_subgrp_pp,
                   side = "r", colour = NA, alpha = 0.8, aes(fill = phy)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = phy),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_FG_pp_AscBas_plot_subdata1,
            aes(x = Guild3, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "red",
            size = 5) +
  geom_text(data = FRRN_FG_pp_AscBas_plot_subdata2,
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
FRRN_FG_pp_AscBas


# Fungal traits database ------------------------
Fungaltrait_database <- read_excel("./2.database/Fungaltrait_table.xlsx", sheet = 1)
FT_Ee <- Fungaltrait_database %>% select(GENUS, Ectomycorrhiza_exploration_type_template...9) %>% drop_na()

FRRN_taxa_FGanno_subgrp_AMFECM
FRRN_taxa_FGanno_subgrp_ECM <- FRRN_taxa_FGanno_subgrp_AMFECM %>% 
  filter(Guild1 == "Ectomycorrhizal")

colnames(FT_Ee)[1] <- "gen"
colnames(FT_Ee)[2] <- "Ectomycorrhiza_exploration_type_template"

FRRN_ECM_FG_FT_Ee <- 
  FRRN_taxa_FGanno_subgrp_ECM %>% left_join(FT_Ee, by = "gen") %>% drop_na()

FRRN_ECM_FG_FT_Ee$Ectomycorrhiza_exploration_type_template %>% table()

FRRN_ECM_FG_FT_Ee_filterd <- 
  FRRN_ECM_FG_FT_Ee %>% 
  filter(Ectomycorrhiza_exploration_type_template != "mat") %>% 
  filter(Ectomycorrhiza_exploration_type_template != "unknown")

FRRN_ECM_FG_FT_Ee_filterd <- 
  FRRN_ECM_FG_FT_Ee_filterd %>% 
  mutate(Ect_exp_type = if_else(str_detect(Ectomycorrhiza_exploration_type_template, "short-distance"), "short-distance", 
                                Ectomycorrhiza_exploration_type_template))

FRRN_ECM_FG_FT_Ee_filterd$Ect_exp_type %>% table()

FRRN_ECM_FG_FT_Ee_filterd$Ect_exp_type <- 
  factor(FRRN_ECM_FG_FT_Ee_filterd$Ect_exp_type,
         levels = c("contact",
                    "short-distance",
                    "medium-distance_fringe","medium-distance_smooth",
                    "long-distance"))

# fig2c, ECM forage type -----------
kru_ECM <- kruskal(FRRN_ECM_FG_FT_Ee_filterd$Both_ITS_LSU, FRRN_ECM_FG_FT_Ee_filterd$Ect_exp_type, p.adj = "fdr")
kru_ECM$groups

FRRN_ECM_FG_FT_Ee_filterd

FRRN_ECM_maxrDNA <- 
  FRRN_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% summarise(max_rDNA = max(Both_ITS_LSU))


FRRN_ECM_FG_FT_Ee_filterd_subdata1 <- 
  kru_ECM$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type = rownames(.), .before = groups) %>% 
  select(Ect_exp_type, groups) %>% left_join(FRRN_ECM_maxrDNA, by = "Ect_exp_type")

FRRN_ECM_FG_FT_Ee_filterd_MeSd <- FRRN_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% 
  summarise(m_rDNA = mean(Both_ITS_LSU),
            sd_rDNA = sd(Both_ITS_LSU))

FRRN_ECM_FG_FT_Ee_filterd_subdata1 <- 
  FRRN_ECM_FG_FT_Ee_filterd_subdata1 %>% left_join(FRRN_ECM_FG_FT_Ee_filterd_MeSd, by = "Ect_exp_type") %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

FRRN_ECM_FG_FT_Ee_filterd_subdata2 <- 
  FRRN_ECM_FG_FT_Ee_filterd %>% group_by(Ect_exp_type) %>% count(Ect_exp_type) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

kruskal.test(Both_ITS_LSU ~ Ect_exp_type, data = FRRN_ECM_FG_FT_Ee_filterd)
# Kruskal-Wallis chi-squared = 31.114, df = 4, p-value = 2.902e-06

fig2c_FRRN_ECM_Ee <- 
  ggplot(FRRN_ECM_FG_FT_Ee_filterd, aes(Ect_exp_type, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_ECM_FG_FT_Ee_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 5, label = lab, vjust = -0.5),
            colour = "black",
            size = 4) +
  geom_text(data = FRRN_ECM_FG_FT_Ee_filterd_subdata1,
            aes(x = Ect_exp_type, y = max_rDNA + 20, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = FRRN_ECM_FG_FT_Ee_filterd_subdata2,
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
fig2c_FRRN_ECM_Ee

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2c_pdf <- str_c("fig2c_", "FRRN_ECM_Ee_", tm, ".pdf", sep = "")
fig2c_jpg <- str_c("fig2c_", "FRRN_ECM_Ee_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_2c_pdf <- str_c(fig_path, fig2c_pdf)
fig_fullpath_2c_jpg <- str_c(fig_path, fig2c_jpg)

ggsave(fig_fullpath_2c_pdf, fig2c_FRRN_ECM_Ee, width = 7.9, height = 4.74)
ggsave(fig_fullpath_2c_jpg, fig2c_FRRN_ECM_Ee, width = 7.9, height = 4.74)


# figS4, forage type between Asc and Bas ---------------------
FRRN_ECM_FG_FT_Ee_filterd_phy <-
  FRRN_ECM_FG_FT_Ee_filterd %>% 
  filter(phy != "Mucoromycota") %>% 
  filter(Ect_exp_type == "contact" | Ect_exp_type == "short-distance") %>% 
  mutate(Ect_exp_type1 = str_c(phy, Ect_exp_type, sep = "-"))

kru_ECM_phy <- kruskal(FRRN_ECM_FG_FT_Ee_filterd_phy$Both_ITS_LSU,
                       FRRN_ECM_FG_FT_Ee_filterd_phy$Ect_exp_type1, p.adj = "fdr")
kru_ECM_phy

FRRN_ECM_phy_maxrDNA <- 
  FRRN_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% summarise(max_rDNA = max(Both_ITS_LSU))

FRRN_ECM_FG_FT_Ee_filterd_phy_subdata1 <- 
  kru_ECM_phy$groups %>% as.data.frame() %>% 
  mutate(Ect_exp_type1 = rownames(.), .before = groups) %>% 
  select(Ect_exp_type1, groups) %>% left_join(FRRN_ECM_phy_maxrDNA, by = "Ect_exp_type1")

FRRN_ECM_FG_FT_Ee_filterd_phy_subdata2 <- 
  FRRN_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% count(Ect_exp_type1)

FRRN_ECM_FG_FT_Ee_filterd_phy_MeSd <- 
  FRRN_ECM_FG_FT_Ee_filterd_phy %>% group_by(Ect_exp_type1) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                                              sd_rDNA = sd(Both_ITS_LSU))
FRRN_ECM_FG_FT_Ee_filterd_phy_subdata2 <-
  FRRN_ECM_FG_FT_Ee_filterd_phy_subdata2 %>% left_join(FRRN_ECM_FG_FT_Ee_filterd_phy_MeSd, by = "Ect_exp_type1") %>% 
  mutate(lab = str_c("n = ", n, " ", "(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

FRRN_ECM_FG_FT_Ee_filterd_phy$phy %>% table()
FRRN_ECM_FG_FT_Ee_filterd_phy$Ect_exp_type1 %>% table()

scheirerRayHare(Both_ITS_LSU ~ phy + Ect_exp_type, data = FRRN_ECM_FG_FT_Ee_filterd_phy)

figS4_FRRN_ECM_Ee_AscBas <- 
  ggplot(FRRN_ECM_FG_FT_Ee_filterd_phy, aes(Ect_exp_type1, Both_ITS_LSU)) +
  geom_half_violin(side = "r", colour = NA, alpha = 0.8, aes(fill = Ect_exp_type)) +
  geom_half_boxplot(side = "r", errorbar.draw = F, width = 0.1, outlier.shape = NA) +
  geom_half_point_panel(side = "l", transformation = position_jitter(width = 0.3, seed = 100), 
                        range_scale = 1, aes(colour = Ect_exp_type),
                        size = 1,
                        alpha = 0.8) +
  geom_text(data = FRRN_ECM_FG_FT_Ee_filterd_phy_subdata1,
            aes(x = Ect_exp_type1, y = max_rDNA, label = groups, vjust = -0.5),
            colour = "blue",
            size = 6) +
  geom_text(data = FRRN_ECM_FG_FT_Ee_filterd_phy_subdata2,
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
                                                                  Forage ~~ type ~~ (F): H == "10.954;" ~~ df == "1;" ~~ italic(p) == "0.000933" ~~
                                                                  "P x Ft:" ~~ H == "10.757;" ~~ df == "1;" ~~ italic(p) == "0.00104"),
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
figS4_FRRN_ECM_Ee_AscBas

tm <- now() %>% str_split_i(pattern = " ", 1)
figS4_pdf <- str_c("figS4_", "FRRN_ECM_Ee_AscBas_", tm, ".pdf", sep = "")
figS4_jpg <- str_c("figS4_", "FRRN_ECM_Ee_AscBas_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S4_pdf <- str_c(fig_path, figS4_pdf)
fig_fullpath_S4_jpg <- str_c(fig_path, figS4_jpg)

ggsave(fig_fullpath_S4_pdf, figS4_FRRN_ECM_Ee_AscBas, width = 11, height = 5.42)
ggsave(fig_fullpath_S4_jpg, figS4_FRRN_ECM_Ee_AscBas, width = 11, height = 5.42)


# done...



# extra ---------
FRRN_taxa_FGanno_full1

FRRN_rDNA_GS_GN_Taxa_FG <- 
  FRRN_rlt_taxa_spl_GSGN %>% select(-Name.y) %>% left_join(FRRN_taxa_FGanno_full1 %>% 
                                                           select(project, Trophic_Mode, Guild, Guild1, 
                                                                  Growth_Morphology, Trait, Confidence_Ranking),
                                                         by = "project")


FRRN_rDNA_GS_GN_Taxa_FG %>% filter(is.na(Assembly_Length))


# busco completeness
# busco_summary %>% dim()

FRRN_rDNA_GS_GN_Taxa_FG_busco <- FRRN_rDNA_GS_GN_Taxa_FG %>% left_join(busco_summary, by = c("project" = "project_id"))

FRRN_rDNA_GS_GN_Taxa_FG_busco %>% colnames()

FRRN_db_v1 <- FRRN_rDNA_GS_GN_Taxa_FG_busco %>% select(project, Both_ITS_LSU, Assembly_Length, Genes, completeness.x, 
                                                       depth_avg, det_single, det_multi, diff,
                                                       Name.x, phy, cla, ord, fam, Gen, Classification,
                                                       Trophic_Mode, Guild, Growth_Morphology, Trait, Confidence_Ranking)
colnames(FRRN_db_v1) <- c("project", "RRN", "Assembly_Length", "Genes", "completeness",
                          "SingleCopy_gene_SeqencingDepth(avg)", "single_copy_seqDepth", "rDNA_seqDepth", "rDNA_dif",
                          "Species", "Phylum", "Class", "Order", "Family", "Genus",
                          "Classification", "Trophic_Mode", "Guild", "Growth_Morphology", "Trait", "Confidence_Ranking(Guild)")



# write_xlsx(FRRN_db_v1, "FRRN_db_v1.xlsx")




