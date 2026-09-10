

####### fig1b, fig1c #######
# rDNA copy number variation in different fungal taxa group
# by Qiushi-Li, IMCAS
# 2024.11.10


# packages
library(tidyverse)
library(agricolae)

library(readxl)
library(writexl)

# library(ggsignif)
library(ggsci)

library(ggh4x)
library(ggrepel)
library(ggbreak)
library(gghalves)
library(patchwork)

# library(nlme)
# library(lme4)
# library(lmerTest)
# library(MuMIn)


# data -----------------
# FRRN_rlt_taxa_spl %>% view()



# FRRN_rlt_taxa_spl1 %>% filter(project == "Endsp1") %>% view()

# remove Endsp1 because of low completeness and diff-seq-depth among each marker gene
FRRN_rlt_taxa_spl1a <- FRRN_rlt_taxa_spl1 %>% filter(project != "Endsp1")

FRRN_rlt_taxa_spl1a$Both_ITS_LSU %>% range()
FRRN_rlt_taxa_spl1a$Both_ITS_LSU %>% mean()
FRRN_rlt_taxa_spl1a$Both_ITS_LSU %>% median()
# range 1 to 1914
# mean 92.34745
# median 69

# 92 / 5.4

# FRRN_rlt_taxa_spl1a$phy %>% table()

FRRN_rlt_taxa_spl1a %>% filter(
  phy == "Glomeromycota"
) %>% select(phy, Both_ITS_LSU) %>% 
  group_by(phy) %>% 
  summarise(
    Glo_rrn_mean = mean(Both_ITS_LSU),
    Glo_rrn_median = median(Both_ITS_LSU),
    Glo_rrn_max = max(Both_ITS_LSU),
    Glo_rrn_min = min(Both_ITS_LSU)
  )

# FRRN_rlt_taxa_spl1a$phy %>% table()

# fig1c phylum rDNA copy number ------------------

FRRN_rlt_taxa_spl_grp <- 
  FRRN_rlt_taxa_spl1a %>% group_by(phy) %>% count(phy) %>% 
  mutate(grp = if_else(n >= 3, "Main_taxa", "Other_taxa")) %>% 
  select(-n)

FRRN_rlt_taxa_spl_grp

FRRN_rlt_taxa_spl1a <- 
  FRRN_rlt_taxa_spl1a %>% left_join(FRRN_rlt_taxa_spl_grp, by = "phy")

FRRN_rlt_taxa_spl1a <- FRRN_rlt_taxa_spl1a %>% 
  mutate(phy = str_sub(phy, start = 3))


FRRN_rlt_taxa_main <- 
  FRRN_rlt_taxa_spl1a %>% filter(grp == "Main_taxa")


main_taxa_kru <- kruskal(FRRN_rlt_taxa_main$Both_ITS_LSU, FRRN_rlt_taxa_main$phy, p.adj = "fdr")
main_taxa_kru


main_taxa_kru_rlt <- 
  main_taxa_kru$groups %>% as.data.frame() %>% 
  mutate(phy = rownames(.), .before = groups) %>%
  select(phy, groups)
main_taxa_kru_rlt

main_taxa_rDNA_max <- FRRN_rlt_taxa_main %>% group_by(phy) %>% summarise(max_rDNA = max(Both_ITS_LSU))
main_taxa_rDNA_max

main_taxa_rDNA_MeSd <- FRRN_rlt_taxa_main %>% group_by(phy) %>% summarise(m_rDNA = mean(Both_ITS_LSU), sd_rDNA = sd(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(",round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

main_taxa_kru_plotdata <- main_taxa_kru_rlt %>% left_join(main_taxa_rDNA_max, by = "phy") %>% left_join(main_taxa_rDNA_MeSd, by = "phy")
main_taxa_kru_plotdata

plot_subdata_Phy <- FRRN_rlt_taxa_spl1a %>% count(phy) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

FRRN_taxa_ord1 <- main_taxa_kru$groups %>% rownames() %>% rev()

FRRN_taxa_ord1 <- c(FRRN_taxa_ord1, "Monoblepharomycota", "Basidiobolomycota")

FRRN_rlt_taxa_spl1a$phy <- factor(FRRN_rlt_taxa_spl1a$phy, levels = FRRN_taxa_ord1)

dim(FRRN_rlt_taxa_spl1a)

FRRN_rlt_taxa_spl1a$project %>% unique() %>% length()

#rCNV_rlt_taxa_spl1 %>% filter(phy == "Neocallimastigomycota")


# fix(main_taxa_kru_plotdata)

# > 500 were not shown
FRRN_rlt_taxa_spl1a %>% filter(Both_ITS_LSU > 500)


FRRN_rlt_taxa_spl1a %>% filter(phy == "Glomeromycota") %>% 
  select(project) %>% unique() %>% pull() %>% length()
# 45

fig1c_FRRN_Phy_1 <- 
  FRRN_rlt_taxa_spl1a %>% 
  filter(grp == "Main_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                                 "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                                 "Blastocladiomycota" = "red")) +
  scale_x_discrete(limits = FRRN_taxa_ord1[1:11]) +
  stat_boxplot(geom = "errorbar", aes(colour = phy), width = 0.3, linewidth = 0.6) +
  geom_boxplot(aes(colour = phy), linewidth = 1, width = 0.6, outliers = F, alpha = 1) +
  geom_jitter(aes(colour = phy), size = 1.5, width = 0.15, alpha = 0.2) +
  geom_text(data = plot_subdata_Phy,
            aes(x = phy, y = 0, label = lab), colour = "black", show.legend = F, vjust = 2, hjust = 0.5, size = 4) +
  geom_text(data = main_taxa_kru_plotdata, aes(phy, max_rDNA + 20, label = groups), vjust = -0.5, colour = "blue", size = 5.5) +
  geom_text(data = main_taxa_kru_plotdata, aes(phy, max_rDNA + 5, label = lab), vjust = -0.5, colour = "black", size = 3.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(-40, 510), breaks = seq(0, 510, 100), labels = seq(0, 510, 100)) +
  labs(x = NULL, y = "rDNA copy number", fill = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 13, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold", size = 13), 
        axis.text = element_text(colour = "black"),
        panel.background = element_rect(colour = "black"),
        axis.title = element_text(face = "bold", size = 15),
        legend.position = "none")
fig1c_FRRN_Phy_1

# rCNV_rlt_taxa_spl1$grp %>% table()

###### ing ######
fig1c_FRRN_Phy_2 <- 
  FRRN_rlt_taxa_spl1a %>% 
  filter(grp == "Other_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_x_discrete(limits = FRRN_taxa_ord1[c(12, 13)]) +
  stat_boxplot(geom = "errorbar", colour = "grey", width = 0.3, size = 0.6) +
  geom_boxplot(colour = "grey", size = 1, width = 0.6, outliers = F, alpha = 1) +
  geom_jitter(colour = "grey", size = 1.5, width = 0.15, alpha = 0.2) +
  geom_text(data = plot_subdata_Phy,
            aes(x = phy, y = 0, label = lab), colour = "black", show.legend = F, vjust = 2, hjust = 0.5, size = 4) +
  geom_text(data = main_taxa_kru_plotdata, aes(phy, max_rDNA, label = groups), vjust = -0.5, colour = "blue", size = 5.5) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(x = NULL, y = NULL, fill = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 13, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold", size = 13), 
        axis.text = element_text(colour = "black"),
        panel.background = element_rect(colour = "black"),
        axis.title = element_text(face = "bold", size = 15),
        legend.position = "none")
fig1c_FRRN_Phy_2

# help("plot_layout")

fig1c_dis <- "AAAAAAAAAABB"

fig1c_FRRN_Phy1 <- fig1c_FRRN_Phy_1 + fig1c_FRRN_Phy_2 + plot_layout(design = fig1c_dis)
fig1c_FRRN_Phy1

tm <- now() %>% str_split_i(pattern = " ", 1)
fig1c_pdf <- str_c("fig1c_", "FRRN_phy_", tm, ".pdf", sep = "")
fig1c_jpg <- str_c("fig1c_", "FRRN_phy_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_1c_pdf <- str_c(fig_path, fig1c_pdf)
fig_fullpath_1c_jpg <- str_c(fig_path, fig1c_jpg)

ggsave(fig_fullpath_1c_pdf, fig1c_FRRN_Phy1, width = 9.86, height = 6.73)
ggsave(fig_fullpath_1c_jpg, fig1c_FRRN_Phy1, width = 9.86, height = 6.73)


# fig1b, rDNA copy number distribution among main taxa groups --------------- 
FRRN_rlt_taxa_spl2 <- FRRN_rlt_taxa_spl1a %>% 
  add_count(phy) %>% 
  mutate(phylum = if_else(n >= 3, phy, "Others"))

FRRN_rlt_taxa_spl2$phylum <- factor(FRRN_rlt_taxa_spl2$phylum, levels = c(
  "Ascomycota", "Basidiomycota", "Blastocladiomycota", "Chytridiomycota", "Entomophthoromycota","Glomeromycota",  "Kickxellomycota",
  "Mortierellomycota", "Mucoromycota", "Neocallimastigomycota", "Zoopagomycota", "Others"
))

#fig1a_total_distribution <- 
#  ggplot(rCNV_rlt_taxa_spl2, aes(Both_ITS_LSU)) +
#  geom_histogram(binwidth = 10, aes(fill = phylum, group = phylum), position = "stack", alpha = 0.8, colour = "transparent") +
#  geom_vline(aes(xintercept = median(Both_ITS_LSU)), colour = "#FF007F", linetype = 2) +
#  guides(fill = guide_legend(ncol = 2), alpha = NULL) +
#  scale_fill_d3(palette = "category20") +
#  scale_x_continuous(limits = c(0, 2000), expand = c(0.03, 0.03)) +
#  labs(x = "rDNA copy number", y = "No. of fungi", fill = "Phylum") +
#  annotate("text", x = 980, y = 75,
#           label = "rDNA copy number distribution \n among 1104 fungal sequencing projects in MycoCosm",
#           fontface = "bold") +
#  theme_bw() +
#  theme(axis.title = element_text(face = "bold", size = 14, colour = "black"),
#        axis.text = element_text(size = 12, colour = "black"),
#        legend.title = element_text(face = "bold"),
#        legend.text = element_text(face = "italic"),
#        legend.position = "inside",
#        legend.position.inside = c(0.7, 0.35),
#        legend.background = element_rect(colour = "grey"),
#        aspect.ratio = 0.5)
#fig1a_total_distribution

#ggsave("fig1a_total_distribution_20250129_F.pdf", fig1a_total_distribution, units = "cm", width = 21.764, height = 11.781)
#ggsave("fig1a_total_distribution_20250129_F.jpg", fig1a_total_distribution, units = "cm", width = 21.764, height = 11.781)

# "#ff00ff","#00ff00", "deepskyblue", "gold", "red", "navy", "darkgreen","maroon3", "black", "bisque", "grey"

FRRN_rlt_taxa_spl2$Both_ITS_LSU %>% median()
# 69

fig1b_total_distribution1 <- 
  ggplot(FRRN_rlt_taxa_spl2, aes(Both_ITS_LSU)) +
  geom_histogram(binwidth = 10, aes(fill = phylum, group = phylum), position = "stack") +
  geom_vline(aes(xintercept = median(Both_ITS_LSU)), colour = "black", linetype = 2, linewidth = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                               "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                               "Blastocladiomycota" = "red", "Others" = "grey"),
                    ) +
  scale_x_continuous(limits = c(0, 300), expand = c(0.03, 0.03)) +
  labs(x = "rDNA copy number", y = "No. of fungi", fill = "Phylum", title = "rDNA copy number distribution among 1,157 fungal sequencing projects") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.key.size = unit(0.9, "cm"),
        legend.key.spacing.y = unit(0.3, units = "cm"),
        legend.position = "right")
fig1b_total_distribution1

fig1b_total_dis_sub <- 
  ggplot(FRRN_rlt_taxa_spl2, aes(Both_ITS_LSU)) +
  geom_histogram(binwidth = 10, aes(fill = phylum, group = phylum), position = "stack") +
  geom_vline(aes(xintercept = median(Both_ITS_LSU)), colour = "black", linetype = 2) +
  scale_fill_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                               "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                               "Others" = "grey")) +
  scale_x_continuous(limits = c(0, 2000), expand = c(0.03, 0.03)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 10, colour = "black"),
        axis.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        legend.position = "none")
fig1b_total_dis_sub


fig1b_total_dis_fin <- fig1b_total_distribution1 + inset_element(fig1b_total_dis_sub, left = 0.4, right = 0.99, bottom = 0.55, top = 0.99)
fig1b_total_dis_fin


tm <- now() %>% str_split_i(pattern = " ", 1)
fig1b_pdf <- str_c("fig1b_", "total_distribution_", tm, ".pdf", sep = "")
fig1b_jpg <- str_c("fig1b_", "total_distribution_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_1b_pdf <- str_c(fig_path, fig1b_pdf)
fig_fullpath_1b_jpg <- str_c(fig_path, fig1b_jpg)


ggsave(fig_fullpath_1b_pdf, fig1b_total_dis_fin, width = 10.4, height = 6.73)
ggsave(fig_fullpath_1b_jpg, fig1b_total_dis_fin, width = 10.4, height = 6.73)

# FigS2, distrubution between fungi and bateria ---------------------
FRRN_rlt_taxa_spl2

rrnDB_5.9 <- read_tsv("./2.database/rrnDB-5.9.tsv")
rrnDB_1 <- 
  rrnDB_5.9 %>% select(5, 6, 12, 13)

colnames(rrnDB_1) <- c("Data_source_organism_name", "NCBI_scientific_name", "gene_16S_count", "gene_23S_count")

RRN_B <- 
  rrnDB_1 %>% select(Data_source_organism_name, gene_16S_count) %>% drop_na() %>% 
  mutate(
    taxa = "Bacteria"
  )
colnames(RRN_B) <- c("ID", "rDNA_GCN", "taxa")

RRN_F <- 
 FRRN_rlt_taxa_spl2 %>% select(project, Both_ITS_LSU) %>% mutate(
  taxa = "Fungi"
)

colnames(RRN_F) <- c("ID", "rDNA_GCN", "taxa")
RRN_FB_dist <- rbind(RRN_B, RRN_F)

# wow
RRN_B %>% filter(rDNA_GCN == 21)

figS2_RRN_FB_density_plot <- 
  ggplot(RRN_FB_dist, aes(rDNA_GCN, colour = taxa)) +
  geom_density(adjust = 5, bw = 1) +
  annotate("text", x = 26, y = 0.04, label = "bacteria", colour = "red") +
  annotate("text", x = 58, y = 0.01, label = "fungi", colour = "blue") +
  scale_x_break(c(250, 1500), scales = 0.3) +
  scale_colour_manual(values = c("red", "blue")) +
  labs(x = "rDNA copy number", y = "Density", colour = "Taxonomy") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none"
    )
figS2_RRN_FB_density_plot

tm <- now() %>% str_split_i(pattern = " ", 1)
figS2_pdf <- str_c("figS2_", "FB_density_", tm, ".pdf", sep = "")
figS2_jpg <- str_c("figS2_", "FB_density_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S2_pdf <- str_c(fig_path, figS2_pdf)
fig_fullpath_S2_jpg <- str_c(fig_path, figS2_jpg)


ggsave(fig_fullpath_S2_pdf, figS2_RRN_FB_density_plot, width = 8.08, height = 3.51)
ggsave(fig_fullpath_S2_jpg, figS2_RRN_FB_density_plot, width = 8.08, height = 3.51)


# fig1a ------------------------
RRN_FB_dist

RRN_FB_dist_kru <- kruskal.test(rDNA_GCN ~ taxa, data = RRN_FB_dist)
RRN_FB_dist_kru$statistic
# Kruskal-Wallis chi-squared = 3233.461, df = 1, p-value < 2.2e-16

RRN_FB_dist_subdata1 <- RRN_FB_dist %>% 
  group_by(taxa) %>% summarise(m_rDNA = mean(rDNA_GCN),
                                 sd_rDNA = sd(rDNA_GCN),
                                 max_rDNA = max(rDNA_GCN)) %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))
RRN_FB_dist_subdata1

RRN_FB_dist_subdata2 <- RRN_FB_dist %>% 
  group_by(taxa) %>% count(taxa) %>% mutate(lab = str_c("n =", n, sep = " ")) %>% 
  left_join(RRN_FB_dist_subdata1, by = "taxa") %>% 
  mutate(lab = str_c(lab.x, lab.y, sep = "\n"))

fig1a_FB_kru <- 
  ggplot(RRN_FB_dist, aes(taxa, rDNA_GCN)) +
  stat_boxplot(geom = "errorbar", aes(colour = taxa), width = 0.3, linewidth = 0.6) +
  geom_boxplot(aes(colour = taxa), linewidth = 1, width = 0.6, 
               outliers = F, alpha = 1, outlier.size = 0.5, outlier.colour = "black") +
  geom_text(data = RRN_FB_dist_subdata2,
            aes(x = taxa, y = -15, label = lab), colour = "black", size = 3.5) + 
  annotate(geom = "text", x = 1.5, y = 265, label = expression("chi-square" == "3233.461"), size = 5) +
  annotate(geom = "text", x = 1.5, y = 250, label = expression("df" == "1;" ~~ italic(p) < "2.2e-16"), size = 5) +
  scale_colour_manual(values = c("red", "blue")) +
  scale_y_continuous(limits = c(-20, 270)) +
  labs(x = NULL, y = "rDNA copy number") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
fig1a_FB_kru

# help("geom_boxplot")

tm <- now() %>% str_split_i(pattern = " ", 1)
fig1a_pdf <- str_c("fig1a_", "RRN_FB_kru_", tm, ".pdf", sep = "")
fig1a_jpg <- str_c("fig1a_", "RRN_FB_kru_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_1a_pdf <- str_c(fig_path, fig1a_pdf)
fig_fullpath_1a_jpg <- str_c(fig_path, fig1a_jpg)


ggsave(fig_fullpath_1a_pdf, fig1a_FB_kru, width = 3.21, height = 6.73)
ggsave(fig_fullpath_1a_jpg, fig1a_FB_kru, width = 3.21, height = 6.73)


# fig1c, genome size ------------------

# fungi list
FRRN_proj

fungi_GS_GN <- FRRN_proj %>% select(Project_ID, Name, Assembly_Length, Genes) %>% 
  drop_na() %>% distinct_all()
fungi_GS_GN


# rDNA copy number ----------------------------
FRRN_rlt_taxa_spl2

# new
FRRN_rlt_taxa_spl_GSGN <- FRRN_rlt_taxa_spl2 %>% left_join(fungi_GS_GN, by = c("project" = "Project_ID"))
# FRRN_rlt_taxa_spl_GSGN %>% filter(is.na(Assembly_Length)) %>% view()

# fig ------
phy_lst <- FRRN_rlt_taxa_spl_GSGN %>% group_by(phy) %>% count(phy) %>% filter(n >= 10) %>% select(phy) %>% pull()
phy_n <- FRRN_rlt_taxa_spl_GSGN %>% group_by(phy) %>% count(phy) %>% filter(n >= 10) %>%
  mutate(lab_n = str_c("n = ", n, sep = ""))




FRRN_rlt_taxa_spl_GSGN1 <- FRRN_rlt_taxa_spl_GSGN %>% filter(phy %in% phy_lst)

#c(1, 4, 15, 1, 16, 17, 18, 19)

FRRN_GS_lm <- FRRN_rlt_taxa_spl_GSGN1 %>% select(Both_ITS_LSU, Assembly_Length, phy) %>% 
  mutate(Both_ITS_LSU = log10(Both_ITS_LSU),
         Assembly_Length = log10(Assembly_Length))


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


colnames(FRRN_GS_lm)[3] <- "group"

phy_list <- FRRN_GS_lm$group %>% unique()

FRRN_GS_cor_subdata <- 
  subplot_data_corr(DF = FRRN_GS_lm, Y_val = "Assembly_Length", X_val = "Both_ITS_LSU",
                    GROUP_list = phy_list, METHOD = "spearman", adj_METHOD = "fdr")
FRRN_GS_cor_subdata



#rCNV_GS_lm_fit

#rDNA_GS_lmer <- lmer(log10(Assembly_Length) ~ log10(Both_ITS_LSU) + (1 + log10(Both_ITS_LSU) | phy), data = rCNV_rlt_taxa_spl_GSGN1)
#rDNA_GS_lmer %>% summary()

#anova(rDNA_GS_lmer)
#r.squaredGLMM(rDNA_GS_lmer)

#rCNV_rlt_taxa_spl_GSGN1 %>% view()


# facet -------------

colnames(FRRN_GS_cor_subdata)[1] <- "phy"

FRRN_rlt_taxa_spl_GSGN1$phy <- factor(FRRN_rlt_taxa_spl_GSGN1$phy,
                                      levels = c(
                                        "Glomeromycota",
                                        "Ascomycota",
                                        "Basidiomycota",
                                        "Mucoromycota",
                                        "Chytridiomycota",
                                        "Mortierellomycota"))

FRRN_rlt_taxa_spl_GSGN2 <- 
  FRRN_rlt_taxa_spl_GSGN1 %>% mutate(
    Both_ITS_LSU = log10(Both_ITS_LSU),
    Assembly_Length = log10(Assembly_Length)
  )


# help("geom_smooth")

# rCNV_rlt_taxa_spl_GSGN2 %>% count(phy)



FRRN_GS_cor_subdata1 <- FRRN_GS_cor_subdata %>% 
  left_join(phy_n, by = "phy")

fig1d_facet_FRRN_GS <- 
  ggplot(FRRN_rlt_taxa_spl_GSGN2, aes(Both_ITS_LSU, Assembly_Length, colour = phy)) +
  geom_point(size = 3.5, alpha = 0.3) +
  geom_smooth(se = T, linewidth = 1, show.legend = F, alpha = 0.3, colour = "red", span = 5) +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                                 "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                                 "Blastocladiomycota" = "grey",
                                 "Monoblepharomycota" = "grey",
                                 "Basidiobolomycota" = "grey"),
  ) +
  geom_text(data = FRRN_GS_cor_subdata1,
            aes(x = anno_x1, y = Inf, label = rlt_adj_sig),
            parse = T,
            colour = "blue",
            size = 5.5,
            vjust = 3) +
  geom_text(data = FRRN_GS_cor_subdata1,
            aes(x = anno_x1, y = -Inf, label = lab_n),
            colour = "blue",
            size = 5.5,
            vjust = -2) +
  # scale_shape_manual(values = c(15, 0, 16, 1, 17, 2, 3, 4)) +
  facet_wrap(~ phy, scales = "free", nrow = 1) +
  guides(colour = guide_legend(title = "Phylum", ncol = 1, override.aes = list(alpha = 1, size = 2)),
         linetype = guide_legend(title = "Phylum"),
         shape = guide_legend(title = "Phylum")) +
  labs(x = "rDNA copy number (log10-transformed)", y = "Assembly genome size\n(Mbp, log10-transformed)") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect("white"),
    strip.text = element_text(colour = "black", size = 20, face = "bold.italic"),
    axis.title = element_text(face = "bold", size = 16, colour = "black"),
    axis.text = element_text(size = 14, colour = "black")
  ) +
  facetted_pos_scales(
    x = list(
      phy == "Glomeromycota" ~ scale_x_continuous(labels = floor(10 ^ seq(0, 1.5, 0.5)),
                                                  breaks = seq(0, 1.5, 0.5)),
      phy == "Ascomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1, 3, 1)),
                                               breaks = seq(1, 3, 1)),
      phy == "Basidiomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1, 3, 1)),
                                                  breaks = seq(1, 3, 1)),
      phy == "Mucoromycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1.5, 2.5, 0.5)),
                                                 breaks = seq(1.5, 2.5, 0.5)),
      phy == "Chytridiomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1.6, 2.4, 0.4)),
                                                    breaks = seq(1.6, 2.4, 0.4)),
      phy == "Mortierellomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1.9, 2.3, 0.2)),
                                                      breaks = seq(1.9, 2.3, 0.2))
    ),
    y = list(
      phy == "Glomeromycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.5, 8.5, 0.5) / 1000000),
                                                  breaks = seq(7.5, 8.5, 0.5)),
      phy == "Ascomycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.2, 8.0, 0.4) / 1000000),
                                               breaks = seq(7.2, 8.0, 0.4)),
      phy == "Basidiomycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.5, 9, 0.5) / 1000000),
                                                  breaks = seq(7.5, 9, 0.5)),
      phy == "Mucoromycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.5, 8.4, 0.3) / 1000000),
                                                 breaks = seq(7.5, 8.4, 0.3)),
      phy == "Chytridiomycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.2, 7.8, 0.2) / 1000000),
                                                    breaks = seq(7.2, 7.8, 0.2)),
      phy == "Mortierellomycota" ~ scale_y_continuous(labels = floor(10 ^ seq(7.5, 7.9, 0.1) / 1000000),
                                                      breaks = seq(7.5, 7.9, 0.1))
    )
  )
fig1d_facet_FRRN_GS



tm <- now() %>% str_split_i(pattern = " ", 1)
fig1d_facet_pdf <- str_c("fig1d_", "FRRN_GS_facet_", tm, ".pdf", sep = "")
fig1d_facet_jpg <- str_c("fig1d_", "FRRN_GS_facet_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_1d_pdf <- str_c(fig_path, fig1d_facet_pdf)
fig_fullpath_1d_jpg <- str_c(fig_path, fig1d_facet_jpg)


ggsave(fig_fullpath_1d_pdf, fig1d_facet_FRRN_GS, width = 21.1, height = 3.89)
ggsave(fig_fullpath_1d_jpg, fig1d_facet_FRRN_GS, width = 21.1, height = 3.89)



# done.
# FRRN_rlt_taxa_spl_GSGN1 %>% filter(phylum == "Mucoromycota") %>% filter(Assembly_Length > 125000000) %>% view()


# ---------------------------------
# rCNV_rlt_taxa_spl_GSGN1
# 
# fig1d_rCNV_GS <- 
#   ggplot(rCNV_rlt_taxa_spl_GSGN1, aes(Both_ITS_LSU %>% log10(), Assembly_Length %>% log10(), colour = phy)) +
#   geom_point(shape = 1, size = 2, alpha = 0.7) +
#   geom_smooth(aes(colour = phy), method = "lm", se = T, linewidth = 1.5, show.legend = F, alpha = 0.1) +
#   scale_colour_manual(values = c("Glomeromycota" = "#ff00ff",
#                                  "Ascomycota" = "deepskyblue",
#                                  "Kickxellomycota" = "gold",
#                                  "Basidiomycota" = "pink", 
#                                  "Mucoromycota" = "tomato",
#                                  "Chytridiomycota" = "yellowgreen",
#                                  "Neocallimastigomycota" = "purple",
#                                  "Mortierellomycota" = "black"),
#                       labels = c("Glomeromycota (n = 45, r = 0.446, p = 2.111e-03**)",
#                                  "Ascomycota (n = 473, r = 0.204, p = 8.321e-06***)",
#                                  "Kickxellomycota (n = 7, r = 0.725, p = 0.0651 NS)",
#                                  "Basidiomycota (n = 449, r = 0.222, p = 2.033e-06***)",
#                                  "Mucoromycota (n = 70, r = -0.276, p = 0.0207*)",
#                                  "Chytridiomycota (n = 20, r = 0.478, p = 0.0331*)",
#                                  "Neocallimastigomycota (n = 8, r = 0.229, p = 0.586 NS)",
#                                  "Mortierellomycota (n = 53, r = 0.038, p = 0.787 NS)")
#   ) +
#   annotate(geom = "text", x = 1.5, y = 9,
#            label = expression(F == "4.339," ~~ italic(p) == "0.0843;" ~~ Conditional ~~ R^2 == "0.763"),
#            size = 4.5,
#            colour = "blue",
#            parse = T) +
#   scale_x_continuous(labels = floor(10 ^ seq(0, 3, 1))) +
#   scale_y_continuous(labels = floor(10 ^ seq(7, 9, 0.5) / 1000000)) +
#   guides(colour = guide_legend(title = "Phylum",
#                                ncol = 1, 
#                                keyspacing = 1, 
#                                override.aes = list(alpha = 1, size = 6, shape = 16))
#   ) +
#   labs(x = "rDNA copy number (log10-transformed)", y = "Assembly genome size (Mbp, log10-transformed)") +
#   theme_bw() +
#   theme(
#     legend.position = "right",
#     legend.text = element_text(size = 14, face = "italic"),
#     legend.key.spacing.y = unit(0.5, units = "cm"),
#     legend.title = element_text(size = 16, face = "bold"),
#     axis.title = element_text(face = "bold", size = 14, colour = "black"),
#     axis.text = element_text(size = 12, colour = "black")
#   )
# fig1d_rCNV_GS
# 
# rCNV_GS_lm_fit
# rCNV_rlt_taxa_spl_GSGN1 %>% group_by(phy) %>% count(phy)
# 
# fig1d_rCNV_GS <- 
#   ggplot(rCNV_rlt_taxa_spl_GSGN1, aes(Both_ITS_LSU %>% log10(), Assembly_Length %>% log10(), colour = phy)) +
#   geom_point(shape = 1, size = 2, alpha = 0.7) +
#   geom_smooth(aes(colour = phy), method = "lm", se = T, linewidth = 1.5, show.legend = F, alpha = 0.1) +
#   scale_colour_manual(values = c("Glomeromycota" = "#ff00ff",
#                                  "Ascomycota" = "deepskyblue",
#                                  "Kickxellomycota" = "gold",
#                                  "Basidiomycota" = "pink", 
#                                  "Mucoromycota" = "tomato",
#                                  "Chytridiomycota" = "yellowgreen",
#                                  "Neocallimastigomycota" = "purple",
#                                  "Mortierellomycota" = "black"),
#                       labels = c("Glomeromycota (n = 45, r = 0.446, p = 2.111e-03**)",
#                                  "Ascomycota (n = 473, r = 0.204, p = 8.321e-06***)",
#                                  "Kickxellomycota (n = 7, r = 0.725, p = 0.0651 NS)",
#                                  "Basidiomycota (n = 449, r = 0.222, p = 2.033e-06***)",
#                                  "Mucoromycota (n = 70, r = -0.276, p = 0.0207*)",
#                                  "Chytridiomycota (n = 20, r = 0.478, p = 0.0331*)",
#                                  "Neocallimastigomycota (n = 8, r = 0.229, p = 0.586 NS)",
#                                  "Mortierellomycota (n = 53, r = 0.038, p = 0.787 NS)")
#   ) +
#   annotate(geom = "text", x = 1.5, y = 9,
#            label = expression(F == "4.339," ~~ italic(p) == "0.0843;" ~~ Conditional ~~ R^2 == "0.763"),
#            size = 4.5,
#            colour = "blue",
#            parse = T) +
#   scale_x_continuous(labels = floor(10 ^ seq(0, 3, 1))) +
#   scale_y_continuous(labels = floor(10 ^ seq(7, 9, 0.5) / 1000000)) +
#   guides(colour = guide_legend(title = "Phylum",
#                                ncol = 1, 
#                                keyspacing = 1, 
#                                override.aes = list(alpha = 1, size = 6, shape = 16))
#   ) +
#   labs(x = "rDNA copy number (log10-transformed)", y = "Assembly genome size (Mbp, log10-transformed)") +
#   theme_bw() +
#   theme(
#     legend.position = "right",
#     legend.text = element_text(size = 14, face = "italic"),
#     legend.key.spacing.y = unit(0.5, units = "cm"),
#     legend.title = element_text(size = 16, face = "bold"),
#     axis.title = element_text(face = "bold", size = 14, colour = "black"),
#     axis.text = element_text(size = 12, colour = "black")
#   )
# fig1d_rCNV_GS
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# fig1d_pdf <- str_c("fig1d_", "rCNV_GS_", tm, ".pdf", sep = "")
# fig1d_jpg <- str_c("fig1d_", "rCNV_GS_", tm, ".jpg", sep = "")
# 
# ggsave(fig1d_pdf, fig1d_rCNV_GS, width = 11, height = 5.23)
# ggsave(fig1d_jpg, fig1d_rCNV_GS, width = 11, height = 5.23)



# done.


