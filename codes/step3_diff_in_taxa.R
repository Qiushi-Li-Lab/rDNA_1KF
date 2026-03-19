

####### fig1b, fig1c #######
# rDNA copy number variation in different fungal taxa group
# by Qiushi-Li, IMCAS
# 2024.11.10

# load
# load("rCNV_version_1.0_20250219.RData")

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

library(nlme)
library(lme4)
library(lmerTest)
library(MuMIn)


# data -----------------
rCNV_rlt_taxa_spl %>% view()

rCNV_rlt_taxa_spl$Both_ITS_LSU %>% range()
rCNV_rlt_taxa_spl$Both_ITS_LSU %>% mean()
rCNV_rlt_taxa_spl$Both_ITS_LSU %>% median()

rCNV_rlt_taxa_spl %>% filter(Both_ITS_LSU == 1) %>% view()

93 / 5.4
# fig1c phylum rDNA copy number ------------------

rCNV_rlt_taxa_spl_grp <- 
  rCNV_rlt_taxa_spl %>% group_by(phy) %>% count(phy) %>% 
  mutate(grp = if_else(n >= 3, "Main_taxa", "Other_taxa")) %>% 
  select(-n)

#rCNV_taxa_ord <- rCNV_rlt_taxa_spl_grp %>% arrange(grp, phy) %>% select(phy) %>% pull()
#rCNV_taxa_ord1 <- str_sub(rCNV_taxa_ord, start = 3)

rCNV_rlt_taxa_spl_grp

rCNV_rlt_taxa_spl1 <- 
  rCNV_rlt_taxa_spl %>% left_join(rCNV_rlt_taxa_spl_grp, by = "phy")

rCNV_rlt_taxa_spl1 <- rCNV_rlt_taxa_spl1 %>% 
  mutate(phy = str_sub(phy, start = 3))


rCNV_rlt_taxa_spl_main <- 
  rCNV_rlt_taxa_spl1 %>% filter(grp == "Main_taxa")


main_taxa_kru <- kruskal(rCNV_rlt_taxa_spl_main$Both_ITS_LSU, rCNV_rlt_taxa_spl_main$phy, p.adj = "fdr")
main_taxa_kru


main_taxa_kru_rlt <- 
  main_taxa_kru$groups %>% as.data.frame() %>% 
  mutate(phy = rownames(.), .before = groups) %>%
  select(phy, groups)
main_taxa_kru_rlt

main_taxa_rDNA_max <- rCNV_rlt_taxa_spl_main %>% group_by(phy) %>% summarise(max_rDNA = max(Both_ITS_LSU))
main_taxa_rDNA_max

main_taxa_rDNA_MeSd <- rCNV_rlt_taxa_spl_main %>% group_by(phy) %>% summarise(m_rDNA = mean(Both_ITS_LSU), sd_rDNA = sd(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(",round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

main_taxa_kru_plotdata <- main_taxa_kru_rlt %>% left_join(main_taxa_rDNA_max, by = "phy") %>% left_join(main_taxa_rDNA_MeSd, by = "phy")
main_taxa_kru_plotdata

plot_subdata_Phy <- rCNV_rlt_taxa_spl1 %>% count(phy) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

rCNV_taxa_ord1 <- main_taxa_kru$groups %>% rownames() %>% rev()

rCNV_taxa_ord1 <- c(rCNV_taxa_ord1, "Blastocladiomycota", "Monoblepharomycota", "Basidiobolomycota")

rCNV_rlt_taxa_spl1$phy <- factor(rCNV_rlt_taxa_spl1$phy, levels = rCNV_taxa_ord1)

dim(rCNV_rlt_taxa_spl1)

rCNV_rlt_taxa_spl1$project %>% unique() %>% length()

#rCNV_rlt_taxa_spl1 %>% filter(phy == "Neocallimastigomycota")

#fix(main_taxa_kru_plotdata)

# > 500 were not shown
rCNV_rlt_taxa_spl1 %>% filter(Both_ITS_LSU > 500)


fig1c_rCNV_Phy_1 <- 
  rCNV_rlt_taxa_spl1 %>% 
  filter(grp == "Main_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                                 "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy")) +
  scale_x_discrete(limits = rCNV_taxa_ord1[1:10]) +
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
fig1c_rCNV_Phy_1

# rCNV_rlt_taxa_spl1$grp %>% table()

fig1c_rCNV_Phy_2 <- 
  rCNV_rlt_taxa_spl1 %>% 
  filter(grp == "Other_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_x_discrete(limits = rCNV_taxa_ord1[c(12, 11, 13)]) +
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
fig1c_rCNV_Phy_2

# help("plot_layout")

fig1c_dis <- "AAAAAAAAAABBB"

fig1c_rCNV_Phy1 <- fig1c_rCNV_Phy_1 + fig1c_rCNV_Phy_2 + plot_layout(design = fig1c_dis)
fig1c_rCNV_Phy1


tm <- now() %>% str_split_i(pattern = " ", 1)
fig1c_pdf <- str_c("fig1c_", "rCNV_phy_", tm, ".pdf", sep = "")
fig1c_jpg <- str_c("fig1c_", "rCNV_phy_", tm, ".jpg", sep = "")

ggsave(fig1c_pdf, fig1c_rCNV_Phy1, width = 9.86, height = 6.73)
ggsave(fig1c_jpg, fig1c_rCNV_Phy1, width = 9.86, height = 6.73)


# fig1b, rDNA copy number distribution among main taxa groups --------------- 
rCNV_rlt_taxa_spl2 <- rCNV_rlt_taxa_spl1 %>% 
  add_count(phy) %>% 
  mutate(phylum = if_else(n >= 3, phy, "Others"))

rCNV_rlt_taxa_spl2$phylum <- factor(rCNV_rlt_taxa_spl2$phylum, levels = c(
  "Ascomycota", "Basidiomycota", "Chytridiomycota", "Entomophthoromycota","Glomeromycota",  "Kickxellomycota",
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



rCNV_rlt_taxa_spl2$Both_ITS_LSU %>% median()

fig1b_total_distribution1 <- 
  ggplot(rCNV_rlt_taxa_spl2, aes(Both_ITS_LSU)) +
  geom_histogram(binwidth = 10, aes(fill = phylum, group = phylum), position = "stack") +
  geom_vline(aes(xintercept = median(Both_ITS_LSU)), colour = "black", linetype = 2, linewidth = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                               "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                               "Others" = "grey"),
                    ) +
  scale_x_continuous(limits = c(0, 300), expand = c(0.03, 0.03)) +
  labs(x = "rDNA copy number", y = "No. of fungi", fill = "Phylum", title = "rDNA copy number distribution among 1138 fungal sequencing projects") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 14, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "italic", size = 12),
        legend.key.size = unit(0.9, "cm"),
        legend.key.spacing.y = unit(0.35, units = "cm"),
        legend.position = "right")
fig1b_total_distribution1

fig1b_total_dis_sub <- 
  ggplot(rCNV_rlt_taxa_spl2, aes(Both_ITS_LSU)) +
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


ggsave(fig1b_pdf, fig1b_total_dis_fin, width = 10.4, height = 6.73)
ggsave(fig1b_jpg, fig1b_total_dis_fin, width = 10.4, height = 6.73)

# FigS2, distrubution between fungi and bateria ---------------------
rCNV_rlt_taxa_spl2

rrnDB_5.9 <- read_tsv("./database/rrnDB-5.9.tsv")
rrnDB_1 <- 
  rrnDB_5.9 %>% select(5, 6, 12, 13)

colnames(rrnDB_1) <- c("Data_source_organism_name", "NCBI_scientific_name", "gene_16S_count", "gene_23S_count")

rCNV_B <- 
  rrnDB_1 %>% select(Data_source_organism_name, gene_16S_count) %>% drop_na() %>% 
  mutate(
    taxa = "Bacteria"
  )
colnames(rCNV_B) <- c("ID", "rDNA_GCN", "taxa")

rCNV_F <- 
 rCNV_rlt_taxa_spl2 %>% select(project, Both_ITS_LSU) %>% mutate(
  taxa = "Fungi"
)

colnames(rCNV_F) <- c("ID", "rDNA_GCN", "taxa")
rCNV_FB_dist <- rbind(rCNV_B, rCNV_F)


# wow
rCNV_B %>% filter(rDNA_GCN == 21)

figS2_rCNV_FB_density_plot <- 
  ggplot(rCNV_FB_dist, aes(rDNA_GCN, colour = taxa)) +
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
figS2_rCNV_FB_density_plot

# fig1a ------------------------
rCNV_FB_dist

rCNV_FB_dist_kru <- kruskal.test(rDNA_GCN ~ taxa, data = rCNV_FB_dist)
rCNV_FB_dist_kru$statistic
# Kruskal-Wallis chi-squared = 3181.3, df = 1, p-value < 2.2e-16

rCNV_FB_dist_subdata1 <- rCNV_FB_dist %>% 
  group_by(taxa) %>% summarise(m_rDNA = mean(rDNA_GCN),
                                 sd_rDNA = sd(rDNA_GCN),
                                 max_rDNA = max(rDNA_GCN)) %>% 
  mutate(lab = str_c("(", round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))
rCNV_FB_dist_subdata1

rCNV_FB_dist_subdata2 <- rCNV_FB_dist %>% 
  group_by(taxa) %>% count(taxa) %>% mutate(lab = str_c("n =", n, sep = " ")) %>% 
  left_join(rCNV_FB_dist_subdata1, by = "taxa") %>% 
  mutate(lab = str_c(lab.x, lab.y, sep = "\n"))

fig1a_FB_kru <- 
  ggplot(rCNV_FB_dist, aes(taxa, rDNA_GCN)) +
  stat_boxplot(geom = "errorbar", aes(colour = taxa), width = 0.3, linewidth = 0.6) +
  geom_boxplot(aes(colour = taxa), linewidth = 1, width = 0.6, 
               outliers = F, alpha = 1, outlier.size = 0.5, outlier.colour = "black") +
  geom_text(data = rCNV_FB_dist_subdata2,
            aes(x = taxa, y = -15, label = lab), colour = "black", size = 3.5) + 
  annotate(geom = "text", x = 1.5, y = 265, label = expression("chi-square" == "3181.303"), size = 5) +
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

help("geom_boxplot")

tm <- now() %>% str_split_i(pattern = " ", 1)
fig1a_pdf <- str_c("fig1a_", "rCNV_FB_kru_", tm, ".pdf", sep = "")
fig1a_jpg <- str_c("fig1a_", "rCNV_FB_kru_", tm, ".jpg", sep = "")

ggsave(fig1a_pdf, fig1a_FB_kru, width = 3.21, height = 6.73)
ggsave(fig1a_jpg, fig1a_FB_kru, width = 3.21, height = 6.73)


# IQR-fig ----------------------
# fig1c, IQR ------------------
IQR_outlier <- function(val) {
  
  Q1 <- quantile(val, 0.25)
  Q3 <- quantile(val, 0.75)
  IQR <- Q3 - Q1
  
  # filtered
  outliers_ind <- (val < (Q1 - 1.5 * IQR) | val > (Q3 + 1.5 * IQR))
  
  return(outliers_ind)
  
}


rCNV_rlt_taxa_spl_IQR <- 
  rCNV_rlt_taxa_spl %>% group_by(phy) %>% 
  filter(!IQR_outlier(Both_ITS_LSU)) %>% 
  ungroup()


rCNV_rlt_taxa_spl_IQR_grp <- 
  rCNV_rlt_taxa_spl_IQR %>% group_by(phy) %>% count(phy) %>% 
  mutate(grp = if_else(n >= 3, "Main_taxa", "Other_taxa")) %>% 
  select(-n)

rCNV_rlt_taxa_spl_IQR_grp

rCNV_rlt_taxa_spl_IQR1 <- 
  rCNV_rlt_taxa_spl_IQR %>% left_join(rCNV_rlt_taxa_spl_IQR_grp, by = "phy")

rCNV_rlt_taxa_spl_IQR1 <- rCNV_rlt_taxa_spl_IQR1 %>% 
  mutate(phy = str_sub(phy, start = 3))


rCNV_rlt_taxa_spl_IQR_main <- 
  rCNV_rlt_taxa_spl_IQR1 %>% filter(grp == "Main_taxa")


main_taxa_kru_IQR <- kruskal(rCNV_rlt_taxa_spl_IQR_main$Both_ITS_LSU, rCNV_rlt_taxa_spl_IQR_main$phy, p.adj = "fdr")
main_taxa_kru_IQR

main_taxa_kru_IQR_rlt <- 
  main_taxa_kru_IQR$groups %>% as.data.frame() %>% 
  mutate(phy = rownames(.), .before = groups) %>%
  select(phy, groups)
main_taxa_kru_IQR_rlt

main_taxa_IQR_rDNA_max <- rCNV_rlt_taxa_spl_IQR_main %>% group_by(phy) %>% summarise(max_rDNA = max(Both_ITS_LSU))
main_taxa_IQR_rDNA_max

main_taxa_IQR_rDNA_MeSd <- rCNV_rlt_taxa_spl_IQR_main %>% group_by(phy) %>% summarise(m_rDNA = mean(Both_ITS_LSU), sd_rDNA = sd(Both_ITS_LSU)) %>% 
  mutate(lab = str_c("(",round(m_rDNA, 0), "±", round(sd_rDNA, 0), ")"))

main_taxa_kru_IQR_plotdata <- main_taxa_kru_IQR_rlt %>% left_join(main_taxa_IQR_rDNA_max, by = "phy") %>% 
  left_join(main_taxa_IQR_rDNA_MeSd, by = "phy")
main_taxa_kru_IQR_plotdata

plot_subdata_Phy_IQR <- rCNV_rlt_taxa_spl_IQR1 %>% count(phy) %>% 
  mutate(lab = str_c("n =", n, sep = " "))

rCNV_taxa_IQR_ord1 <- main_taxa_kru_IQR$groups %>% rownames() %>% rev()

rCNV_taxa_IQR_ord1 <- c(rCNV_taxa_IQR_ord1, "Blastocladiomycota", "Monoblepharomycota", "Basidiobolomycota")

rCNV_rlt_taxa_spl_IQR1$phy <- factor(rCNV_rlt_taxa_spl_IQR1$phy, levels = rCNV_taxa_IQR_ord1)

dim(rCNV_rlt_taxa_spl_IQR1)

#rCNV_rlt_taxa_spl_IQR1$project %>% unique() %>% length()
#rCNV_rlt_taxa_spl_IQR1 %>% filter(phy == "Neocallimastigomycota")


fig1c_rCNV_Phy_1_IQR <- 
  rCNV_rlt_taxa_spl_IQR1 %>% 
  filter(grp == "Main_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_colour_manual(values = c("#ff00ff","turquoise", "deepskyblue", "gold", "tomato", "pink", "yellowgreen","purple", "black", "navy")) +
  scale_x_discrete(limits = rCNV_taxa_IQR_ord1[1:10]) +
  stat_boxplot(geom = "errorbar", aes(colour = phy), width = 0.3, size = 0.6) +
  geom_boxplot(aes(colour = phy), size = 1, width = 0.6, outliers = F, alpha = 1) +
  geom_jitter(aes(colour = phy), size = 1.5, width = 0.15, alpha = 0.2) +
  geom_text(data = plot_subdata_Phy_IQR,
            aes(x = phy, y = 0, label = lab), colour = "black", show.legend = F, vjust = 2, hjust = 0.5, size = 4) +
  geom_text(data = main_taxa_kru_IQR_plotdata, aes(phy, max_rDNA + 20, label = groups), vjust = -0.5, colour = "blue", size = 5.5) +
  geom_text(data = main_taxa_kru_IQR_plotdata, aes(phy, max_rDNA + 5, label = lab), vjust = -0.5, colour = "black", size = 3.5) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(x = NULL, y = "rDNA copy number", fill = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 13, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold", size = 13), 
        axis.text = element_text(colour = "black"),
        panel.background = element_rect(colour = "black"),
        axis.title = element_text(face = "bold", size = 15),
        legend.position = "none")
fig1c_rCNV_Phy_1_IQR


fig1c_rCNV_Phy_2_IQR <- 
  rCNV_rlt_taxa_spl_IQR1 %>% 
  filter(grp == "Other_taxa") %>% 
  ggplot(aes(phy, Both_ITS_LSU)) +
  #scale_colour_simpsons() +
  scale_x_discrete(limits = rCNV_taxa_IQR_ord1[11:13]) +
  stat_boxplot(geom = "errorbar", colour = "grey", width = 0.3, size = 0.6) +
  geom_boxplot(colour = "grey", size = 1, width = 0.6, outliers = F, alpha = 1) +
  geom_jitter(colour = "grey", size = 1.5, width = 0.15, alpha = 0.2) +
  geom_text(data = plot_subdata_Phy_IQR,
            aes(x = phy, y = 0, label = lab), colour = "black", show.legend = F, vjust = 2, hjust = 0.5, size = 4) +
  geom_text(data = main_taxa_kru_IQR_plotdata, aes(phy, max_rDNA, label = groups), vjust = -0.5, colour = "blue", size = 5.5) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(x = NULL, y = NULL, fill = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", size = 13, angle = 50, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold", size = 13), 
        axis.text = element_text(colour = "black"),
        panel.background = element_rect(colour = "black"),
        axis.title = element_text(face = "bold", size = 15),
        legend.position = "none")
fig1c_rCNV_Phy_2_IQR

# help("plot_layout")

fig1c_dis <- "AAAAAAAAAABBB"

fig1c_rCNV_Phy1_IQR <- fig1c_rCNV_Phy_1_IQR + fig1c_rCNV_Phy_2_IQR + plot_layout(design = fig1c_dis)
fig1c_rCNV_Phy1_IQR

tm <- now() %>% str_split_i(pattern = " ", 1)
fig1c_IQR_pdf <- str_c("fig1c_", "rCNV_phy_IQR_", tm, ".pdf", sep = "")
fig1c_IQR_jpg <- str_c("fig1c_", "rCNV_phy_IQR_", tm, ".jpg", sep = "")


ggsave(fig1c_IQR_pdf, fig1c_rCNV_Phy1_IQR, width = 10.9, height = 7.16)
ggsave(fig1c_IQR_jpg, fig1c_rCNV_Phy1_IQR, width = 10.9, height = 7.16)


# fig1c, genome size ------------------

# fungi list
fungi_lst <- read_excel("./JGI_F_proj/JGI_totalstatus_20250223.xlsx", sheet = 1)

fungi_GS_GN <- fungi_lst %>% select(Project_ID, Name, Assembly_Length, Genes) %>% 
  drop_na() %>% distinct_all()
fungi_GS_GN


# rDNA copy number ----------------------------
rCNV_rlt_taxa_spl2

# new
rCNV_rlt_taxa_spl_GSGN <- rCNV_rlt_taxa_spl2 %>% left_join(fungi_GS_GN, by = c("project" = "Project_ID"))
rCNV_rlt_taxa_spl_GSGN

# fig ------
phy_lst <- rCNV_rlt_taxa_spl_GSGN %>% group_by(phy) %>% count(phy) %>% filter(n >= 10) %>% select(phy) %>% pull()

rCNV_rlt_taxa_spl_GSGN1 <- rCNV_rlt_taxa_spl_GSGN %>% filter(phy %in% phy_lst)


#c(1, 4, 15, 1, 16, 17, 18, 19)

rCNV_GS_lm <- rCNV_rlt_taxa_spl_GSGN1 %>% select(Both_ITS_LSU, Assembly_Length, phy) %>% 
  mutate(Both_ITS_LSU = log10(Both_ITS_LSU),
         Assembly_Length = log10(Assembly_Length))

lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(phy == subgrp)
    
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
      phy = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + 
        range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(log10(df_tmp[[yval]]), na.rm = T)[2] - range(log10(df_tmp[[yval]]), na.rm = T)[1]) * 200 +
        range(log10(df_tmp[[yval]]), na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p,
      pval1 = lmR_P_tmp
    )
    
  } else {
    
    df_rlt <- data.frame(
      phy = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.9 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p,
      pval1 = lmR_P_tmp
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  return(df_rlt_F)
  
}

phy_list <- rCNV_GS_lm$phy %>% unique()


rCNV_GS_lm_fit1 <- 
  map_dfr(phy_list, ~ lm_rlt_Spcs(rCNV_GS_lm,
                                           yval = "Assembly_Length",
                                           xval = "Both_ITS_LSU",
                                           subgrp = .,
                                  trans = "no_trans"))

#rCNV_GS_lm_fit

#rDNA_GS_lmer <- lmer(log10(Assembly_Length) ~ log10(Both_ITS_LSU) + (1 + log10(Both_ITS_LSU) | phy), data = rCNV_rlt_taxa_spl_GSGN1)
#rDNA_GS_lmer %>% summary()

#anova(rDNA_GS_lmer)
#r.squaredGLMM(rDNA_GS_lmer)

#rCNV_rlt_taxa_spl_GSGN1 %>% view()


# facet -------------


rCNV_GS_lm_fit1

rCNV_rlt_taxa_spl_GSGN1$phy <- factor(rCNV_rlt_taxa_spl_GSGN1$phy,
                                      levels = c(
                                        "Glomeromycota",
                                        "Ascomycota",
                                        "Basidiomycota",
                                        "Mucoromycota",
                                        "Chytridiomycota",
                                        "Mortierellomycota"))

rCNV_rlt_taxa_spl_GSGN2 <- 
  rCNV_rlt_taxa_spl_GSGN1 %>% mutate(
    Both_ITS_LSU = log10(Both_ITS_LSU),
    Assembly_Length = log10(Assembly_Length)
  )

fig1d_facet_rCNV_GS <- 
  ggplot(rCNV_rlt_taxa_spl_GSGN2, aes(Both_ITS_LSU, Assembly_Length, colour = phy)) +
  geom_point(size = 3.5, alpha = 0.3) +
  geom_smooth(method = "lm", se = T, linewidth = 1, show.legend = F, alpha = 0.3, colour = "blue") +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff", "Zoopagomycota" = "turquoise", "Ascomycota" = "deepskyblue", "Kickxellomycota" = "gold", "Mucoromycota" = "tomato", 
                                 "Basidiomycota" = "pink", "Chytridiomycota" = "yellowgreen","Neocallimastigomycota" = "purple", "Mortierellomycota" = "black", "Entomophthoromycota" = "navy",
                                 "Blastocladiomycota" = "grey",
                                 "Monoblepharomycota" = "grey",
                                 "Basidiobolomycota" = "grey"),
  ) +
  geom_text(data = rCNV_GS_lm_fit1,
            aes(x = anno_x1, y = anno_y1 * 1, label = rlt_sig),
            parse = T,
            colour = "red",
            size = 5.5) +
  scale_shape_manual(values = c(15, 0, 16, 1, 17, 2, 3, 4)) +
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
      phy == "Mucoromycota" ~ scale_x_continuous(labels = floor(10 ^ seq(0, 2, 1)),
                                                 breaks = seq(0, 2, 1)),
      phy == "Chytridiomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1.6, 2.4, 0.4)),
                                                    breaks = seq(1.6, 2.4, 0.4)),
      phy == "Mortierellomycota" ~ scale_x_continuous(labels = floor(10 ^ seq(1, 2.5, 0.5)),
                                                      breaks = seq(1, 2.5, 0.5))
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
fig1d_facet_rCNV_GS



tm <- now() %>% str_split_i(pattern = " ", 1)
fig1d_facet_pdf <- str_c("fig1d_", "rCNV_GS_facet_", tm, ".pdf", sep = "")
fig1d_facet_jpg <- str_c("fig1d_", "rCNV_GS_facet_", tm, ".jpg", sep = "")

ggsave(fig1d_facet_pdf, fig1d_facet_rCNV_GS, width = 21.1, height = 3.89)
ggsave(fig1d_facet_jpg, fig1d_facet_rCNV_GS, width = 21.1, height = 3.89)


# done.







# ---------------------------------
rCNV_rlt_taxa_spl_GSGN
fig1d_rCNV_GS <- 
  ggplot(rCNV_rlt_taxa_spl_GSGN1, aes(Both_ITS_LSU %>% log10(), Assembly_Length %>% log10(), colour = phy)) +
  geom_point(shape = 1, size = 2, alpha = 0.7) +
  geom_smooth(aes(colour = phy), method = "lm", se = T, linewidth = 1.5, show.legend = F, alpha = 0.1) +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff",
                                 "Ascomycota" = "deepskyblue",
                                 "Kickxellomycota" = "gold",
                                 "Basidiomycota" = "pink", 
                                 "Mucoromycota" = "tomato",
                                 "Chytridiomycota" = "yellowgreen",
                                 "Neocallimastigomycota" = "purple",
                                 "Mortierellomycota" = "black"),
                      labels = c("Glomeromycota (n = 45, r = 0.446, p = 2.111e-03**)",
                                 "Ascomycota (n = 473, r = 0.204, p = 8.321e-06***)",
                                 "Kickxellomycota (n = 7, r = 0.725, p = 0.0651 NS)",
                                 "Basidiomycota (n = 449, r = 0.222, p = 2.033e-06***)",
                                 "Mucoromycota (n = 70, r = -0.276, p = 0.0207*)",
                                 "Chytridiomycota (n = 20, r = 0.478, p = 0.0331*)",
                                 "Neocallimastigomycota (n = 8, r = 0.229, p = 0.586 NS)",
                                 "Mortierellomycota (n = 53, r = 0.038, p = 0.787 NS)")
  ) +
  annotate(geom = "text", x = 1.5, y = 9,
           label = expression(F == "4.339," ~~ italic(p) == "0.0843;" ~~ Conditional ~~ R^2 == "0.763"),
           size = 4.5,
           colour = "blue",
           parse = T) +
  scale_x_continuous(labels = floor(10 ^ seq(0, 3, 1))) +
  scale_y_continuous(labels = floor(10 ^ seq(7, 9, 0.5) / 1000000)) +
  guides(colour = guide_legend(title = "Phylum",
                               ncol = 1, 
                               keyspacing = 1, 
                               override.aes = list(alpha = 1, size = 6, shape = 16))
  ) +
  labs(x = "rDNA copy number (log10-transformed)", y = "Assembly genome size (Mbp, log10-transformed)") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14, face = "italic"),
    legend.key.spacing.y = unit(0.5, units = "cm"),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
fig1d_rCNV_GS

rCNV_GS_lm_fit
rCNV_rlt_taxa_spl_GSGN1 %>% group_by(phy) %>% count(phy)

fig1d_rCNV_GS <- 
  ggplot(rCNV_rlt_taxa_spl_GSGN1, aes(Both_ITS_LSU %>% log10(), Assembly_Length %>% log10(), colour = phy)) +
  geom_point(shape = 1, size = 2, alpha = 0.7) +
  geom_smooth(aes(colour = phy), method = "lm", se = T, linewidth = 1.5, show.legend = F, alpha = 0.1) +
  scale_colour_manual(values = c("Glomeromycota" = "#ff00ff",
                                 "Ascomycota" = "deepskyblue",
                                 "Kickxellomycota" = "gold",
                                 "Basidiomycota" = "pink", 
                                 "Mucoromycota" = "tomato",
                                 "Chytridiomycota" = "yellowgreen",
                                 "Neocallimastigomycota" = "purple",
                                 "Mortierellomycota" = "black"),
                      labels = c("Glomeromycota (n = 45, r = 0.446, p = 2.111e-03**)",
                                 "Ascomycota (n = 473, r = 0.204, p = 8.321e-06***)",
                                 "Kickxellomycota (n = 7, r = 0.725, p = 0.0651 NS)",
                                 "Basidiomycota (n = 449, r = 0.222, p = 2.033e-06***)",
                                 "Mucoromycota (n = 70, r = -0.276, p = 0.0207*)",
                                 "Chytridiomycota (n = 20, r = 0.478, p = 0.0331*)",
                                 "Neocallimastigomycota (n = 8, r = 0.229, p = 0.586 NS)",
                                 "Mortierellomycota (n = 53, r = 0.038, p = 0.787 NS)")
  ) +
  annotate(geom = "text", x = 1.5, y = 9,
           label = expression(F == "4.339," ~~ italic(p) == "0.0843;" ~~ Conditional ~~ R^2 == "0.763"),
           size = 4.5,
           colour = "blue",
           parse = T) +
  scale_x_continuous(labels = floor(10 ^ seq(0, 3, 1))) +
  scale_y_continuous(labels = floor(10 ^ seq(7, 9, 0.5) / 1000000)) +
  guides(colour = guide_legend(title = "Phylum",
                               ncol = 1, 
                               keyspacing = 1, 
                               override.aes = list(alpha = 1, size = 6, shape = 16))
  ) +
  labs(x = "rDNA copy number (log10-transformed)", y = "Assembly genome size (Mbp, log10-transformed)") +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14, face = "italic"),
    legend.key.spacing.y = unit(0.5, units = "cm"),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
fig1d_rCNV_GS

tm <- now() %>% str_split_i(pattern = " ", 1)
fig1d_pdf <- str_c("fig1d_", "rCNV_GS_", tm, ".pdf", sep = "")
fig1d_jpg <- str_c("fig1d_", "rCNV_GS_", tm, ".jpg", sep = "")

ggsave(fig1d_pdf, fig1d_rCNV_GS, width = 11, height = 5.23)
ggsave(fig1d_jpg, fig1d_rCNV_GS, width = 11, height = 5.23)
