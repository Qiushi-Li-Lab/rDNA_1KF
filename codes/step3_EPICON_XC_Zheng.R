


##### rDNA copy number of fungi in EPICON
# 
# by Qiushi-Li, IM-CAS

library(tidyverse)
library(vegan)
library(broom)
library(readxl)
library(patchwork)
library(gghalves)
library(agricolae)
library(ggsci)

# load("rCNV_20241211.Rdata")

# EPICON
load("./database/EPICON.data.preparation.RC.bNTI.ted.2019.04.19.Rdata")

# Annotated
EPICON_fungi_anno <- read.csv("./database/All_otu_blast.csv", header = T)
EPICON_fungi_anno %>% colnames()

EPICON_fungi_PL <- EPICON_fungi_anno %>% 
  mutate(Class = str_sub(Class, start = 4),
         Order = str_sub(Order, start = 4),
         Family = str_sub(Family, start = 4),
         Genus = str_sub(Genus, start = 4),
         Subphylum = str_sub(Subphylum, start = 4))

# rDNA table
rCNV_rlt_taxa_spl
rCNV_taxa_FG

# explore --------------
dim(fung0)
fung0$Fungi %>% table()

EPICON_fungal <- fung0 %>% filter(Fungi == "Fungi")
dim(EPICON_fungal)


EPICON_fungal_otutab <- EPICON_fungal %>% select(1:1251)
dim(EPICON_fungal_otutab)

EPICON_fungal_otuExtra <- EPICON_fungal %>% select(1252:1302)

EPICON_fungal_otuExtra1 <- 
  EPICON_fungal_otuExtra %>% select(OTU, ID)
EPICON_fungal_otuExtra1

# spc tab
rCNV_spc_tab <- rCNV_rlt_taxa_spl %>% mutate(
  spc = str_c(str_split_i(Name, pattern = " ", 1), str_split_i(Name, pattern = " ", 2), sep = " ")
  ) %>% group_by(spc) %>% summarise(rDNA_m = mean(ITS_only))
rCNV_spc_tab

# gen tab
rCNV_gen_tab <- rCNV_rlt_taxa_spl %>% select(Gen, ITS_only) %>%
  group_by(Gen) %>% summarise(rDNA_m = mean(ITS_only)) %>%
  mutate(Gen = str_sub(Gen, start = 3))
rCNV_gen_tab

# fam tab
rCNV_fam_tab <- rCNV_rlt_taxa_spl %>% select(fam, ITS_only) %>%
  group_by(fam) %>% summarise(rDNA_m = mean(ITS_only)) %>%
  mutate(fam = str_sub(fam, start = 3))
rCNV_fam_tab

# ord tab
rCNV_ord_tab <- rCNV_rlt_taxa_spl %>% select(ord, ITS_only) %>%
  group_by(ord) %>% summarise(rDNA_m = mean(ITS_only)) %>%
  mutate(ord = str_sub(ord, start = 3))
rCNV_ord_tab

# cla tab
rCNV_cla_tab <- rCNV_rlt_taxa_spl %>% select(cla, ITS_only) %>%
  group_by(cla) %>% summarise(rDNA_m = mean(ITS_only)) %>%
  mutate(cla = str_sub(cla, start = 3))
rCNV_cla_tab

# phy tab
rCNV_phy_tab <- rCNV_rlt_taxa_spl %>% select(phy, ITS_only) %>%
  group_by(phy) %>% summarise(rDNA_m = mean(ITS_only)) %>%
  mutate(phy = str_sub(phy, start = 3))
rCNV_phy_tab


################### ---------

# FGanno -----------------------
# EPICON_fungi_PL %>% select(Order) %>% view()

FunGuild_database1

FGanno_gen <- FunGuild_database1 %>% filter(Taxon_Level == FG_index_df$Tax_Lev[1]) %>% select(-Taxon_Level)
colnames(FGanno_gen)[1] <- FG_index_df$Tax_Ind[1]

FGanno_EPICON_tmp <- EPICON_fungi_PL %>% 
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

# frDNA ---------------------------
EPICON_frDNA_gen <- EPICON_fungi_PL %>% left_join(rCNV_gen_tab, by = c("Genus" = "Gen")) %>% 
  select(ID, Genus, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% mutate(tax_lev = "Gen", .after = ID)
#EPICON_frDNA_gen
colnames(EPICON_frDNA_gen)[3] <- "tax"

EPICON_frDNA_anno <- EPICON_frDNA_gen

EPICON_frDNA_fam <- EPICON_fungi_PL %>% filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(rCNV_fam_tab, by = c("Family" = "fam")) %>% 
  select(ID, Family, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% mutate(tax_lev = "Fam", .after = ID)
#EPICON_frDNA_fam
colnames(EPICON_frDNA_fam)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_fam)

EPICON_frDNA_ord <- EPICON_fungi_PL %>% filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(rCNV_ord_tab, by = c("Order" = "ord")) %>% 
  select(ID, Order, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% mutate(tax_lev = "Ord", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_ord)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_ord)


EPICON_frDNA_cla <- EPICON_fungi_PL %>% filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(rCNV_cla_tab, by = c("Class" = "cla")) %>% 
  select(ID, Class, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% mutate(tax_lev = "Cla", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_cla)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_cla)


EPICON_frDNA_phy <- EPICON_fungi_PL %>% filter(!ID %in% EPICON_frDNA_anno$ID) %>% 
  left_join(rCNV_phy_tab, by = c("Subphylum" = "phy")) %>% 
  select(ID, Subphylum, rDNA_m) %>% filter(!is.na(rDNA_m)) %>% mutate(tax_lev = "Phy", .after = ID)
#EPICON_frDNA_ord
colnames(EPICON_frDNA_phy)[3] <- "tax"

EPICON_frDNA_anno <- rbind(EPICON_frDNA_anno, EPICON_frDNA_phy)
EPICON_frDNA_anno

# done ...

# otutab?
EPICON_fungal_otutab
colnames(EPICON_fungal_otutab)

# env
env

# rDNA anno
EPICON_frDNA_anno

# ok~ start ~
EPICON_fungal_otutab1 <- EPICON_fungal_otutab %>% select(env$aa)
EPICON_fungal_otutab1
# min abd is 362
EPICON_fungal_otutab1_r <- rrarefy(EPICON_fungal_otutab1 %>% t(), 362) %>% t()


apply(EPICON_fungal_otutab1_r, sum, MARGIN = 2)



rel_per <- function(val) {
  
  tmp_col <- val
  tmp_col_sum <- sum(tmp_col)
  
  re_per_col <- tmp_col / tmp_col_sum
  
  return(re_per_col)
  
}


EPICON_fungal_otutab1_per <- apply(EPICON_fungal_otutab1_r, rel_per, MARGIN = 2) %>% as.data.frame()
colnames(EPICON_fungal_otutab1_per)[1]

EPICON_fungal_otutab1_per1 <- 
  EPICON_fungal_otutab1_per %>% mutate(OTU_ID = rownames(EPICON_fungal_otutab1_per), .before = TP01L01) %>% 
  mutate(OTU_ID = str_split_i(OTU_ID, pattern = "_", 1)) %>% left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))


sample_id <- colnames(EPICON_fungal_otutab1_per)


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab1_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)

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

EPICON_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)

env_timepoint <- 
  data.frame(
    Timepiont = c("TP00", "TP01", "TP02", "TP03", "TP04", "TP05",
                  "TP06", "TP07", "TP08", "TP09", "TP10", "TP11",
                  "TP12", "TP13", "TP14", "TP15", "TP16", "TP17"),
    timepoint = c(0:17)
  )

env_timepoint

EPICON_rDNAm_env <- 
  EPICON_rDNAm %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")

EPICON_rDNAm_env$Timepiont %>% table()
EPICON_rDNAm_env$Habitat %>% table()
EPICON_rDNAm_env$Treatment %>% table()



# figS8, rCNV diff compartment -------

EPICON_com_kru <- 
  kruskal(EPICON_rDNAm_env$rDNAm, EPICON_rDNAm_env$Habitat, p.adj = "fdr")
EPICON_com_kru

EPICON_plot_subdata <- 
  EPICON_com_kru$groups %>% as.data.frame() %>% 
  mutate(Habitat = rownames(.), .before = groups) %>%
  select(Habitat, groups)
EPICON_plot_subdata

EPICON_com_maxrCNV <- EPICON_rDNAm_env %>% group_by(Habitat) %>% summarise(max_rDNA = max(rDNAm))

EPICON_plot_subdata1 <- EPICON_plot_subdata %>% left_join(EPICON_com_maxrCNV, by = "Habitat")

kruskal.test(rDNAm ~ Habitat, data = EPICON_rDNAm_env) %>% tidy()

figS8_EPICON_com <- 
  ggplot(EPICON_rDNAm_env, aes(Habitat, rDNAm)) +
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
  annotate(geom = "text", x = 2.5, y = 125, label = "chi-square = 414.42, df = 3, p = 1.67e-89", size = 5) +
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
figS8_EPICON_com

tm <- now() %>% str_split_i(pattern = " ", 1)
figS8_pdf <- str_c("figS8_", "EPICON_com_", tm, ".pdf", sep = "")
figS8_jpg <- str_c("figS8_", "EPICON_com_", tm, ".jpg", sep = "")

ggsave(figS8_pdf, figS8_EPICON_com, width = 5.62, height = 5.42)
ggsave(figS8_jpg, figS8_EPICON_com, width = 5.62, height = 5.42)


# envs
EPICON_rDNAm_env1 <- 
  EPICON_rDNAm_env %>% 
  mutate(grp = str_c(Habitat, Treatment, sep = "_"))
# EPICON_rDNAm_env1$Timepiont

EPICON_rDNAm_env1$Habitat %>% table()
EPICON_rDNAm_env1$Treatment %>% table()
EPICON_rDNAm_env1$grp %>% table()

EPICON_rDNAm_plot <- 
  ggplot(EPICON_rDNAm_env1) +
  geom_point(
             aes(Timepiont, rDNAm, colour = grp),
             position = position_dodge(width = 0.7), size = 1) +
  geom_boxplot(
               aes(Timepiont, rDNAm, colour = grp),
               position = position_dodge(width = 0.7), width = 0.4) +
  geom_smooth(aes(timepoint + 1, rDNAm, colour = grp), se = F) +
  geom_vline(xintercept = seq(1.5, 17.5, 1), linewidth = 0.05, linetype = "dashed") +
  #annotate(geom = "text", x = 9.5, y = 96, label = "Habitat(H): df = 3, p = 2.52e-89") +
  #annotate(geom = "text", x = 9.5, y = 94, label = "Treatment(T): df = 2, p = 1.11e-8") +
  #annotate(geom = "text", x = 9.5, y = 92, label = "H * T: df = 6, p = 4.61e-42") +
  scale_colour_manual(values = c("lightgreen", "green", "darkgreen",
                                 "pink", "red", "darkred",
                                 "lightblue", "blue", "darkblue",
                                 "grey80", "grey40", "black")) +
  guides(colour = guide_legend(
    title = NULL,
    ncol = 4)) +
  labs(x = "Time", y = "rDNA copy number") + 
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, 0.9),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 20, colour = "black", face = "bold"),
        legend.title = element_text(face = "bold", size = 18),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black", linewidth = 1.5))
EPICON_rDNAm_plot 

# ggsave("EPICON_rDNAm_plot_20241211_new.jpg", EPICON_rDNAm_plot)

# facet ------------
EPICON_rDNAm_env1 <- 
  EPICON_rDNAm_env1 %>% mutate(
    Treatment = case_when(
      Treatment == "Control" ~ "Control",
      Treatment == "Post_flowering_drought" ~ "Post-flowering drought",
      Treatment == "Pre_flowering_drought" ~ "Pre-flowering drought"
    )
  )

EPICON_rDNAm_env1$Treatment <- factor(EPICON_rDNAm_env1$Treatment, levels = c("Control", "Pre-flowering drought", "Post-flowering drought"))

EPICON_rDNAm_env1$Timepiont

EPICON_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_rDNAm_env1)
EPICON_rDNAm_aov_sum <- summary(EPICON_rDNAm_aov)

EPICON_rDNAm_aov %>% tidy()

EPICON_rDNAm_plot_facet_Treatment <- 
  ggplot(EPICON_rDNAm_env1) +
  geom_point(
    aes(Timepiont, rDNAm, colour = Treatment),
    position = position_dodge(width = 0.8), size = 2) +
  geom_boxplot(
    aes(Timepiont, rDNAm, colour = Treatment),
    position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_rDNAm_env1 %>% filter(Timepiont != "TP00"),
              aes(timepoint + 1, rDNAm, colour = Treatment), se = F) +
  facet_wrap(~ Habitat, ncol = 1) +
  scale_x_discrete(labels = seq(0, 17, 1)) +
  scale_colour_manual(values = c("black",
                                 "blue",
                                 "red")) +
  guides(colour = guide_legend(
    title = "Treatment",
    ncol = 1)) +
  labs(x = "Week",
       y = "Community-weighted rDNA copy number",
       subtitle = "Compartment (C): df = 3, p = 1.64e-208; Treatment (T): df = 2, p = 4.08e-19; Week (W): df = 17, p = 7.04e-53\n
       C * T: df = 6, p = 8.21e-67; C * W: df = 48, p = 8.94e-85; T * W: df = 25, p = 0.0181\n
       C * T * W: df = 69, p = 2.76e-7") + 
  theme_bw() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.91),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        strip.text = element_text(size = 12),
        plot.subtitle = element_text(face = "bold", vjust = 0.5, hjust = 0.5))
EPICON_rDNAm_plot_facet_Treatment


EPICON_rDNAm_env1 <- 
  EPICON_rDNAm_env1 %>% mutate(
    show_type = if_else(timepoint == 0, 1, 0)
  )



# fig3c, EPICON ------------------------------
EPICON_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_rDNAm_env1)
EPICON_rDNAm_aov_sum <- summary(EPICON_rDNAm_aov)

EPICON_rDNAm_aov %>% tidy()

fig3c_EPICON_rDNAm_plot_Treatment_Compartment <- 
  ggplot(EPICON_rDNAm_env1) +
  geom_point(
    aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_rDNAm_env1 %>% filter(Timepiont != "TP00"),
              aes(timepoint + 1, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 9.5, y = 85, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "1.31e-211;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "6.63e-30;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "6.19e-49")) +
  annotate(geom = "text", x = 9.5, y = 81, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "5.86e-82;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "1.72e-101;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "0.0107")) +
  annotate(geom = "text", x = 9.5, y = 77, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "5.86e-8")) +
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
fig3c_EPICON_rDNAm_plot_Treatment_Compartment

tm <- now() %>% str_split_i(pattern = " ", 1)
fig3c_pdf <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".pdf", sep = "")
fig3c_jpg <- str_c("fig3c_", "EPICON_rDNAm_plot_Treatment_Compartment_", tm, ".jpg", sep = "")


ggsave(fig3c_pdf, fig3c_EPICON_rDNAm_plot_Treatment_Compartment, width = 11.5, height = 5.8)
ggsave(fig3c_jpg, fig3c_EPICON_rDNAm_plot_Treatment_Compartment, width = 11.5, height = 5.8)



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


EPICON_fungal_otutab1_1 <- EPICON_fungal_otutab1 %>% mutate(
  OTU_ID = str_split_i(rownames(EPICON_fungal_otutab1), i = 1, pattern = "_"), .before = TP01L01
)

FGanno_EPICON_FFF1$Guild2 %>% table()

# EPICON PP otutab -------------
# 163 PP OTUs
EPICON_fungal_otutab_PP <- EPICON_fungal_otutab1_1 %>% 
  left_join(FGanno_EPICON_FFF1 %>% select(ID, Guild2), by = c("OTU_ID" = "ID")) %>%
  filter(Guild2 == "Plant pathogen fungi") %>% select(-Guild2)

EPICON_fungal_otutab_PP1 <- EPICON_fungal_otutab_PP %>% 
  mutate(sum_abd = apply(EPICON_fungal_otutab_PP %>% select(-OTU_ID), sum, MARGIN = 1)) %>% 
  filter(sum_abd != 0) %>% select(-sum_abd)

EPICON_fungal_otutab_PP1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% hist()
EPICON_fungal_otutab_PP1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% min()


EPICON_fungal_otutab_PP_r <- rrarefy(EPICON_fungal_otutab_PP1 %>% select(-OTU_ID) %>% t() , 105) %>% t() %>% as.data.frame()

EPICON_fungal_otutab_PP_r %>% apply(sum, MARGIN = 2)


EPICON_fungal_otutab_PP_per <- apply(EPICON_fungal_otutab_PP_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = EPICON_fungal_otutab_PP1$OTU_ID, .before = TP01L01)

EPICON_fungal_otutab_PP_per1 <- 
  EPICON_fungal_otutab_PP_per %>% 
  left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))

sample_id <- colnames(EPICON_fungal_otutab_PP_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab_PP_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
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

EPICON_fungal_otutab_PP_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
EPICON_fungal_otutab_PP_rDNAm


EPICON_PP_rDNAm <- 
  EPICON_fungal_otutab_PP_rDNAm %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")


EPICON_PP_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_PP_rDNAm)
EPICON_PP_rDNAm_aov_sum <- summary(EPICON_PP_rDNAm_aov)

EPICON_PP_rDNAm_aov %>% tidy()

# figS10, plant pathogen --------------------------------
figS10_EPICON_PP <- 
  ggplot(EPICON_PP_rDNAm) +
  #geom_point(
  #  aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_PP_rDNAm,
              aes(timepoint, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 8.5, y = 90, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "1.47e-134;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "1.03e-27;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "1.77e-23")) +
  annotate(geom = "text", x = 8.5, y = 87, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "7.53e-53;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "1.39e-53;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "7.07e-6")) +
  annotate(geom = "text", x = 8.5, y = 84, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "4.69e-3")) +
  #annotate(geom = "text", x = 9, y = 100,
  #         label = "Compartment (C): df = 3, p = 1.46e-164; Treatment (T): df = 2, p = 3.19e-19; Week (W): df = 17, p = 1.01e-77\n
  #         C * T: df = 6, p = 8.23e-90; C * W: df = 48, p = 9.80e-95; T * W: df = 25, p = 2.72e-6\n
  #         C * T * W: df = 69, p = 7.71e-7", hjust = 0.5, vjust = 0.5,
  #         fontface = "bold",
  #         size = 4) +
  # facet_wrap(~ Treatment, ncol = 1) +
  scale_x_continuous(labels = seq(0, 17, 1), breaks = seq(0, 17, 1), limits = c(0, 17)) +
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
figS10_EPICON_PP

tm <- now() %>% str_split_i(pattern = " ", 1)
figS10_pdf <- str_c("figS10_", "EPICON_PP_", tm, ".pdf", sep = "")
figS10_jpg <- str_c("figS10_", "EPICON_PP_", tm, ".jpg", sep = "")

ggsave(figS10_pdf, figS10_EPICON_PP, width = 11.5, height = 5.8)
ggsave(figS10_jpg, figS10_EPICON_PP, width = 11.5, height = 5.8)


# EPICON Sap otutab ------
# 544 Sap OTUs
EPICON_fungal_otutab_Sap <- EPICON_fungal_otutab1_1 %>% 
  left_join(FGanno_EPICON_FFF1 %>% select(ID, Guild2), by = c("OTU_ID" = "ID")) %>%
  filter(Guild2 == "Saprotroph fungi") %>% select(-Guild2)

EPICON_fungal_otutab_Sap1 <- EPICON_fungal_otutab_Sap %>% 
  mutate(sum_abd = apply(EPICON_fungal_otutab_Sap %>% select(-OTU_ID), sum, MARGIN = 1)) %>% 
  filter(sum_abd != 0) %>% select(-sum_abd)

#EPICON_fungal_otutab_Sap1 %>% select(-OTU_ID) %>% apply(sum, MARGIN = 2) %>% as.data.frame() %>% view()


EPICON_fungal_otutab_Sap_r <- rrarefy(EPICON_fungal_otutab_Sap1 %>% select(-OTU_ID, -TP05L01) %>% t() , 214) %>% t() %>% as.data.frame()

EPICON_fungal_otutab_Sap_r %>% apply(sum, MARGIN = 2)


EPICON_fungal_otutab_Sap_per <- apply(EPICON_fungal_otutab_Sap_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = EPICON_fungal_otutab_Sap1$OTU_ID, .before = TP01L01)

EPICON_fungal_otutab_Sap_per1 <- 
  EPICON_fungal_otutab_Sap_per %>% 
  left_join(EPICON_frDNA_anno %>% select(ID, rDNA_m), by = c("OTU_ID" = "ID"))

sample_id <- colnames(EPICON_fungal_otutab_Sap_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id


rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- EPICON_fungal_otutab_Sap_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
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

EPICON_fungal_otutab_Sap_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
EPICON_fungal_otutab_Sap_rDNAm

EPICON_Sap_rDNAm <- 
  EPICON_fungal_otutab_Sap_rDNAm %>% 
  left_join(env %>% select(aa, Timepiont, Habitat, Treatment), by = c("sample_id" = "aa")) %>% 
  left_join(env_timepoint, by = "Timepiont")


EPICON_Sap_rDNAm_aov <- aov(rDNAm ~ Habitat * Treatment * Timepiont, data = EPICON_Sap_rDNAm)
EPICON_Sap_rDNAm_aov_sum <- summary(EPICON_Sap_rDNAm_aov)

EPICON_Sap_rDNAm_aov %>% tidy()

# figS9, saprotroph -------------------------------
figS9_EPICON_Sap <- 
  ggplot(EPICON_Sap_rDNAm) +
  #geom_point(
  #  aes(Timepiont, rDNAm, colour = Habitat), alpha = EPICON_rDNAm_env1$show_type, size = 1.5) +
  #geom_boxplot(
  #  aes(Timepiont, rDNAm, colour = Habitat),
  #  position = position_dodge(width = 0.8), width = 0.5, outliers = T, outlier.colour = "grey") +
  geom_smooth(data = EPICON_Sap_rDNAm,
              aes(timepoint, rDNAm, colour = Habitat, linetype = Treatment), se = T, alpha = 0.05) +
  annotate(geom = "text", x = 8.5, y = 83, label = expression(
    "Compartment (C):" ~~ df == "3," ~~ italic(p) == "5.71e-167;" ~~ "Treatment (T):" ~~ df == "2," ~~ italic(p) == "3.04e-7;" ~~ "Week (W):" ~~ df == "17," ~~ italic(p) == "2.31e-46")) +
  annotate(geom = "text", x = 8.5, y = 80, label = expression(
    "C x T:" ~~ df == "6," ~~ italic(p) == "3.40e-22;" ~~ "C x W:" ~~ df == "48," ~~ italic(p) == "6.08e-64;" ~~ "T x W:" ~~ df == "25," ~~ italic(p) == "0.155")) +
  annotate(geom = "text", x = 8.5, y = 77, label = expression(
    "C x T x W:" ~~ df == "69," ~~ italic(p) == "6.85e-12")) +
  #annotate(geom = "text", x = 9, y = 100,
  #         label = "Compartment (C): df = 3, p = 1.46e-164; Treatment (T): df = 2, p = 3.19e-19; Week (W): df = 17, p = 1.01e-77\n
  #         C * T: df = 6, p = 8.23e-90; C * W: df = 48, p = 9.80e-95; T * W: df = 25, p = 2.72e-6\n
  #         C * T * W: df = 69, p = 7.71e-7", hjust = 0.5, vjust = 0.5,
  #         fontface = "bold",
  #         size = 4) +
  # facet_wrap(~ Treatment, ncol = 1) +
  scale_x_continuous(labels = seq(0, 17, 1), breaks = seq(0, 17, 1), limits = c(0, 17)) +
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
figS9_EPICON_Sap


tm <- now() %>% str_split_i(pattern = " ", 1)
figS9_pdf <- str_c("figS9_", "EPICON_Sap_", tm, ".pdf", sep = "")
figS9_jpg <- str_c("figS9_", "EPICON_Sap_", tm, ".jpg", sep = "")

ggsave(figS9_pdf, figS9_EPICON_Sap, width = 11.5, height = 5.8)
ggsave(figS9_jpg, figS9_EPICON_Sap, width = 11.5, height = 5.8)

# XC-Li, leat and root fungi -----------------------

# meta data
XC_metadata <- read.csv("./database/fun.env_510.csv", header = T)

# root
XC_root_otutab <- read.csv("./database/root.fun.otu_510.csv", header = T)
XC_root_otutab_tax <- read.csv("./database/root.fun.tax_510.csv", header = T)

# leaf
XC_leaf_otutab <- read.csv("./database/leaf.fun.otu_510.csv", header = T)
XC_leaf_otutab_tax <- data.frame(
  OTU = XC_leaf_otutab$OTU.ID,
  taxonomy = XC_leaf_otutab$X
)

XC_leaf_otu_ind <- data.frame(
  OTU_ID = str_c("OTU", 1:15742),
  OTU.ID = XC_leaf_otutab$OTU.ID
)

XC_leaf_otutab1 <- XC_leaf_otu_ind %>% left_join(XC_leaf_otutab, by = "OTU.ID") %>% 
  select(-OTU.ID) %>% 
  select(-(2:9))

XC_leaf_otutab_tax <- 
  XC_leaf_otutab_tax %>% 
  mutate(
    kingdom = str_split_i(taxonomy, pattern = ",", 1)
  ) %>% 
  mutate(
    kingdom = str_extract(kingdom, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    phylum = str_split_i(taxonomy, pattern = ",", 2)
  ) %>% 
  mutate(
    phylum = str_extract(phylum, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    class = str_split_i(taxonomy, pattern = ",", 3)
  ) %>% 
  mutate(
    class = str_extract(class, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    order = str_split_i(taxonomy, pattern = ",", 4)
  ) %>% 
  mutate(
    order = str_extract(order, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    family = str_split_i(taxonomy, pattern = ",", 5)
  ) %>% 
  mutate(
    family = str_extract(family, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    genus = str_split_i(taxonomy, pattern = ",", 6)
  ) %>% 
  mutate(
    genus = str_extract(genus, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    species = str_split_i(taxonomy, pattern = ",", 7)
  ) %>% 
  mutate(
    species = str_extract(species, pattern = "[^\\(]+") %>% str_sub(start = 3)
  ) %>% 
  mutate(
    species = str_c(str_split_i(species, pattern = "_", 1), str_split_i(species, pattern = "_", 2), sep = " ")
  )


XC_leaf_otutab_tax1 <- XC_leaf_otu_ind %>% left_join(XC_leaf_otutab_tax, by = c("OTU.ID" = "OTU")) %>% select(-OTU.ID)

# data we needed
XC_leaf_otutab1
XC_leaf_otutab_tax1


# rCNV tab
rCNV_phy_tab
rCNV_cla_tab
rCNV_ord_tab
rCNV_fam_tab
rCNV_gen_tab
rCNV_spc_tab

# spc
rCNV_XCleaf_spc <- XC_leaf_otutab_tax1 %>% left_join(rCNV_spc_tab, by = c("species" = "spc")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_spc

rCNV_XCleaf_spc_na <- XC_leaf_otutab_tax1 %>% left_join(rCNV_spc_tab, by = c("species" = "spc")) %>% 
  filter(is.na(rDNA_m))

# gen
rCNV_XCleaf_gen <- rCNV_XCleaf_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("genus" = "Gen")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_gen

rCNV_XCleaf_gen_na <- rCNV_XCleaf_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("genus" = "Gen")) %>% 
  filter(is.na(rDNA_m))

# fam
rCNV_XCleaf_fam <- rCNV_XCleaf_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("family" = "fam")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_fam

rCNV_XCleaf_fam_na <- rCNV_XCleaf_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("family" = "fam")) %>% 
  filter(is.na(rDNA_m))


# ord
rCNV_XCleaf_ord <- rCNV_XCleaf_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("order" = "ord")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_ord

rCNV_XCleaf_ord_na <- rCNV_XCleaf_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("order" = "ord")) %>% 
  filter(is.na(rDNA_m))


# cla
rCNV_XCleaf_cla <- rCNV_XCleaf_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("class" = "cla")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_cla

rCNV_XCleaf_cla_na <- rCNV_XCleaf_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("class" = "cla")) %>% 
  filter(is.na(rDNA_m))


# phy
rCNV_XCleaf_phy <- rCNV_XCleaf_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("phylum" = "phy")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCleaf_phy

rCNV_XCleaf_phy_na <- rCNV_XCleaf_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("phylum" = "phy")) %>% 
  filter(is.na(rDNA_m))
# rCNV_XCleaf_phy_na %>% dim()
# rCNV_XCleaf_phy_na %>% view()

# rCNV XCleaf anno --------------------------
rCNV_XCleaf_anno <- rbind(rCNV_XCleaf_phy, rCNV_XCleaf_cla, rCNV_XCleaf_ord, rCNV_XCleaf_fam, rCNV_XCleaf_gen, rCNV_XCleaf_spc)
# dim(rCNV_XCleaf_anno)
# done

# XC_leaf_otutab1 %>% view()
XCleaf_otutab1_per <- apply(XC_leaf_otutab1 %>% select(-1), rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = XC_leaf_otutab1$OTU_ID, .before = X1)


XCleaf_otutab1_per1 <- 
  XCleaf_otutab1_per %>% 
  left_join(rCNV_XCleaf_anno %>% select(OTU_ID, rDNA_m), by = "OTU_ID")

sample_id <- colnames(XCleaf_otutab1_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id

rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- XCleaf_otutab1_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
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

XCleaf_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
XCleaf_rDNAm

XCleaf_rDNAm1 <- 
  XCleaf_rDNAm %>% mutate(
    grp = if_else(str_detect(sample_id, "X"), "BS", "NS")
  )



XCleaf_kru <- kruskal(XCleaf_rDNAm1$rDNAm, XCleaf_rDNAm1$grp)
XCleaf_kru

# view(XC_metadata)

XCleaf_rDNAm2 <- XCleaf_rDNAm1 %>% left_join(XC_metadata, by = c("sample_id" = "Sample.No."))

# XCleaf_rDNAm2 %>% select(grp, group)

XCleaf_rDNAm_lat <- 
  ggplot(XCleaf_rDNAm2, aes(Latitude, rDNAm)) +
  geom_point(size = 1, colour = "blue", alpha = 0.3) +
  geom_smooth(colour = "red") +
  scale_colour_npg() +
  facet_wrap(~ group) +
  labs(y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.title = element_text(colour = "black", face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 12),
    strip.text = element_text(size = 14)
  )
XCleaf_rDNAm_lat


XCleaf_rDNAm_Tem <- 
  ggplot(XCleaf_rDNAm2, aes(Temperature.seasonality, rDNAm)) +
  geom_point(size = 1, colour = "blue", alpha = 0.3) +
  geom_smooth(colour = "red") +
  scale_colour_npg() +
  facet_wrap(~ group) +
  labs(y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.title = element_text(colour = "black", face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 12),
    strip.text = element_text(size = 14)
  )
XCleaf_rDNAm_Tem

# root -----------
XC_metadata
XC_root_otutab
XC_root_otutab_tax

XC_root_otutab_tax1 <- XC_root_otutab_tax %>% select(OTU, kingdom, phylum, class, order, family, genus, species) %>% 
  mutate(species = str_sub(species, start = 4)) %>% 
  mutate(species = str_replace(species, pattern = "\\.", " "))

XC_root_otutab_tax1


# rCNV tab
rCNV_phy_tab
rCNV_cla_tab
rCNV_ord_tab
rCNV_fam_tab
rCNV_gen_tab
rCNV_spc_tab

# spc
rCNV_XCroot_spc <- XC_root_otutab_tax1 %>% left_join(rCNV_spc_tab, by = c("species" = "spc")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_spc

rCNV_XCroot_spc_na <- XC_root_otutab_tax1 %>% left_join(rCNV_spc_tab, by = c("species" = "spc")) %>% 
  filter(is.na(rDNA_m))

# gen
rCNV_XCroot_gen <- rCNV_XCroot_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("genus" = "Gen")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_gen

rCNV_XCroot_gen_na <- rCNV_XCroot_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("genus" = "Gen")) %>% 
  filter(is.na(rDNA_m))

# fam
rCNV_XCroot_fam <- rCNV_XCroot_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("family" = "fam")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_fam

rCNV_XCroot_fam_na <- rCNV_XCroot_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("family" = "fam")) %>% 
  filter(is.na(rDNA_m))


# ord
rCNV_XCroot_ord <- rCNV_XCroot_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("order" = "ord")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_ord

rCNV_XCroot_ord_na <- rCNV_XCroot_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("order" = "ord")) %>% 
  filter(is.na(rDNA_m))


# cla
rCNV_XCroot_cla <- rCNV_XCroot_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("class" = "cla")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_cla

rCNV_XCroot_cla_na <- rCNV_XCroot_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("class" = "cla")) %>% 
  filter(is.na(rDNA_m))


# phy
rCNV_XCroot_phy <- rCNV_XCroot_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("phylum" = "phy")) %>% 
  filter(!is.na(rDNA_m))
rCNV_XCroot_phy

rCNV_XCroot_phy_na <- rCNV_XCroot_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("phylum" = "phy")) %>% 
  filter(is.na(rDNA_m))
# rCNV_XCroot_phy_na %>% dim()
# rCNV_XCroot_phy_na %>% view()

# rCNV XCroot anno --------------------------
rCNV_XCroot_anno <- rbind(rCNV_XCroot_phy, rCNV_XCroot_cla, rCNV_XCroot_ord, rCNV_XCroot_fam, rCNV_XCroot_gen, rCNV_XCroot_spc)
# dim(rCNV_XCroot_anno)
# done

colnames(XC_root_otutab)

XCroot_otutab1_per <- apply(XC_root_otutab %>% select(-1), rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = XC_root_otutab$OTU, .before = bet001)


XCroot_otutab1_per1 <- 
  XCroot_otutab1_per %>% 
  left_join(rCNV_XCroot_anno %>% select(OTU, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(XCroot_otutab1_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id

rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- XCroot_otutab1_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
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

XCroot_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
XCroot_rDNAm

XCroot_rDNAm1 <- XCroot_rDNAm %>% left_join(XC_metadata, by = c("sample_id" = "Sample.No."))

#XCroot_rDNAm1 %>% select(group)

XCroot_rDNAm_Tem <- 
  ggplot(XCroot_rDNAm1, aes(Temperature.seasonality, rDNAm)) +
  geom_point(size = 1, colour = "blue", alpha = 0.3) +
  geom_smooth(colour = "red") +
  scale_colour_npg() +
  facet_wrap(~ group) +
  labs(y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.title = element_text(colour = "black", face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 12),
    strip.text = element_text(size = 14)
  )
XCroot_rDNAm_Tem


XCroot_rDNAm_lm_rlt <- lm_rlt_Spcs(df = XCroot_rDNAm1, yval = "rDNAm", xval = "Latitude", subgrp = "all")
XCroot_rDNAm_lm_rlt

XCroot_rDNAm_lm <- lm(rDNAm ~ Latitude, data = XCroot_rDNAm1)


library(nlme)
library(lme4)
library(lmerTest)

# lme4 packages results
XCroot_rDNAm_lat_lmer <- lmerTest::lmer(rDNAm ~ Latitude + (1 | Family), data = XCroot_rDNAm1)
XCroot_rDNAm_lat_lmer_sum <- summary(XCroot_rDNAm_lat_lmer)


# Prof. Gao codes
XCroot_rDNAm_lat_lme <- lme(rDNAm ~ Latitude, random = ~1 | Family, data = XCroot_rDNAm1)
XCroot_rDNAm_lat_lme_sum <- summary(XCroot_rDNAm_lat_lme)
XCroot_rDNAm_lat_lme_sum
anova(XCroot_rDNAm_lat_lme)



library(MuMIn)
r.squaredGLMM(XCroot_rDNAm_lat_lme)[2]
r.squaredGLMM(XCroot_rDNAm_lat_lmer)[2]

XCroot_rDNAm_lat <- 
  ggplot(XCroot_rDNAm1, aes(Latitude, rDNAm, colour = Family)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text(data = XCroot_rDNAm_lm_rlt, aes(x = anno_x1, y = anno_y1, label = rlt_sig), parse = T, colour = "red", size = 5) +
  geom_smooth(colour = "blue", method = "lm") +
  scale_colour_nejm() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.title = element_text(colour = "black", face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 12),
    strip.text = element_text(size = 14),
    plot.subtitle = element_text(face = "bold", hjust = 0.5, vjust = 0.5),
    legend.title = element_text(face = "bold")
  )
XCroot_rDNAm_lat



# root EMF ----------
# XC_root_otutab_tax$EMF %>% table()

XC_root_otutab_EMFind <- XC_root_otutab_tax %>% filter(EMF == "emf") %>% select(OTU) %>% pull()
XC_root_otutab_EMF <- XC_root_otutab %>% filter(OTU %in% XC_root_otutab_EMFind)

EMF_sample_sub <- apply(XC_root_otutab_EMF %>% select(-1), sum, MARGIN = 2) %>%
  as.data.frame() 
  
colnames(EMF_sample_sub) <- "Freq"
EMF_sample_sub$Freq %>% hist()
EMF_sample_sub1 <- EMF_sample_sub %>% mutate(sample_id = rownames(EMF_sample_sub), .before = Freq) %>% 
  filter(Freq > 300)


XC_root_otutab_EMF1 <- XC_root_otutab_EMF %>% select(OTU, all_of(EMF_sample_sub1$sample_id))

apply(XC_root_otutab_EMF1 %>% select(-OTU), sum, MARGIN = 2) %>% min()

XC_root_otutab_EMF1_r <- rrarefy(XC_root_otutab_EMF1 %>% select(-OTU) %>% t(), 301) %>% t() %>% as.data.frame()

XC_root_otutab_EMF1_r <- XC_root_otutab_EMF1_r %>% mutate(OTU = XC_root_otutab_EMF1$OTU, .before = bet001)


XCroot_otutab_EMF1_per <- apply(XC_root_otutab_EMF1_r %>% select(-OTU), rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = XC_root_otutab_EMF1_r$OTU, .before = bet001)

XCroot_otutab_EMF1_per1 <- 
  XCroot_otutab_EMF1_per %>% 
  left_join(rCNV_XCroot_anno %>% select(OTU, rDNA_m), by = c("OTU_ID" = "OTU"))

sample_id <- colnames(XCroot_otutab_EMF1_per1 %>% select(-OTU_ID, -rDNA_m))
sample_id

rDNAm_calcu <- function(sample_id) {
  
  tmp_df <- XCroot_otutab_EMF1_per1 %>% select(OTU_ID, all_of(sample_id), rDNA_m)
  
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

XCroot_EMF_rDNAm <- map_dfr(sample_id, ~ rDNAm_calcu(.), .progress = T)
XCroot_EMF_rDNAm

XCroot_EMF_rDNAm1 <- XCroot_EMF_rDNAm %>% left_join(XC_metadata, by = c("sample_id" = "Sample.No."))


XCroot_EMF_rDNAm_lat_lme <- lme(rDNAm ~ Latitude, random = ~1 | Family, data = XCroot_EMF_rDNAm1)
XCroot_EMF_rDNAm_lat_lme_sum <- summary(XCroot_EMF_rDNAm_lat_lme)
XCroot_EMF_rDNAm_lat_lme_sum
anova(XCroot_rDNAm_lat_lme)
anova(XCroot_EMF_rDNAm_lat_lme)


library(MuMIn)
r.squaredGLMM(XCroot_rDNAm_lat_lme)[2]
r.squaredGLMM(XCroot_EMF_rDNAm_lat_lme)[2]

lm(rDNAm ~ Latitude, data = XCroot_EMF_rDNAm1) %>% summary()

XCroot_EMF_rDNAm_lm_rlt <- lm_rlt_Spcs(df = XCroot_EMF_rDNAm1, yval = "rDNAm", xval = "Latitude", subgrp = "all")
XCroot_EMF_rDNAm_lm_rlt
# italic(r) == 0.198 ~~ italic(p) == 4.31e-13

#fix(XCroot_EMF_rDNAm_lm_rlt)

# figS6, Community-weighted rCNV ~ latitude of ECM ----------------------------
figS6_XCroot_EMF_rDNAm_lat <- 
  ggplot(XCroot_EMF_rDNAm1, aes(Latitude, rDNAm, colour = Family)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 0.4) +
  geom_text(data = XCroot_EMF_rDNAm_lm_rlt, aes(x = anno_x1, y = 140, label = rlt_sig), parse = T, colour = "red", size = 4.5) +
  geom_smooth(colour = "red", method = "lm") +
  scale_colour_nejm() +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal")) +
  labs(y = "Community-weighted rDNA copy number of ectomycorrhizal fungi") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.95),
    axis.title.x = element_text(colour = "black", face = "bold", size = 12),
    axis.title.y = element_text(colour = "black", face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 12),
    strip.text = element_text(size = 14),
    plot.subtitle = element_text(face = "bold", hjust = 0.5, vjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )
figS6_XCroot_EMF_rDNAm_lat

tm <- now() %>% str_split_i(pattern = " ", 1)
figS6_pdf <- str_c("figS6_", "XCroot_EMF_rDNAm_lat_", tm, ".pdf", sep = "")
figS6_jpg <- str_c("figS6_", "XCroot_EMF_rDNAm_lat_", tm, ".jpg", sep = "")

ggsave(figS6_pdf, figS6_XCroot_EMF_rDNAm_lat, width = 6.6, height = 6.43)
ggsave(figS6_jpg, figS6_XCroot_EMF_rDNAm_lat, width = 6.6, height = 6.43)

#XCroot_rDNAm_lat_p <- XCroot_rDNAm_lat + XCroot_EMF_rDNAm_lat + plot_layout(guides = "collect")
#ggsave("XCroot_rDNAm_lat_p_20241218.pdf", XCroot_rDNAm_lat_p)
#ggsave("XCroot_rDNAm_lat_p_20241218.jpg", XCroot_rDNAm_lat_p)



# Zheng data another version ----------
Z_soil_meta_cfb <- read.csv("./database/Zheng_CFB.env.2022.03.28.csv", header = T)
Z_soil_taxa_cfb <- read.csv("./database/Zheng_CFB.fung.ID.2022.03.28.csv", header = T)
Z_soil_otutab_cfb <- read.csv("./database/Zheng_CFB.fung.2022.03.28.csv", header = T)


# Z another anno
rCNV_phy_tab
rCNV_cla_tab
rCNV_ord_tab
rCNV_fam_tab
rCNV_gen_tab
rCNV_spc_tab


# taxa tab
Z_soil_taxa_cfb1 <- Z_soil_taxa_cfb %>% select(OTU, Phylum, Class, Order, Family, Genus, Species, Lifestyle2) %>% 
  mutate(Species = str_replace(Species, pattern = "_", replacement = " "))
#view(Z_soil_taxa_cfb1)


# spc
rCNV_Zsoil_cfb_spc <- Z_soil_taxa_cfb1 %>% left_join(rCNV_spc_tab, by = c("Species" = "spc")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_spc

rCNV_Zsoil_cfb_spc_na <- Z_soil_taxa_cfb1 %>% left_join(rCNV_spc_tab, by = c("Species" = "spc")) %>% 
  filter(is.na(rDNA_m))


# gen
rCNV_Zsoil_cfb_gen <- rCNV_Zsoil_cfb_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("Genus" = "Gen")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_gen

rCNV_Zsoil_cfb_gen_na <- rCNV_Zsoil_cfb_spc_na %>% select(-rDNA_m) %>% left_join(rCNV_gen_tab, by = c("Genus" = "Gen")) %>% 
  filter(is.na(rDNA_m))

# fam
rCNV_Zsoil_cfb_fam <- rCNV_Zsoil_cfb_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("Family" = "fam")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_fam

rCNV_Zsoil_cfb_fam_na <- rCNV_Zsoil_cfb_gen_na %>% select(-rDNA_m) %>% left_join(rCNV_fam_tab, by = c("Family" = "fam")) %>% 
  filter(is.na(rDNA_m))


# ord
rCNV_Zsoil_cfb_ord <- rCNV_Zsoil_cfb_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("Order" = "ord")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_ord

rCNV_Zsoil_cfb_ord_na <- rCNV_Zsoil_cfb_fam_na %>% select(-rDNA_m) %>% left_join(rCNV_ord_tab, by = c("Order" = "ord")) %>% 
  filter(is.na(rDNA_m))


# cla
rCNV_Zsoil_cfb_cla <- rCNV_Zsoil_cfb_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("Class" = "cla")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_cla

rCNV_Zsoil_cfb_cla_na <- rCNV_Zsoil_cfb_ord_na %>% select(-rDNA_m) %>% left_join(rCNV_cla_tab, by = c("Class" = "cla")) %>% 
  filter(is.na(rDNA_m))


# phy
rCNV_Zsoil_cfb_phy <- rCNV_Zsoil_cfb_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(!is.na(rDNA_m))
rCNV_Zsoil_cfb_phy

rCNV_Zsoil_cfb_phy_na <- rCNV_Zsoil_cfb_cla_na %>% select(-rDNA_m) %>% left_join(rCNV_phy_tab, by = c("Phylum" = "phy")) %>% 
  filter(is.na(rDNA_m))
#rCNV_Zsoil_cfb_phy_na %>% dim()
#rCNV_Zsoil_cfb_phy_na %>% view()

# rCNV Zsoil anno --------------------------
rCNV_Zsoil_cfb_anno <- rbind(rCNV_Zsoil_cfb_phy, rCNV_Zsoil_cfb_cla, rCNV_Zsoil_cfb_ord, 
                             rCNV_Zsoil_cfb_fam, rCNV_Zsoil_cfb_gen, rCNV_Zsoil_cfb_spc)
dim(rCNV_Zsoil_cfb_anno)
# done


# calcu ---------- 
colnames(rCNV_Zsoil_cfb_anno)

Z_soil_otutab_cfb %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min()
# min 5276

Z_soil_otutab_cfb_r <- rrarefy(Z_soil_otutab_cfb %>% select(-X) %>% t() , 5276) %>% t() %>% as.data.frame()
Z_soil_otutab_cfb_r %>% apply(sum, MARGIN = 2)


colnames(Z_soil_otutab_cfb_r) 

Zsoil_otutab_cfb_per <- apply(Z_soil_otutab_cfb_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb$X, .before = X301)

Zsoil_otutab_cfb_per1 <- 
  Zsoil_otutab_cfb_per %>% 
  left_join(rCNV_Zsoil_cfb_anno %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

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


lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Factors == subgrp)
    
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
      Factors = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  } else {
    
    df_rlt <- data.frame(
      Factors = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


Zheng_soil_rDNAm1_subdata1 <- 
  map_dfr(.x = c("Latitude1", "Longitude1", "Altitude1",
  "TN", "TC", "ACa", "AFe", "AK", "AMg", "TP", "C.N", "C.P", "N.P", "soil_D"),
          .f = ~ lm_rlt_Spcs(df = Zsoil_cfb_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = ., trans = "no_trans"))
Zheng_soil_rDNAm1_subdata1

lm_rlt_Spcs(df = Zsoil_cfb_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = "all", trans = "no_trans")

Zsoil_cfb_rDNAm2$Factors

ggplot(Zsoil_cfb_rDNAm2, aes(Value, rDNAm, colour = )) +
  geom_point() +
  geom_text(data = Zheng_soil_rDNAm1_subdata1, aes(x = anno_x1, y = anno_y1, label = rlt_sig), parse = T) +
  geom_smooth(method = "lm") +
  facet_wrap(~ Factors, scales = "free") +
  theme_bw()


rCNV_Zsoil_cfb_anno_ECM <- rCNV_Zsoil_cfb_anno %>% filter(Lifestyle2 == "EcM")
Z_soil_otutab_cfb_ECM <- Z_soil_otutab_cfb %>% filter(X %in% rCNV_Zsoil_cfb_anno_ECM$OTU)

rCNV_Zsoil_cfb_anno_PP <- rCNV_Zsoil_cfb_anno %>% filter(Lifestyle2 == "Plant pathogen")
Z_soil_otutab_cfb_PP <- Z_soil_otutab_cfb %>% filter(X %in% rCNV_Zsoil_cfb_anno_PP$OTU)

rCNV_Zsoil_cfb_anno_Sap <- rCNV_Zsoil_cfb_anno %>% filter(Lifestyle2 == "Saprotroph")
Z_soil_otutab_cfb_Sap <- Z_soil_otutab_cfb %>% filter(X %in% rCNV_Zsoil_cfb_anno_Sap$OTU)


Z_soil_otutab_cfb_ECM %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min() # ok
Z_soil_otutab_cfb_ECM$X %>% length()

# ECM calcu -----------------
Z_soil_otutab_cfb_ECM_r <- rrarefy(Z_soil_otutab_cfb_ECM %>% select(-X) %>% t() , 257) %>% t() %>% as.data.frame()
Z_soil_otutab_cfb_ECM_r %>% apply(sum, MARGIN = 2)

colnames(Z_soil_otutab_cfb_ECM_r) 

Z_soil_otutab_cfb_ECM_per <- apply(Z_soil_otutab_cfb_ECM_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb_ECM$X, .before = X301)

Z_soil_otutab_cfb_ECM_per1 <- 
  Z_soil_otutab_cfb_ECM_per %>% 
  left_join(rCNV_Zsoil_cfb_anno_ECM %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

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



lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Factors == subgrp)
    
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
      Factors = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  } else {
    
    df_rlt <- data.frame(
      Factors = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


Zheng_soil_ECM_rDNAm1_subdata1 <- 
  map_dfr(.x = c("Latitude1", "Longitude1", "Altitude1",
                 "TN", "TC", "ACa", "AFe", "AK", "AMg", "TP", "C.N", "C.P", "N.P", "soil_D"),
          .f = ~ lm_rlt_Spcs(df = Zsoil_cfb_ECM_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = ., trans = "no_trans"))
Zheng_soil_ECM_rDNAm1_subdata1

lm_rlt_Spcs(df = Zsoil_cfb_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = "all", trans = "no_trans")

Zsoil_cfb_rDNAm2$Factors

ggplot(Zsoil_cfb_ECM_rDNAm2, aes(Value, rDNAm, colour = )) +
  geom_point() +
  geom_text(data = Zheng_soil_ECM_rDNAm1_subdata1, aes(x = anno_x1, y = anno_y1, label = rlt_sig), parse = T) +
  geom_smooth(method = "lm") +
  facet_wrap(~ Factors, scales = "free") +
  theme_bw()



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
  left_join(rCNV_Zsoil_cfb_anno_PP %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

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



lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Factors == subgrp)
    
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
      Factors = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  } else {
    
    df_rlt <- data.frame(
      Factors = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


Zheng_soil_PP_rDNAm1_subdata1 <- 
  map_dfr(.x = c("Latitude1", "Longitude1", "Altitude1",
                 "TN", "TC", "ACa", "AFe", "AK", "AMg", "TP", "C.N", "C.P", "N.P", "soil_D"),
          .f = ~ lm_rlt_Spcs(df = Zsoil_cfb_PP_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = ., trans = "no_trans"))
Zheng_soil_PP_rDNAm1_subdata1



ggplot(Zsoil_cfb_PP_rDNAm2, aes(Value, rDNAm, colour = )) +
  geom_point() +
  geom_text(data = Zheng_soil_PP_rDNAm1_subdata1, aes(x = anno_x1, y = anno_y1, label = rlt_sig), parse = T) +
  geom_smooth(method = "lm") +
  facet_wrap(~ Factors, scales = "free") +
  theme_bw()



# Sap calcu -----------------
Z_soil_otutab_cfb_Sap %>% select(-X) %>% apply(sum, MARGIN = 2) %>% min()

Z_soil_otutab_cfb_Sap_r <- rrarefy(Z_soil_otutab_cfb_Sap %>% select(-X) %>% t() , 2656) %>% t() %>% as.data.frame()

Z_soil_otutab_cfb_Sap_r %>% apply(sum, MARGIN = 2)


Z_soil_otutab_cfb_Sap_per <- apply(Z_soil_otutab_cfb_Sap_r, rel_per, MARGIN = 2) %>% as.data.frame() %>% 
  mutate(OTU_ID = Z_soil_otutab_cfb_Sap$X, .before = X301)

Z_soil_otutab_cfb_Sap_per1 <- 
  Z_soil_otutab_cfb_Sap_per %>% 
  left_join(rCNV_Zsoil_cfb_anno_Sap %>% select(OTU, Lifestyle2, rDNA_m), by = c("OTU_ID" = "OTU"))

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



lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all"), trans = c("no_trans")) {
  
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Factors == subgrp)
    
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
      Factors = subgrp,
      anno_x1 = (range(log10(df_tmp[[xval]]), na.rm = T)[2] - range(log10(df_tmp[[xval]]), na.rm = T)[1]) * 0.5 + range(log10(df_tmp[[xval]]), na.rm = T)[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  } else {
    
    df_rlt <- data.frame(
      Factors = subgrp,
      anno_x1 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.5 + range(df_tmp[[xval]])[1],
      anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 1.2 +
        range(df_tmp[[yval]], na.rm = T)[1],
      rsig = str_c("italic(r) == ", lmR_r),
      psig = Rp_sig,
      anno_x2 = (range(df_tmp[[xval]], na.rm = T)[2] - range(df_tmp[[xval]], na.rm = T)[1]) * 0.6 + range(df_tmp[[xval]])[1],
      anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
        range(df_tmp[[yval]], na.rm = T)[1],
      pval = lmR_p
    )
    
  }
  
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  
  return(df_rlt_F)
  
}


Zheng_soil_Sap_rDNAm1_subdata1 <- 
  map_dfr(.x = c("Latitude1", "Longitude1", "Altitude1",
                 "TN", "TC", "ACa", "AFe", "AK", "AMg", "TP", "C.N", "C.P", "N.P", "soil_D"),
          .f = ~ lm_rlt_Spcs(df = Zsoil_cfb_Sap_rDNAm2, yval = "rDNAm", xval = "Value", subgrp = ., trans = "no_trans"))
Zheng_soil_Sap_rDNAm1_subdata1



ggplot(Zsoil_cfb_Sap_rDNAm2, aes(Value, rDNAm, colour = )) +
  geom_point() +
  geom_text(data = Zheng_soil_Sap_rDNAm1_subdata1, aes(x = anno_x1, y = anno_y1, label = rlt_sig), parse = T) +
  geom_smooth(method = "lm") +
  facet_wrap(~ Factors, scales = "free") +
  theme_bw()


Zsoil_cfb_EPS_rDNAm <- rbind(Zsoil_cfb_ECM_rDNAm2 %>% mutate(Guild = "Ectomycorrhizal fungi"), 
                             Zsoil_cfb_PP_rDNAm2 %>% mutate(Guild = "Plant pathogen"), 
                             Zsoil_cfb_Sap_rDNAm2 %>% mutate(Guild = "Saprotroph fungi"))

Zheng_soil_EPS_rDNAm1_subdata <- rbind(Zheng_soil_ECM_rDNAm1_subdata1 %>% mutate(Guild = "Ectomycorrhizal fungi"),
                                       Zheng_soil_PP_rDNAm1_subdata1 %>% mutate(Guild = "Plant pathogen") %>% 
                                         mutate(anno_y1 = anno_y1 + 10),
                                       Zheng_soil_Sap_rDNAm1_subdata1 %>% mutate(Guild = "Saprotroph fungi") %>% 
                                         mutate(anno_y1 = anno_y1 - 35))




ggplot(Zsoil_cfb_EPS_rDNAm, aes(Value, rDNAm, colour = Guild)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_text(data = Zheng_soil_EPS_rDNAm1_subdata,
            aes(x = anno_x1, y = anno_y1, label = rlt_sig),
            parse = T, show.legend = F,
            fontface = "bold") +
  geom_smooth(method = "lm") +
  scale_colour_aaas() +
  facet_wrap(~ Factors, scales = "free") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black", face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.1)
  )

Zsoil_cfb_EPS_rDNAm_Lat <- 
  Zsoil_cfb_EPS_rDNAm %>% filter(Factors == "Latitude1") %>% 
  mutate(Factors = str_sub(Factors, end = -2))

Zsoil_cfb_EPS_subdata <- Zheng_soil_EPS_rDNAm1_subdata %>% filter(Factors == "Latitude1")

lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm %>% filter(Guild == "Ectomycorrhizal fungi") %>% filter(Factors == "Latitude1")) %>% 
  summary() %>% tidy()

lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm %>% filter(Guild == "Plant pathogen") %>% filter(Factors == "Latitude1")) %>% 
  summary() %>% tidy()

lm(rDNAm ~ Value,
   data = Zsoil_cfb_EPS_rDNAm %>% filter(Guild == "Saprotroph fungi") %>% filter(Factors == "Latitude1")) %>% 
  summary() %>% tidy()

#fix(Zsoil_cfb_EPS_subdata)

# figS7, soil, EPS rDNAm ~ Latitude -----------------------
figS7_Zsoil_EPS_Lat <- 
  ggplot(Zsoil_cfb_EPS_rDNAm_Lat, aes(Value, rDNAm, colour = Guild)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 0.6) +
  geom_text(data = Zsoil_cfb_EPS_subdata,
            aes(x = anno_x1, y = anno_y1,
                  label = rlt_sig),
            parse = T, show.legend = F,
            fontface = "bold",
            size = 4) +
  geom_smooth(method = "lm") +
  guides(colour = guide_legend(
    title = NULL,
    nrow = 1
  )) +
  scale_x_continuous(limits = c(20, 53)) +
  scale_y_continuous(limits = c(30, 150)) +
  scale_colour_aaas() +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", face = "bold", size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.05),
    aspect.ratio = 1
  )
figS7_Zsoil_EPS_Lat

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7_pdf <- str_c("figS7_", "Zsoil_EPS_Lat_", tm, ".pdf", sep = "")
figS7_jpg <- str_c("figS7_", "Zsoil_EPS_Lat_", tm, ".jpg", sep = "")

ggsave(figS7_pdf, figS7_Zsoil_EPS_Lat, width = 5.4, height = 5.12)
ggsave(figS7_jpg, figS7_Zsoil_EPS_Lat, width = 5.4, height = 5.12)

# combine XC and Zheng, root and soil ------------------

XCroot_rDNAm_lat

XCroot_rDNAm1_comb_tmp <- XCroot_rDNAm1 %>% select(sample_id, rDNAm, Latitude) %>% 
  mutate(Compartment = "Root")

colnames(XCroot_rDNAm1_comb_tmp)


Zsoil_cfb_rDNAm1_comb_tmp <- Zsoil_cfb_rDNAm1 %>% select(sample_id, rDNAm, Latitude1) %>% 
  mutate(Compartment = "Soil")

colnames(Zsoil_cfb_rDNAm1_comb_tmp)[3] <- "Latitude"

RS_rDNAm <- rbind(XCroot_rDNAm1_comb_tmp, Zsoil_cfb_rDNAm1_comb_tmp)

RS_rDNAm

XCroot_rDNAm_lm_rlt
Zheng_soil_rDNAm1_lat_lm_rlt <- Zheng_soil_rDNAm1_subdata1 %>% filter(Factors == "Latitude1")


colnames(XCroot_rDNAm_lm_rlt)[1] <- "Compartment"
colnames(Zheng_soil_rDNAm1_lat_lm_rlt)[1] <- "Compartment"

XCroot_rDNAm_lm_rlt$Compartment <- "Root"
Zheng_soil_rDNAm1_lat_lm_rlt$Compartment <- "Soil"

RS_rDNAm_subdata <- rbind(XCroot_rDNAm_lm_rlt, Zheng_soil_rDNAm1_lat_lm_rlt)

lm(rDNAm ~ Latitude, data = XCroot_rDNAm1) %>% tidy()
# italic(r) == 0.198 ~~ italic(p) == 1.60e-14

lm(rDNAm ~ Latitude1, data = Zsoil_cfb_rDNAm1) %>% tidy()
# italic(r) == 0.135 ~~ italic(p) == 0.037


#fix(RS_rDNAm_subdata)

# fig3b, RS rDNAm ~ Lat --------------------------
fig3b_RS_rDNAm_lat <- 
  ggplot(RS_rDNAm, aes(Latitude, rDNAm)) +
  geom_jitter(aes(colour = Compartment), size = 1.8, alpha = 0.4, width = 0.4) +
  geom_smooth(aes(colour = Compartment), method = "lm", alpha = 0.1) +
  geom_text(data = RS_rDNAm_subdata,
            aes(x = anno_x1, y = anno_y1,
                label = rlt_sig, 
                colour = Compartment),
            parse = T,
            show.legend = F,
            size = 4.5) +
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
    legend.position.inside = c(0.27, 0.9),
    aspect.ratio = 1
  )
fig3b_RS_rDNAm_lat

tm <- now() %>% str_split_i(pattern = " ", 1)
fig3b_pdf <- str_c("fig3b_", "RS_rDNAm_lat_", tm, ".pdf", sep = "")
fig3b_jpg <- str_c("fig3b_", "RS_rDNAm_lat_", tm, ".jpg", sep = "")

ggsave(fig3b_pdf, fig3b_RS_rDNAm_lat, width = 6.1, height = 5.8)
ggsave(fig3b_jpg, fig3b_RS_rDNAm_lat, width = 6.1, height = 5.8)


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


EPICON_leaf_lm <- 
  map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
          .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_env1 %>% filter(Habitat == "Leaf"),
                             yval = "rDNAm", xval = "timepoint",
                             subgrp = ., trans = "no_trans"))
EPICON_leaf_lm

EPICON_rDNAm_env1 %>% filter(Habitat == "Leaf") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_med = median(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_soil_lm <- 
  map_dfr(.x = c("Control", "Pre-flowering drought", "Post-flowering drought"),
          .f = ~ lm_rlt_Spcs(df = EPICON_rDNAm_env1 %>% filter(Habitat == "Soil"),
                             yval = "rDNAm", xval = "timepoint",
                             subgrp = ., trans = "no_trans"))
EPICON_soil_lm

EPICON_rDNAm_env1 %>% filter(Habitat == "Soil") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_med = median(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_rDNAm_env1 %>% filter(Habitat == "Root") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_med = median(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()


EPICON_rDNAm_env1 %>% filter(Habitat == "Rhizosphere") %>% group_by(Treatment, Timepiont) %>% summarise(
  rDNAm_med = median(rDNAm),
  rDNAm_sd = sd(rDNAm)
) %>% view()

rCNV_phy_tab
EPICON_rDNAm_env1 %>% dim()

# done.
