

# rCNV ~ Fungal traits ----------
# fig2d-e #
# Fungal traits database were obtained from Camenzind et al (2024) # 
# by Qiushi-Li, IM-CAS, 2025.02.20

##### packages --------------
library(tidyverse)
library(readxl)
library(ggsci)
library(vegan)
library(ggrepel)
library(ggh4x)
library(patchwork)

# RData
# load("rCNV_version_1.0_20250219.RData")

# Camenzind data
NC_func <- read_excel("./database/Camenzind_trait data_NatComm.xlsx", sheet = 2)

NC_func_F <- 
  NC_func %>% mutate(Genus = str_split_i(Species, pattern = "_", 1)) %>%
  mutate(spc = str_c(
    str_split_i(Species, pattern = "_", 1),
    str_split_i(Species, pattern = "_", 2),
    sep = " ")) %>% 
  drop_na()

# Camenzind proj
NC_proj <- read_excel("./database/NC_project.xlsx", sheet = 1)

NC_proj_list <- NC_proj %>% select(Name)
NC_proj_list

# rCNV data
rCNV_rlt_taxa_spl %>% view()


NC_rCNV_tab <- 
  rCNV_rlt_taxa_spl %>% filter(Name %in% NC_proj_list$Name) %>% 
  mutate(
    phy = str_sub(phy, start = 3),
    cla = str_sub(cla, start = 3),
    ord = str_sub(ord, start = 3),
    fam = str_sub(fam, start = 3),
    gen = str_sub(Gen, start = 3)
  ) %>% 
  mutate(spc = str_c(
    str_split_i(Name, pattern = " ", 1),
    str_split_i(Name, pattern = " ", 2),
    sep = " ")) %>% select(project, Both_ITS_LSU, gen, spc)

NC_rCNV_tab_spc <- NC_rCNV_tab %>% group_by(spc) %>% summarise(m_rDNA = mean(Both_ITS_LSU))

NC_func_F_rDNA_spc <- NC_func_F %>% left_join(NC_rCNV_tab_spc, by = "spc") %>% filter(!is.na(m_rDNA))
NC_func_F_rDNA_spc_na <- NC_func_F %>% left_join(NC_rCNV_tab_spc, by = "spc") %>% filter(is.na(m_rDNA))

NC_rCNV_tab_gen <- NC_rCNV_tab %>% group_by(gen) %>% summarise(m_rDNA = mean(Both_ITS_LSU))
#NC_rCNV_tab_gen

NC_func_F_rDNA_gen <- NC_func_F_rDNA_spc_na %>% select(-m_rDNA) %>% left_join(NC_rCNV_tab_gen, by = c("Genus" = "gen")) %>% filter(!is.na(m_rDNA))
#NC_func_F_rDNA_gen


NC_func_F1 <- rbind(NC_func_F_rDNA_spc, NC_func_F_rDNA_gen)


NC_func_F_p <- NC_func_F1 %>% gather(key = "Fungal_trait", value = "Fungal_value", extension:CUE)
NC_func_F_p <- NC_func_F_p %>% filter(Fungal_value != "NA")

NC_func_F_p$Fungal_value <- as.numeric(NC_func_F_p$Fungal_value)

NC_func_F_p1 <- NC_func_F_p %>% select(-strainID,-newID)
#NC_func_F_p1
#write_xlsx(NC_func_F_p1, "NC_func_rCNV.xlsx")

#NC_func_F_p1 <- read_excel("./output/NC_func_rCNV.xlsx", sheet = 1)
#NC_func_F_p1$Species %>% unique()

NC_func_F_p1_wide <-
  NC_func_F_p1 %>% spread(key = Fungal_trait, value = Fungal_value) %>% drop_na()

view(NC_func_F_p1_wide)

NC_func_pca <- NC_func_F_p1_wide %>% select(5:45)
NC_func_pca


NC_func_sal <- scale(NC_func_pca, center = T, scale = T)
NC_func_sal_pca <- rda(NC_func_sal)
NC_pca_sum <- summary(NC_func_sal_pca)


NC_func_sal_pca_p_main <- data.frame(
  PC1 = NC_pca_sum$sites[, 1],
  PC2 = NC_pca_sum$sites[, 2],
  Phylum = NC_func_F_p1_wide$Phylum,
  Species = NC_func_F_p1_wide$Species,
  Genus = NC_func_F_p1_wide$Genus
)

NC_func_sal_pca_p_main

NC_func_sal_pca_p_fac <- data.frame(
  PC1 = NC_pca_sum$species[, 1],
  PC2 = NC_pca_sum$species[, 2],
  phylo = rownames(NC_pca_sum$species)
)


NC_func_sal_pca_p_fac <-
  NC_func_sal_pca_p_fac %>%
  mutate(phylo = if_else(phylo == "m_rDNA", "rDNA copy number", phylo)) %>%
  mutate(grp = if_else(phylo == "rDNA copy number", "main", "others"))


rCNV_func_pca_p <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    data = NC_func_sal_pca_p_main,
    aes(x = PC1, y = PC2),
    colour = "black",
    size = 2
  ) +
  geom_segment(
    data = NC_func_sal_pca_p_fac,
    aes(
      x = 0,
      xend = PC1 * 3,
      y = 0,
      yend = PC2 * 3,
      colour = grp
    ),
    arrow =
      arrow(
        angle = 30,
        length = unit(0.2, units = "cm"),
        type = "closed"
      )
  ) +
  geom_text_repel(
    data = NC_func_sal_pca_p_fac,
    aes(
      x = PC1 * 3.3,
      y = PC2 * 3.3,
      label = phylo,
      colour = grp
    ),
    size = 3,
    direction = "y"
  ) +
  scale_colour_manual(values = c("red", "black")) +
  labs(x = "PC1 (18.9%)", y = "PC2 (13.5%)") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(
      colour = "black",
      face = "bold",
      size = 15
    ),
    legend.position = "none"
  )
rCNV_func_pca_p


# lm
# lm model
lm_rlt_Spcs <- function(df, yval, xval, subgrp = c("all")) {
  
  if(subgrp != "all") {
    
    df_tmp <- df %>% filter(Fungal_trait == subgrp)
    
  } else {
    
    df_tmp <- df
    
  }
  
  
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
  
  df_rlt <- data.frame(
    Fungal_trait = subgrp,
    anno_x1 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.5 + range(df_tmp[[xval]])[1],
    anno_y1 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.8 +
      range(df_tmp[[yval]], na.rm = T)[1],
    rsig = str_c("italic(r) == ", lmR_r),
    psig = Rp_sig,
    anno_x2 = (range(df_tmp[[xval]])[2] - range(df_tmp[[xval]])[1]) * 0.6 + range(df_tmp[[xval]])[1],
    anno_y2 = (range(df_tmp[[yval]], na.rm = T)[2] - range(df_tmp[[yval]], na.rm = T)[1]) * 0.1 +
      range(df_tmp[[yval]], na.rm = T)[1],
    pval = lmR_p
  )
  
  pval_sig <- str_c("italic(p) == ", lmR_p)
  df_rlt_F <- df_rlt %>% mutate(rlt_sig = str_c(rsig, " ~~ ", pval_sig))
  
  return(df_rlt_F)
  
}


Fungal_trait_list <- NC_func_F_p1$Fungal_trait %>% unique()

NC_func_lm_subdata <- 
  map_dfr(Fungal_trait_list, ~ lm_rlt_Spcs(NC_func_F_p1,
                                           yval = "Fungal_value",
                                           xval = "m_rDNA",
                                           subgrp = .))


NC_func_lm_subdata1 <- NC_func_lm_subdata %>% filter(psig != "NS")

NC_func_F_p1_2 <- NC_func_F_p1 %>% filter(Fungal_trait %in% NC_func_lm_subdata1$Fungal_trait)

NC_func_F_p1_3 <- 
  NC_func_F_p1_2 %>% mutate(
    Fungal_trait1 = case_when(
      str_detect(Fungal_trait, "hyphal_diam") ~ "hyphal diameter",
      str_detect(Fungal_trait, "enz_pho") ~ "acid phosphatase",
      str_detect(Fungal_trait, "density") ~ "mycelial density",
      str_detect(Fungal_trait, "melanin") ~ "melanin content",
    )
  )

NC_func_lm_subdata1 <- 
  NC_func_lm_subdata1 %>% mutate(
    Fungal_trait1 = case_when(
      str_detect(Fungal_trait, "hyphal_diam") ~ "hyphal diameter",
      str_detect(Fungal_trait, "enz_pho") ~ "acid phosphatase",
      str_detect(Fungal_trait, "density") ~ "mycelial density",
      str_detect(Fungal_trait, "melanin") ~ "melanin content",
    )
  )


NC_func_F_p1_3$Fungal_trait1 <- factor(NC_func_F_p1_3$Fungal_trait1,
                                      levels = c("hyphal diameter", "acid phosphatase",
                                                 "mycelial density", "melanin content"
                                                 ))


NC_p1a <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "hyphal diameter"),
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "hyphal diameter"),
            aes(x = anno_x1, y = anno_y1, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_colour_d3(palette = "category20") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
  labs(x = "rDNA copy number", y = "Hyphal diameter") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        aspect.ratio = 1)
NC_p1a

NC_p1b <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "acid phosphatase"),
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "acid phosphatase"),
            aes(x = anno_x1, y = anno_y1, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_colour_d3(palette = "category20") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
  labs(x = "rDNA copy number", y = "Acid phosphatase") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        aspect.ratio = 1)
NC_p1b


NC_p1c <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "mycelial density"),
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus), size = 3.8) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "mycelial density"),
            aes(x = anno_x1, y = anno_y1, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_colour_d3(palette = "category20") +
  labs(x = "rDNA copy number", y = "Mycelial density") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
NC_p1c


NC_p1d <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"),
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus), size = 3.8) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "melanin content"),
            aes(x = anno_x1, y = anno_y1, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_colour_d3(palette = "category20") +
  labs(x = "rDNA copy number", y = "Melanin content") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
NC_p1d

#fit_lm <- lm(m_rDNA ~ Fungal_value, data = NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"))
#fit_lm$residuals
#fitted.values(fit_lm)
#plot(fitted.values(fit_lm), fit_lm$residuals)

# fig2d-e -----------------
fig2de_trait <- NC_p1c + NC_p1d + plot_layout(ncol = 1)
fig2de_trait

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2de_pdf <- str_c("fig2de_", "trait_", tm, ".pdf", sep = "")
fig2de_jpg <- str_c("fig2de_", "trait_", tm, ".jpg", sep = "")


ggsave(fig2de_pdf, fig2de_trait, width = 6.08, height = 8.91)
ggsave(fig2de_jpg, fig2de_trait, width = 6.08, height = 8.91)

# figS5 -----------------
NC_func_lm_subdata1e <- lm_rlt_Spcs(NC_func_F_p1, yval = "Fungal_value", xval = "m_rDNA", subgrp = "extension")
NC_func_F_pe <- NC_func_F_p1 %>% filter(Fungal_trait == "extension")


NC_p1e <- 
  ggplot(NC_func_F_pe,
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus), size = 3.8) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = NC_func_lm_subdata1e,
            aes(x = anno_x1, y = anno_y1 + 0.2, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_colour_d3(palette = "category20") +
  labs(x = "rDNA copy number", y = "Extension") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1)
NC_p1e

figS5_extension <- NC_p1e

tm <- now() %>% str_split_i(pattern = " ", 1)
figS5_pdf <- str_c("figS5_", "extension_", tm, ".pdf", sep = "")
figS5_jpg <- str_c("figS5_", "extension_", tm, ".jpg", sep = "")

ggsave(figS5_pdf, figS5_extension, width = 6.48, height = 4.9)
ggsave(figS5_jpg, figS5_extension, width = 6.48, height = 4.9)

# done...

