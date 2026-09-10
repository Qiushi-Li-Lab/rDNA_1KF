

# FRRN ~ Fungal traits ----------
# fig2d-e #
# Fungal phenotype database were obtained from Camenzind et al (2024) # 
# by Qiushi-Li, IM-CAS, 2025.02.20

##### packages --------------
library(tidyverse)
library(readxl)
library(ggsci)
library(vegan)
library(ggrepel)
library(ggh4x)
library(patchwork)


# Camenzind data
Fungal_phenotype_data <- read_excel("./2.database/Camenzind_trait data_NatComm.xlsx", sheet = 2)

Fungal_phenotype_data1 <- 
  Fungal_phenotype_data %>% mutate(Genus = str_split_i(Species, pattern = "_", 1)) %>%
  mutate(spc = str_c(
    str_split_i(Species, pattern = "_", 1),
    str_split_i(Species, pattern = "_", 2),
    sep = " ")) %>% 
  drop_na()

# Camenzind proj
Fungal_phenotype_proj <- read_excel("./2.database/NC_project.xlsx", sheet = 1)

proj_list <- Fungal_phenotype_proj %>% select(Name)
proj_list

# FRRN data
FRRN_rlt_taxa_spl_GSGN %>% colnames()

NC_FRRN_tab <- 
  FRRN_rlt_taxa_spl_GSGN %>% 
  filter(Name.x %in% proj_list$Name) %>% 
  mutate(
    cla = str_sub(cla, start = 3),
    ord = str_sub(ord, start = 3),
    fam = str_sub(fam, start = 3),
    gen = str_sub(Gen, start = 3)
  ) %>% 
  mutate(spc = str_c(
    str_split_i(Name.x, pattern = " ", 1),
    str_split_i(Name.x, pattern = " ", 2),
    sep = " ")) %>% select(project, Both_ITS_LSU, Assembly_Length, gen, spc)


NC_FRRN_tab_spc <- NC_FRRN_tab %>% group_by(spc) %>% summarise(m_rDNA = mean(Both_ITS_LSU),
                                                               m_GS = mean(Assembly_Length))

NC_func_F_rDNA_spc <- Fungal_phenotype_data1 %>% 
  left_join(NC_FRRN_tab_spc, by = "spc") %>% filter(!is.na(m_rDNA))

NC_func_F_rDNA_spc_na <- Fungal_phenotype_data1 %>% 
  left_join(NC_FRRN_tab_spc, by = "spc") %>% filter(is.na(m_rDNA))

NC_FRRN_tab_gen <- NC_FRRN_tab %>% group_by(gen) %>% 
  summarise(m_rDNA = mean(Both_ITS_LSU),
            m_GS = mean(Assembly_Length))
#NC_rCNV_tab_gen

NC_func_F_rDNA_gen <- NC_func_F_rDNA_spc_na %>% 
  select(-m_rDNA, -m_GS) %>% 
  left_join(NC_FRRN_tab_gen, by = c("Genus" = "gen")) %>% filter(!is.na(m_rDNA))
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

# view(NC_func_F_p1_wide)

# NC_func_pca <- NC_func_F_p1_wide %>% select(5:45)
# NC_func_pca
# 
# 
# NC_func_sal <- scale(NC_func_pca, center = T, scale = T)
# # NC_func_sal_pca <- rda(NC_func_sal)
# # NC_pca_sum <- summary(NC_func_sal_pca)
# 
# 
# NC_func_sal_pca_p_main <- data.frame(
#   PC1 = NC_pca_sum$sites[, 1],
#   PC2 = NC_pca_sum$sites[, 2],
#   Phylum = NC_func_F_p1_wide$Phylum,
#   Species = NC_func_F_p1_wide$Species,
#   Genus = NC_func_F_p1_wide$Genus
# )
# 
# NC_func_sal_pca_p_main
# 
# NC_func_sal_pca_p_fac <- data.frame(
#   PC1 = NC_pca_sum$species[, 1],
#   PC2 = NC_pca_sum$species[, 2],
#   phylo = rownames(NC_pca_sum$species)
# )
# 
# 
# NC_func_sal_pca_p_fac <-
#   NC_func_sal_pca_p_fac %>%
#   mutate(phylo = if_else(phylo == "m_rDNA", "rDNA copy number", phylo)) %>%
#   mutate(grp = if_else(phylo == "rDNA copy number", "main", "others"))
# 
# 
# rCNV_func_pca_p <-
#   ggplot() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_point(
#     data = NC_func_sal_pca_p_main,
#     aes(x = PC1, y = PC2),
#     colour = "black",
#     size = 2
#   ) +
#   geom_segment(
#     data = NC_func_sal_pca_p_fac,
#     aes(
#       x = 0,
#       xend = PC1 * 3,
#       y = 0,
#       yend = PC2 * 3,
#       colour = grp
#     ),
#     arrow =
#       arrow(
#         angle = 30,
#         length = unit(0.2, units = "cm"),
#         type = "closed"
#       )
#   ) +
#   geom_text_repel(
#     data = NC_func_sal_pca_p_fac,
#     aes(
#       x = PC1 * 3.3,
#       y = PC2 * 3.3,
#       label = phylo,
#       colour = grp
#     ),
#     size = 3,
#     direction = "y"
#   ) +
#   scale_colour_manual(values = c("red", "black")) +
#   labs(x = "PC1 (18.9%)", y = "PC2 (13.5%)") +
#   theme_bw() +
#   theme(
#     axis.text = element_text(colour = "black", size = 12),
#     axis.title = element_text(
#       colour = "black",
#       face = "bold",
#       size = 15
#     ),
#     legend.position = "none"
#   )
# rCNV_func_pca_p



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


colnames(NC_func_F_p1)[7] <- "group"
Fungal_trait_list <- NC_func_F_p1$group %>% unique()


FRRN_trait_cor_subdata <- 
  subplot_data_corr(DF = NC_func_F_p1, Y_val = "Fungal_value", X_val = "m_rDNA",
                    GROUP_list = Fungal_trait_list, METHOD = "spearman", adj_METHOD = "fdr")
FRRN_trait_cor_subdata



# NC_func_lm_subdata <- 
#   map_dfr(Fungal_trait_list, ~ lm_rlt_Spcs(NC_func_F_p1,
#                                            yval = "Fungal_value",
#                                            xval = "m_rDNA",
#                                            subgrp = .))



FRRN_trait_cor_subdata1 <- FRRN_trait_cor_subdata %>% filter(pval_sig != "NS")
# NC_func_lm_subdata1

NC_func_F_p1_2 <- NC_func_F_p1 %>% filter(group %in% FRRN_trait_cor_subdata1$group)

NC_func_F_p1_3 <- 
  NC_func_F_p1_2 %>% mutate(
    Fungal_trait1 = case_when(
      str_detect(group, "extension") ~ "extension",
      str_detect(group, "enz_leu") ~ "enz_leu",
      str_detect(group, "density") ~ "mycelial density",
      str_detect(group, "melanin") ~ "melanin content",
    )
  )

FRRN_trait_cor_subdata1 <- 
  FRRN_trait_cor_subdata1 %>% mutate(
    Fungal_trait1 = case_when(
      str_detect(group, "extension") ~ "extension",
      str_detect(group, "enz_leu") ~ "enz_leu",
      str_detect(group, "density") ~ "mycelial density",
      str_detect(group, "melanin") ~ "melanin content",
    )
  )


NC_func_F_p1_3$Fungal_trait1 <- factor(NC_func_F_p1_3$Fungal_trait1,
                                      levels = c("extension", "enz_leu",
                                                 "mycelial density", "melanin content"
                                                 ))

# view(NC_func_F_p1_3)

# NC_p1a <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "extension"),
#          aes(m_rDNA, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "extension"),
#             aes(x = anno_x1, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
#   labs(x = "rDNA copy number", y = "Hyphal diameter") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "none",
#         aspect.ratio = 1)
# NC_p1a

# NC_p1b <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "enz_leu"),
#          aes(m_rDNA, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "enz_leu"),
#             aes(x = anno_x1, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
#   labs(x = "rDNA copy number", y = "Acid phosphatase") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "none",
#         aspect.ratio = 1)
# NC_p1b


FRRN_phenotype_GS_fig2d <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "mycelial density"),
         aes(m_rDNA, Fungal_value)) +
 # geom_point(aes(colour = Genus), size = 3.8) +
  geom_point(aes(colour = Genus, size = m_GS/1000000)) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = FRRN_trait_cor_subdata1 %>% filter(Fungal_trait1 == "mycelial density"),
            aes(x = anno_x1, y = Inf, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5,
            vjust = 3.5) +
  geom_text_repel(aes(label = round(m_GS/1000000, 0)), size = 3) +
  scale_size(range = c(0.5, 5)) +
  # guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3.8),
  #                              label.theme = element_text(face = "italic")),
  #        size = guide_legend(title = "Assembly\ngenome size (Mbp)", override.aes = list(colour = "gray"),
  #                            label.theme = element_text(face = "plain"))) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 8),
                               label.theme = element_text(face = "italic")),
         size = guide_legend(title = "Assembly\ngenome size (Mbp)", override.aes = list(colour = "gray"),
                             label.theme = element_text(face = "plain"))) +
  scale_colour_d3(palette = "category20") +
  labs(x = "rDNA copy number", y = "Mycelial density") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(face = "italic", size = 15),
        legend.key.size = unit(0.9, "cm"),
        legend.key.spacing.y = unit(0.3, units = "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)

FRRN_phenotype_GS_fig2d


FRRN_phenotype_GS_fig2e <- 
  ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"),
         aes(m_rDNA, Fungal_value)) +
  # geom_point(aes(colour = Genus), size = 3.8) +
  geom_point(aes(colour = Genus, size = m_GS/1000000)) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = FRRN_trait_cor_subdata1 %>% filter(Fungal_trait1 == "melanin content"),
            aes(x = anno_x1, y = Inf, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5,
            vjust = 3.5) +
  geom_text_repel(aes(label = round(m_GS/1000000, 0)), size = 3) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3.8),
                               label.theme = element_text(face = "italic")),
         size = guide_legend(title = "Assembly\ngenome size (Mbp)", override.aes = list(colour = "gray"),
                             label.theme = element_text(face = "plain"))) +
  scale_colour_d3(palette = "category20") +
  scale_size(range = c(0.5, 5)) +
  labs(x = "rDNA copy number", y = "Melanin content") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none",
        aspect.ratio = 1)
FRRN_phenotype_GS_fig2e

#fit_lm <- lm(m_rDNA ~ Fungal_value, data = NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"))
#fit_lm$residuals
#fitted.values(fit_lm)
#plot(fitted.values(fit_lm), fit_lm$residuals)

# fig2d-e -----------------
fig2de_trait <- FRRN_phenotype_GS_fig2d + FRRN_phenotype_GS_fig2e + plot_layout(ncol = 1, guides = "collect")
fig2de_trait

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2de_pdf <- str_c("fig2de_", "trait_", tm, ".pdf", sep = "")
fig2de_jpg <- str_c("fig2de_", "trait_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_2de_pdf <- str_c(fig_path, fig2de_pdf)
fig_fullpath_2de_jpg <- str_c(fig_path, fig2de_jpg)

ggsave(fig_fullpath_2de_pdf, fig2de_trait, width = 6.08, height = 8.91)
ggsave(fig_fullpath_2de_jpg, fig2de_trait, width = 6.08, height = 8.91)



# figS5 -----------------
FRRN_trait_cor_subdata_ext <- 
  subplot_data_corr(DF = NC_func_F_p1, Y_val = "Fungal_value", X_val = "m_rDNA",
                    GROUP_list = "extension", METHOD = "spearman", adj_METHOD = "fdr")
FRRN_trait_cor_subdata_ext


NC_func_F_pe <- NC_func_F_p1 %>% filter(group == "extension")


FRRN_phenotype_figS5 <- 
  ggplot(NC_func_F_pe,
         aes(m_rDNA, Fungal_value)) +
  geom_point(aes(colour = Genus, size = m_GS/1000000)) +
  geom_smooth(method = "lm", colour = "red") +
  geom_text(data = FRRN_trait_cor_subdata_ext,
            aes(x = anno_x1, y = Inf, label = rlt_sig),
            parse = T,
            colour = "blue",
            size = 5.5,
            vjust = 3.5) +
  geom_text_repel(aes(label = round(m_GS/1000000, 0)), size = 3) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3.8),
                               label.theme = element_text(face = "italic"),
                               order = 1),
         size = guide_legend(title = "Assembly\ngenome size (Mbp)", override.aes = list(colour = "gray"),
                             label.theme = element_text(face = "plain"),
                             order = 2),
         ) +
  scale_colour_d3(palette = "category20") +
  scale_size(range = c(0.5, 5)) +
  labs(x = "rDNA copy number", y = "Mycelial extension") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        legend.box = "horizontal",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1)
FRRN_phenotype_figS5

tm <- now() %>% str_split_i(pattern = " ", 1)
figS5_pdf <- str_c("figS5_", "extension_", tm, ".pdf", sep = "")
figS5_jpg <- str_c("figS5_", "extension_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S5_pdf <- str_c(fig_path, figS5_pdf)
fig_fullpath_S5_jpg <- str_c(fig_path, figS5_jpg)


ggsave(fig_fullpath_S5_pdf, FRRN_phenotype_figS5, width = 8, height = 4.9)
ggsave(fig_fullpath_S5_jpg, FRRN_phenotype_figS5, width = 8, height = 4.9)

# done...



#### continuous

# Fungal_trait_list <- NC_func_F_p1$Fungal_trait %>% unique()
# 
# NC_func_lm_GS_subdata <- 
#   map_dfr(Fungal_trait_list, ~ lm_rlt_Spcs(NC_func_F_p1,
#                                            yval = "Fungal_value",
#                                            xval = "m_GS",
#                                            subgrp = .))
# 
# NC_func_lm_GS_subdata
# NC_func_lm_GS_subdata1 <- NC_func_lm_GS_subdata %>% filter(psig != "NS")
# 
# # NC_func_F_p1_2 <- NC_func_F_p1 %>% filter(Fungal_trait %in% NC_func_lm_subdata1$Fungal_trait)
# 
# NC_func_F_p1_3 <- 
#   NC_func_F_p1 %>% mutate(
#     Fungal_trait1 = case_when(
#       str_detect(Fungal_trait, "extension") ~ "extension",
#       str_detect(Fungal_trait, "enz_leu") ~ "enz_leu",
#       str_detect(Fungal_trait, "density") ~ "mycelial density",
#       str_detect(Fungal_trait, "melanin") ~ "melanin content",
#       .default = Fungal_trait
#     )
#   )
# 
# NC_func_lm_GS_subdata2 <- 
#   NC_func_lm_GS_subdata %>% mutate(
#     Fungal_trait = case_when(
#       str_detect(Fungal_trait, "extension") ~ "extension",
#       str_detect(Fungal_trait, "enz_leu") ~ "enz_leu",
#       str_detect(Fungal_trait, "density") ~ "mycelial density",
#       str_detect(Fungal_trait, "melanin") ~ "melanin content",
#       .default = Fungal_trait
#     )
#   )


#NC_func_F_p1_3$Fungal_trait1 <- factor(NC_func_F_p1_3$Fungal_trait1,
#                                       levels = c("extension", "enz_leu",
#                                                  "mycelial density", "melanin content"
#                                       ))

# view(NC_func_F_p1_3)

# NC_p1a <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "extension"),
#          aes(m_rDNA, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "extension"),
#             aes(x = anno_x1, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
#   labs(x = "rDNA copy number", y = "Hyphal diameter") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "none",
#         aspect.ratio = 1)
# NC_p1a

# NC_p1b <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "enz_leu"),
#          aes(m_rDNA, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8, alpha = 0.7) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_subdata1 %>% filter(Fungal_trait1 == "enz_leu"),
#             aes(x = anno_x1, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
#   scale_y_continuous(expand = expansion(mult = c(0.1, 0.28))) +
#   labs(x = "rDNA copy number", y = "Acid phosphatase") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "none",
#         aspect.ratio = 1)
# NC_p1b


# NC_p1c_GS <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "mycelial density"),
#          aes(m_GS/1000000, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_GS_subdata2 %>% filter(Fungal_trait == "mycelial density"),
#             aes(x = anno_x1/1000000, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   labs(x = "Assembly genome size (Mbp)", y = "Mycelial density") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "right",
#         aspect.ratio = 1)
# NC_p1c_GS
# 
# 
# NC_p1d_GS <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"),
#          aes(m_GS/1000000, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_GS_subdata2 %>% filter(Fungal_trait == "melanin content"),
#             aes(x = anno_x1/1000000, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   labs(x = "Assembly genome size (Mbp)", y = "Melanin content") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "right",
#         aspect.ratio = 1)
# NC_p1d_GS

#fit_lm <- lm(m_rDNA ~ Fungal_value, data = NC_func_F_p1_3 %>% filter(Fungal_trait1 == "melanin content"))
#fit_lm$residuals
#fitted.values(fit_lm)
#plot(fitted.values(fit_lm), fit_lm$residuals)

# fig2d-e -----------------
# fig2de_trait <- NC_p1c + NC_p1d + plot_layout(ncol = 1)
# fig2de_trait
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# fig2de_pdf <- str_c("fig2de_", "trait_", tm, ".pdf", sep = "")
# fig2de_jpg <- str_c("fig2de_", "trait_", tm, ".jpg", sep = "")
# 
# 
# ggsave(fig2de_pdf, fig2de_trait, width = 6.08, height = 8.91)
# ggsave(fig2de_jpg, fig2de_trait, width = 6.08, height = 8.91)

# figS5 -----------------
# NC_func_lm_subdata1e <- lm_rlt_Spcs(NC_func_F_p1, yval = "Fungal_value", xval = "m_rDNA", subgrp = "extension")
# NC_func_F_pe <- NC_func_F_p1 %>% filter(Fungal_trait == "extension")
# 
# 
# NC_p1e <- 
#   ggplot(NC_func_F_pe,
#          aes(m_rDNA, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_subdata1e,
#             aes(x = anno_x1, y = anno_y1 + 0.2, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   labs(x = "rDNA copy number", y = "Extension") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         aspect.ratio = 1)
# NC_p1e
# 
# figS5_extension <- NC_p1e
# 
# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS5_pdf <- str_c("figS5_", "extension_", tm, ".pdf", sep = "")
# figS5_jpg <- str_c("figS5_", "extension_", tm, ".jpg", sep = "")
# 
# ggsave(figS5_pdf, figS5_extension, width = 6.48, height = 4.9)
# ggsave(figS5_jpg, figS5_extension, width = 6.48, height = 4.9)

# done...


# NC_p1e_GS <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "hyphal_diam"),
#          aes(m_GS/1000000, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_GS_subdata2 %>% filter(Fungal_trait == "hyphal_diam"),
#             aes(x = anno_x1/1000000, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   labs(x = "Assembly genome size (Mbp)", y = "hyphal_diam") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "right",
#         aspect.ratio = 1)
# NC_p1e_GS
# 
# 
# NC_p1f_GS <- 
#   ggplot(NC_func_F_p1_3 %>% filter(Fungal_trait1 == "spore_shape"),
#          aes(m_GS/1000000, Fungal_value)) +
#   geom_point(aes(colour = Genus), size = 3.8) +
#   geom_smooth(method = "lm", colour = "red") +
#   geom_text(data = NC_func_lm_GS_subdata2 %>% filter(Fungal_trait == "spore_shape"),
#             aes(x = anno_x1/1000000, y = anno_y1, label = rlt_sig),
#             parse = T,
#             colour = "blue",
#             size = 5.5) +
#   guides(colour = guide_legend(ncol = 1)) +
#   scale_colour_d3(palette = "category20") +
#   labs(x = "Assembly genome size (Mbp)", y = "spore_shape") +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
#         axis.text = element_text(colour = "black", size = 12),
#         legend.title = element_text(face = "bold"),
#         legend.text = element_text(face = "italic"),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         legend.position = "right",
#         aspect.ratio = 1)
# NC_p1f_GS


# fig2d-e -----------------
FRRN_phenotype_figS5_2f <- FRRN_phenotype_figS5 + theme(legend.position = "none")
fig2def_trait <- FRRN_phenotype_figS5_2f + FRRN_phenotype_GS_fig2d + FRRN_phenotype_GS_fig2e + plot_layout(ncol = 1, guides = "collect")
fig2def_trait

tm <- now() %>% str_split_i(pattern = " ", 1)
fig2def_pdf <- str_c("fig2def_", "trait_", tm, ".pdf", sep = "")
fig2def_jpg <- str_c("fig2def_", "trait_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_2def_pdf <- str_c(fig_path, fig2def_pdf)
fig_fullpath_2def_jpg <- str_c(fig_path, fig2def_jpg)

ggsave(fig_fullpath_2def_pdf, fig2def_trait, width = 7, height = 13.5)
ggsave(fig_fullpath_2def_jpg, fig2def_trait, width = 7, height = 13.5)

# 9+4.5

# new layout style


FRRN_phenotype_GS_fig2da <- FRRN_phenotype_GS_fig2d + 
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5),
                               label.theme = element_text(face = "italic"), order = 1),
         size = guide_legend(title = "Assembly\ngenome size (Mbp)", override.aes = list(colour = "gray"),
                             label.theme = element_text(face = "plain"), order = 2)) +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 15),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "italic", size = 10),
        legend.key.size = unit(0.05, "cm"),
        legend.key.spacing.y = unit(0.05, units = "cm"),
        legend.box = "horizontal",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
FRRN_phenotype_GS_fig2da


FRRN_phenotype_GS_fig2ea <- FRRN_phenotype_GS_fig2e + theme(legend.box = "horizontal")
FRRN_phenotype_figS5_2fa <- FRRN_phenotype_figS5_2f + theme(legend.box = "horizontal")

fig2def_trait_new <- FRRN_phenotype_GS_fig2da + FRRN_phenotype_GS_fig2ea + FRRN_phenotype_figS5_2fa + 
  plot_layout(nrow = 1, guides = "collect")
fig2def_trait_new

ggsave("fig2def.pdf")
ggsave("fig2def.jpg")



fig2def_trait_new1 <- FRRN_phenotype_figS5_2fa + FRRN_phenotype_GS_fig2da + FRRN_phenotype_GS_fig2ea + 
  plot_layout(nrow = 1, guides = "collect")
fig2def_trait_new1
ggsave("fig2def_reviews_20260823.pdf")
ggsave("fig2def_reviews_20260823.jpg")
# Saving 16.9 x 5.69 in image






