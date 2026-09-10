


# corr test among bacteria genome size and rrn
# by Qiushi-Li, IM-CAS 20260401



# packages
library(tidyverse)

# RRN data
bac_rrnDB <- read_tsv("./2.database/rrnDB-5.9.tsv")
bac_rrnDB

colnames(bac_rrnDB) <- str_replace_all(colnames(bac_rrnDB), pattern = " ", replacement = "_")

bac_rrnDB1 <- bac_rrnDB %>% select(Data_source_record_id, 
                                   Data_source_organism_name, NCBI_scientific_name, NCBI_tax_id,
                                   RDP_taxa, 
                                   RDP_taxonomic_lineage, NCBI_taxonomic_lineage,
                                   basecount, `16S_gene_count`, `23S_gene_count`, tRNA_gene_count)
bac_rrnDB1


# gtdb data
gtdb_bac120 <- read_tsv("./2.database/bac120_metadata_r226.tsv")
gtdb_ar53 <- read_tsv("./2.database/ar53_metadata_r226.tsv")

# view(gtdb_bac120)
# head(gtdb_bac120) %>% view()
# colnames(gtdb_bac120)

# head(gtdb_ar53) %>% view()
gtdb_bac120a <- gtdb_bac120 %>% select(ncbi_taxid, genome_size) %>% group_by(ncbi_taxid) %>% summarise(m_GS = mean(genome_size))

bac_rrnDB2 <- bac_rrnDB1 %>% left_join(gtdb_bac120a, by = c("NCBI_tax_id" = "ncbi_taxid"))
bac_rrnDB2

# dim(bac_rrnDB1)
# dim(bac_rrnDB2)

bac_rrnDB3 <- bac_rrnDB2 %>% select(NCBI_tax_id, NCBI_scientific_name, RDP_taxa, RDP_taxonomic_lineage, `16S_gene_count`, m_GS) %>% drop_na()
bac_rrnDB3

colnames(bac_rrnDB3)[5] <- "gene_count_16S"


bac_rrnDB4 <- bac_rrnDB3 %>% 
  separate_wider_delim(RDP_taxonomic_lineage, names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                       delim = "; ", too_few = "debug", too_many = "debug") %>%  
  mutate(
    Phylum = replace_na(Phylum, replace = "p__unknown"),
    Class = replace_na(Class, replace = "c__unknown"),
    Order = replace_na(Order, replace = "o__unknown"),
    Family = replace_na(Family, replace = "f__unknown"),
    Genus = replace_na(Genus, replace = "g__unknown")
  ) %>% filter(
    !str_detect(Phylum, pattern = " ")
  ) %>% filter(
    str_count(Phylum, pattern = "\\|") == 1
  ) %>% filter(Domain != "Archaea|domain") %>% 
  select(NCBI_tax_id, NCBI_scientific_name, RDP_taxa, Phylum, Class, Order, Family, Genus, gene_count_16S, m_GS) %>% 
  add_count(Phylum) %>%
  filter(n > 10) %>% 
  mutate(Phylum = str_split_i(Phylum, pattern = "\\|", 1)) %>% 
  filter(Phylum != "Rhodothermota")

bac_rrnDB_Phy_cont <- bac_rrnDB4 %>% select(Phylum, n) %>% distinct_all()


# help(add_count)
# bac_rrnDB3$RDP_taxonomic_lineage

# lm test



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


Bac_Phy_list <- bac_rrnDB4$Phylum %>% unique()
colnames(bac_rrnDB4)[4] <- "group"

Bac_Phy_lm_subdata <- 
  subplot_data_corr(DF = bac_rrnDB4, Y_val = "m_GS", X_val = "gene_count_16S",
                    GROUP_list = Bac_Phy_list, METHOD = "spearman", adj_METHOD = "fdr")

subplot_data_corr(DF = bac_rrnDB4, Y_val = "m_GS", X_val = "gene_count_16S",
                  METHOD = "spearman", adj_METHOD = "fdr")


Bac_Phy_lm_subdata1 <- Bac_Phy_lm_subdata %>% left_join(bac_rrnDB_Phy_cont, by = c("group" = "Phylum")) %>% 
  mutate(lab_n = str_c("n = ", n, sep = ""))


bac_rrnDB4_order <- bac_rrnDB4 %>% select(group, n) %>% 
  distinct_all() %>% arrange(desc(n)) %>% select(group) %>% pull()

bac_rrnDB4_order

bac_rrnDB4$group <- factor(bac_rrnDB4$group, levels = bac_rrnDB4_order)
Bac_Phy_lm_subdata1$group <- factor(Bac_Phy_lm_subdata1$group, levels = bac_rrnDB4_order)


Bac_Phy_lm_subdata1


figS11_bac_rrn_GS <- 
  ggplot(bac_rrnDB4, aes(gene_count_16S, m_GS/1000000)) +
  geom_jitter(size = 0.4, alpha = 0.3, width = 0.1) +
  # geom_smooth(colour = "red", method = "lm") +
  # geom_hex() +
  geom_smooth(method = "loess", span = 1, colour = "red") +
  geom_text(data = Bac_Phy_lm_subdata1,
            aes(x = 12, y = Inf, label = rlt_adj_sig),
            parse = T,
            colour = "blue",
            size = 4.5,
            vjust = 2) +
  geom_text(data = Bac_Phy_lm_subdata1,
            aes(x = 15, y = -Inf, label = lab_n),
            colour = "blue",
            size = 4.5,
            vjust = -2) +
  #scale_fill_gradient2(low = "#3B4992FF", high = "#008B45FF", mid = "white", midpoint = 2500) +
  #guides(fill = guide_colorbar(override.aes = list(size = 5, alpha = 1),
  #                             nrow = 1,
  #                             direction = "horizontal")) +
  facet_wrap(~ group, scales = "free") +
  scale_x_continuous(limits = c(0, 24), labels = seq(0, 24, 6),
                     breaks = seq(0, 24, 6)) +
  labs(x = "rDNA copy number", y = "Assembly genome size (Mbp)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        #legend.position = "inside",
        #legend.position.inside = c(0.75, 0.1),
        aspect.ratio = 1)
  # facetted_pos_scales(
  #   x = list(
  #     Phylum == "Actinomycetota" ~ scale_x_continuous(labels = seq(0, 12, 2),
  #                                                 breaks = seq(0, 12, 2)),
  #     Phylum == "Fusobacteriota" ~ scale_x_continuous(labels = seq(0, 12, 2),
  #                                                     breaks = seq(0, 12, 2)),
  #     Phylum == "Acidobacteriota" ~ scale_x_continuous(labels = seq(0, 3, 1),
  #                                                     breaks = seq(0, 3, 1)),
  #     Phylum == "Chlorobiota" ~ scale_x_continuous(labels = seq(0, 3, 1),
  #                                                     breaks = seq(0, 3, 1)),
  #     Phylum == "Aquificota" ~ scale_x_continuous(labels = seq(0, 3, 1),
  #                                                     breaks = seq(0, 3, 1)),
  #     Phylum == "Nitrospirota" ~ scale_x_continuous(labels = seq(0, 3, 1),
  #                                                     breaks = seq(0, 3, 1))
  #   )
  # )
figS11_bac_rrn_GS


tm <- now() %>% str_split_i(pattern = " ", 1)
figS11_pdf <- str_c("figS11_", "Bac_rrn_GS_", tm, ".pdf", sep = "")
figS11_jpg <- str_c("figS11_", "Bac_rrn_GS_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S11_pdf <- str_c(fig_path, figS11_pdf)
fig_fullpath_S11_jpg <- str_c(fig_path, figS11_jpg)

ggsave(fig_fullpath_S11_pdf, figS11_bac_rrn_GS, width = 15, height = 12)
ggsave(fig_fullpath_S11_jpg, figS11_bac_rrn_GS, width = 15, height = 12)


# ggsave("bac_Phy_GS_10.pdf", p_bac_rrn_GS, width = 15, height = 12)

# bac_rrnDB2$Data_source_organism_name

p_bac_rrn_GS1 <- p_bac_rrn_GS + facet_wrap(~ group, nrow = 2, scales = "free")
p_bac_rrn_GS1


ggsave("p_bac_rrn_GS1.jpg")

# save.image("FRRN_v0.96.RData")




