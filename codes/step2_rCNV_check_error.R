

# rDNA copynumber variation --- check
# figS1
# 2024.11.10

# packages
library(tidyverse)
library(ggrepel)
library(readxl)
# library(writexl)

# load RData
# load("rCNV_version_1.0_20250219.RData")

rDNA_ME <- read_excel("./rDNA_CNV-master/data/Table_S1_tmp.xlsx", sheet = 1)
rDNA_ME_tmp <- rDNA_ME %>% select(Name, project_id)

# rCNV rlt summary
rCNV_rlt_taxa_spl %>% filter(Quality == "Q20")

rCNV_MEcheck <- rDNA_ME_tmp %>% left_join(rCNV_rlt_taxa_spl, by = c("project_id" = "project")) %>%
  select(project_id, ITS_only, Both_ITS_LSU) %>% drop_na()

# ME - rlt

rDNA_ME_rlt <- read_excel("./rDNA_CNV-master/data/Table_S2.xlsx", sheet = 2)
rDNA_ME_rlt %>% view()
rCNV_MEcheck

rDNA_ME_tmp$rDNA_cpnumber <- rDNA_ME_rlt$`rDNA copy number`

rDNA_MEcheck_p <- 
  rDNA_ME_tmp %>% mutate(rDNA_cpnumber = str_split_i(rDNA_cpnumber, pattern = " ", 1)) %>%
  left_join(rCNV_MEcheck, by = "project_id") %>% drop_na()


rDNA_MEcheck_p$rDNA_cpnumber <- as.numeric(rDNA_MEcheck_p$rDNA_cpnumber)
rDNA_MEcheck_lm <- lm(log10(rDNA_cpnumber) ~ log10(Both_ITS_LSU), data = rDNA_MEcheck_p)
rDNA_MEcheck_lm_sum <- summary(rDNA_MEcheck_lm)

rDNA_MEcheck_lm_sum$r.squared %>% sqrt()
rDNA_MEcheck_lm_sum$coefficients
# r = 0.971
# p = 1.306362e-50

rDNA_MEcheck_p <- rDNA_MEcheck_p %>% mutate(spc = str_c(
  str_split_i(Name, i = 1, pattern = " "),
  str_split_i(Name, i = 2, pattern = " "),
  sep = " ")
)



# figS1, check ----------
figS1_check_p1 <- 
  ggplot(rDNA_MEcheck_p, aes(log10(Both_ITS_LSU), log10(rDNA_cpnumber))) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, colour = "grey70") +
  geom_point(colour = "darkgreen", size = 2.5) +
  geom_smooth(method = "lm", colour = "darkorange") +
  geom_text_repel(aes(label = spc), max.overlaps = 25, fontface = "italic") +
  scale_x_continuous(limits = c(1, 3.3),
                     breaks = seq(0, 3.3, 0.5),
                     labels = floor(10 ^ seq(0, 3.3, 0.5))) +
  scale_y_continuous(limits = c(1, 3.3),
                     breaks = seq(0, 3.3, 0.5),
                     labels = floor(10 ^ seq(0, 3.3, 0.5))) +
  labs(x = "Our results (log10-transformed)", y = "Results of Lofgren et al. (log10-transformed)") +
  annotate("text", x = 1.9, y = 3,
           label = expression(italic(r) == 0.971 ~~ italic(p) == 1.306e-50),
           colour = "red",
           size = 12) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 24),
        axis.text = element_text(size = 20, colour = "black"),
        panel.background = element_rect(colour = "black"),
        aspect.ratio = 1)
figS1_check_p1

tm <- now() %>% str_split_i(pattern = " ", 1)
figS1_pdf <- str_c("figS1_", "check_", tm, ".pdf", sep = "")
figS1_jpg <- str_c("figS1_", "check_", tm, ".jpg", sep = "")

ggsave(figS1_pdf, figS1_check_p1, width = 13.2, height = 12.8)
ggsave(figS1_jpg, figS1_check_p1, width = 13.2, height = 12.8)

# done.

