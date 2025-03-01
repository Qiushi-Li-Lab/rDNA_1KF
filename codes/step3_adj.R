

# rCNV fig3 -----------
# by Peilin-Chen and Qiushi-Li, IM-CAS

# packages
library(tidyverse)
library(vegan)
library(readxl)
library(splitstackshape)
library(broom)
library(reshape2)
library(ape)
library(ggsignif)
library(lme4)
library(lmerTest)
library(colorRamps)


# data
fung0 <- read.csv("./rRNA/fung0.csv", header = T, row.names = 1)
#Fung1 <- subset(fung0, fung0$Kingdom=="Fungi")

Fung1 <- fung0 %>% filter(Kingdom == "Fungi")

row.names(Fung1) <- Fung1$OTU

zg <- read.csv("./rRNA/Group_OTU.csv", header = T, row.names = 1)
nrow(zg)

#zg <- subset(zg, ID_OTU != '#N/A')
#zg <- subset(zg, zg$Habitat == 'Root')
zg <- zg %>% filter(ID_OTU != '#N/A') %>% filter(Habitat == "Root")

nrow(zg)

#fungi <- Fung1[,(1:1251)]
#fungi <- data.frame(t(fungi))
fungi <- Fung1 %>% select(1:1251) %>% t() %>% as.data.frame()

row.names(zg) <- zg$ID_OTU
#fungi_g <- as.data.frame(fungi[row.names(fungi)%in%row.names(zg),])
#dim(fungi_g)
#fungi_g <- data.frame(t(fungi_g))
#fungi_g$ID <- row.names(fungi_g)

fungi_g <- fungi %>% filter(row.names(fungi) %in% row.names(zg)) %>% t() %>% as.data.frame() %>% 
  mutate(ID = row.names(.))
dim(fungi_g)

# tax
tax <- read.csv("./rRNA/All_otu_blast1.csv", header = T)

#fungi_g <- fungi_g %>% left_join(tax, by = 'ID')
#fungi_g <- subset(fungi_g, fungi_g$Kingdom=="Fungi")
fungi_g <- fungi_g %>% left_join(tax, by = "ID") %>% filter(Kingdom == "Fungi")

row.names(fungi_g) <- fungi_g$ID

#fungi_g1 <- data.frame(t(fungi_g[,c(1:196)]))
fungi_g1 <- fungi_g %>% select(1:196) %>% t() %>% as.data.frame()

#otu_Flattening = as.data.frame(t(rrarefy(fungi_g1, min(rowSums(fungi_g1)))))
#otu_Flattening$ID <- row.names(otu_Flattening)
#otu_Flattening <- left_join(otu_Flattening,tax,by='ID')

otu_Flattening <- fungi_g1 %>% rrarefy(min(rowSums(fungi))) %>% t() %>% as.data.frame() %>% 
  mutate(ID = row.names(.)) %>% left_join(tax, by = "ID")

# rRNA ----------------
rRNA <- rCNV_taxa_FG

#rRNA1 <-aggregate(rRNA$Both_ITS_LSU,by=list(rRNA$phy), mean)
rRNA1 <- rRNA %>% group_by(phy) %>% summarise(x = mean(Both_ITS_LSU))

colnames(rRNA1)[1] <- 'Subphylum'

#otu_Flattening <- left_join(otu_Flattening,rRNA1,by='Subphylum')
otu_Flattening <- otu_Flattening %>% left_join(rRNA1, by = "Subphylum")
otu_Flattening$x[is.na(otu_Flattening$x)] <-  mean(rRNA$Both_ITS_LSU)

row.names(otu_Flattening) <- otu_Flattening$ID

# adj ---------
unadj_otu <- otu_Flattening[, c(1:196)]
adj_otu <- otu_Flattening[, c(1:196)] / otu_Flattening$x
adj_otu <- data.frame(t(adj_otu))

total <- apply(adj_otu, 1, sum) 
adj_otu.relabu <- data.frame(lapply(adj_otu, function(x) {  x / total  }) )
adj_otu.relabu <- data.frame(t(adj_otu.relabu))

adj_otu.relabu$ID <- row.names(adj_otu.relabu)
adj_otu.relabu <- left_join(adj_otu.relabu, tax, by = 'ID')

adj_otu.relabu_class <- aggregate(adj_otu.relabu[, c(1:196)], by = list(adj_otu.relabu$Class), sum)

row.names(adj_otu.relabu_class) <- adj_otu.relabu_class$Group.1

adj_otu.relabu_class <- adj_otu.relabu_class[, -1]
adj_otu.relabu_class <- data.frame(t(adj_otu.relabu_class))
adj_otu.relabu_class$ID_OTU <- row.names(adj_otu.relabu_class)

adj_glo <- left_join(zg,adj_otu.relabu_class, by = 'ID_OTU')


pp1 <- read.csv("./rRNA/AMS_otu_trans_rel_content.csv", header = T, row.names = 1)
pp2 <- left_join(pp1, adj_glo[,c(1,21)], by='ID1')

# plot ----------------

# fig4e, correct  --------------
fit <- lm(formula = Arbuscular_mycorrhizal_trans_RA ~ 0+Glomeromycetes, data = pp2)
fit_sum <- summary(fit)
fit_sum$r.squared %>% sqrt()

fit %>% tidy()

#b <- as.character(as.expression( substitute(italic(R) ~ "=" ~ 0.761 "," ~~ italic(p) ~ "=" ~ 1.60e-38)))


fig4e_correct_AMF <- 
  ggplot(pp2, aes(Glomeromycetes, Arbuscular_mycorrhizal_trans_RA, colour = as.factor(Week))) + 
  geom_point(size = 3) +
  stat_function(fun = function(x) 1.264 * x, colour = 'red', size = 1, alpha = 0.5) +
  annotate("text", x = 0.25, y = 0.95, parse = TRUE, label = "y==1.264*x", size = 4.5, colour = 'red') + 
  annotate("text", x = 0.35, y = 0.87, parse = TRUE,
           label = expression(italic(r) == 0.761 ~~ italic(p) == 1.60e-38),
           colour = 'red', size = 4.5) +
  guides(colour = guide_legend(title = "Week", ncol = 1)) +
  scale_colour_manual(values = blue2red(18)[3:18]) +
  ylim(0, 1) + xlim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1.5, colour = 'grey80') +
  theme_bw() +
  theme(strip.text = element_text(size = 12,face="bold"),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        aspect.ratio = 1) + 
  labs(y = "AM fungal relative abundance in metatranscriptome",
       x = "AM fungal relative abundance in DNA amplicon (rrn corrected)")
fig4e_correct_AMF

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4e_pdf <- str_c("fig4e_", "correct_liner_Glomeromycetes_", tm, ".pdf", sep = "")
fig4e_jpg <- str_c("fig4e_", "correct_liner_Glomeromycetes_", tm, ".jpg", sep = "")

ggsave(fig4e_pdf, fig4e_correct_AMF, width = 6.3, height = 5.29)
ggsave(fig4e_jpg, fig4e_correct_AMF, width = 6.3, height = 5.29)


# fig4d, uncorrect -------
fit <-lm(formula = Arbuscular_mycorrhizal_trans_RA ~ 0 + Arbuscular_mycorrhizal_OTU_RA, data = pp2)
fit_sum <- summary(fit)
fit_sum$r.squared %>% sqrt()
fit %>% tidy()

#a <- as.character(as.expression( substitute(italic(R)^2 ~ "=" ~ 0.6167*","~~italic(p) ~ "<" ~ 2.2e-16)))


fig4d_uncorrect_AMF <- 
  ggplot(pp2, aes(Arbuscular_mycorrhizal_OTU_RA, Arbuscular_mycorrhizal_trans_RA, colour = as.factor(Week))) +
  geom_point(size = 3) +
  stat_function(fun = function(x) 6.1132*x, colour = 'red', size = 1, alpha = 0.5) +
  annotate("text", x = 0.36, y = 0.95, parse = TRUE, label = "y==6.1132*x", size = 4.5, colour = 'red') + 
  annotate("text", x = 0.47, y = 0.87, parse = TRUE,
           label = expression(italic(r) == 0.786 ~~ italic(p) == 1.09e-42),
           colour = 'red', size = 4.5) +
  guides(colour = guide_legend(title = "Week", ncol = 1)) +
  scale_colour_manual(values = blue2red(18)[3:18]) +
  ylim(0, 1) + xlim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1.5, colour = 'grey80') +
  theme_bw() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(colour = "black", size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        aspect.ratio = 1) + 
  labs(y = "AM fungal relative abundance in metatranscriptome", x = "AM fungal relative abundance in DNA amplicon (rrn uncorrected)")
fig4d_uncorrect_AMF

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4d_pdf <- str_c("fig4d_", "uncorrect_liner_Glomeromycetes_", tm, ".pdf", sep = "")
fig4d_jpg <- str_c("fig4d_", "uncorrect_liner_Glomeromycetes_", tm, ".jpg", sep = "")

ggsave(fig4d_pdf, fig4d_uncorrect_AMF, height = 5.29)
ggsave(fig4d_jpg, fig4d_uncorrect_AMF, height = 5.29)


# next
unadj_otu <- data.frame(t(unadj_otu))

total <- apply(unadj_otu, 1, sum)
unadj_otu.relabu <- data.frame(lapply(unadj_otu, function(x) {  x / total  }) )
unadj_otu.relabu <- data.frame(t(unadj_otu.relabu))

unadj_otu.relabu$ID <- row.names(unadj_otu.relabu)
unadj_otu.relabu <- left_join(unadj_otu.relabu, tax, by = 'ID')

unadj_otu.relabu_class <-aggregate(unadj_otu.relabu[, c(1:196)], by = list(unadj_otu.relabu$Class), sum)
row.names(unadj_otu.relabu_class) <- unadj_otu.relabu_class$Group.1
unadj_otu.relabu_class <- unadj_otu.relabu_class[, -1]
unadj_otu.relabu_class <- data.frame(t(unadj_otu.relabu_class))
unadj_otu.relabu_class$ID_OTU <- row.names(unadj_otu.relabu_class)

trans_fung <- read.csv("./rRNA/Class.csv", header = T, row.names = 1)
trans_fung <- subset(trans_fung, trans_fung$Phylum=="Fungi")

row.names(trans_fung) <- trans_fung$Class

trans_fung <- trans_fung[, -c(1,2)]
trans_fung <- data.frame(t(trans_fung))

total <- apply(trans_fung, 1, sum) 
trans_fung.relabu <- data.frame(lapply(trans_fung, function(x) {  x / total  }) )
trans_fung.relabu$ID1 <- row.names(trans_fung.relabu)

trans_fung.relabu <- left_join(trans_fung.relabu, zg[, c(1, 5)], by = 'ID1')
trans_fung.relabu <- na.omit(trans_fung.relabu)
trans_fung.relabu$ID_new <- paste0("Trans_", trans_fung.relabu$ID_OTU)
row.names(trans_fung.relabu) <- trans_fung.relabu$ID_new
trans_fung.relabu <- trans_fung.relabu[, c(1:58, 61)]

adj_otu.relabu_class$ID_new <- paste0("AdOTU_", adj_otu.relabu_class$ID_OTU)
row.names(adj_otu.relabu_class) <- adj_otu.relabu_class$ID_new
adj_otu.relabu_class <- adj_otu.relabu_class[, c(1:46, 48)]

unadj_otu.relabu_class$ID_new <- paste0("UadOTU_",unadj_otu.relabu_class$ID_OTU)
row.names(unadj_otu.relabu_class) <- unadj_otu.relabu_class$ID_new
unadj_otu.relabu_class <- unadj_otu.relabu_class[, c(1:46, 48)]

all_class <- as.data.frame(plyr::rbind.fill(unadj_otu.relabu_class, adj_otu.relabu_class, trans_fung.relabu))
row.names(all_class) <- all_class$ID_new
group_all <- all_class[, c(46, 47, 47)]
group_all <- cSplit(group_all, "ID_new.1", "_")


colnames(group_all)[4] <- 'ID_OTU'
group_all <- left_join(group_all, zg[, c(5:10)], by = 'ID_OTU')

all_class <- all_class[, -47]
all_class[is.na(all_class)] <- 0

all_class <- all_class[, order(-colSums(all_class))]

group_all$ave_group <- paste0(group_all$ID_new.1_1, "_", group_all$Timepiont, "_", group_all$Treatment1)
row.names(all_class) == group_all$ID_new

all_class.lev <- aggregate(all_class, by = list(group_all$ave_group), mean)
fung.L1 <- all_class.lev[, c(1, 2:11)]
fung.L1 <- melt(fung.L1, id.vars = "Group.1")
fung.L1 <- cSplit(fung.L1, "Group.1", "_")


names(fung.L1) <- c("Fungi", "Relative_Abundance", "Type", "Week", "Treat", "Treat1", "Treat2")

fung.L2 <- all_class.lev[, c(1, 12:ncol(all_class.lev))]
fung.L2 <- melt(fung.L2, id.vars = "Group.1")
fung.L2 <- cSplit(fung.L2, "Group.1", "_")


names(fung.L2) <- c("Fungi", "Relative_Abundance", "Type", "Week", "Treat", "Treat1", "Treat2")
fung.L2$Fungi <- "Other"
fung.bind <- rbind(fung.L1, fung.L2)

fung.bind$Treat <- factor(fung.bind$Treat, levels = c("Control", "Pre", "Post"), labels = c("Control", "Preflowering", "Postflowering"))
fung.bind$Week <- gsub("TP", "", fung.bind$Week, fixed = TRUE)
fung.bind$Week <- as.numeric(fung.bind$Week)
fung.bind$Type <- factor(fung.bind$Type, levels=c("UadOTU", "Trans", "AdOTU"), labels = c("Uncorrected DNA amplicon", "Metatranscriptome", "rrn corrected DNA amplicon"))

fung.bind_1 <- subset(fung.bind, fung.bind$Treat=="Control")

# fig4a, All bar ----------------------
fig4a_All_bar <- 
  ggplot(fung.bind_1, aes(x = factor(Type), y = Relative_Abundance, fill = Fungi)) +
  geom_bar(stat = 'identity', position = "fill") +  
  labs(x = "Week", y = "Relative abundance") +
  facet_grid(. ~ Week) + 
  scale_fill_manual(values = c(
    "tomato","turquoise", "red", "violet", "yellowgreen",
    "peachpuff", "peru", "pink", "plum2", "purple","grey")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 15,face = "bold"), 
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 10, face = "bold.italic"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 65, vjust = 1, hjust = 1)) +
  guides(fill = guide_legend(title = "", byrow = T, nrow = 2, title.position = "left", title.hjust = 0.5))
fig4a_All_bar

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4a_pdf <- str_c("fig4a_", "All_bar_", tm, ".pdf", sep = "")
fig4a_jpg <- str_c("fig4a_", "All_bar_", tm, ".jpg", sep = "")

ggsave(fig4a_pdf, fig4a_All_bar, width = 10.7, height = 7.07)
ggsave(fig4a_jpg, fig4a_All_bar, width = 10.7, height = 7.07)


# fig4b, all box ----------------
dist <- vegdist(decostand(all_class, "hellinger"), method = 'bray')
dist_matrix <- as.matrix(dist)
dist_matrix_long <- melt(dist_matrix)
dist_matrix_long <-cSplit(dist_matrix_long, "Var1", "_")
dist_matrix_long <-cSplit(dist_matrix_long, "Var2", "_")

dist_matrix_long <- subset(dist_matrix_long, dist_matrix_long$Var1_2 == dist_matrix_long$Var2_2)
dist_matrix_long$Group_combind <- paste0(dist_matrix_long$Var1_1, "_", dist_matrix_long$Var2_1)
dist_matrix_long <- subset(dist_matrix_long, dist_matrix_long$Group_combind == "UadOTU_Trans" | dist_matrix_long$Group_combind == "AdOTU_Trans")
colnames(dist_matrix_long)[3] <- 'ID_OTU'

dist_matrix_long <- left_join(dist_matrix_long, group_all[, c(4, 7, 8)], by = 'ID_OTU')

dist_matrix_long$Timepiont <- gsub("TP", "", dist_matrix_long$Timepiont, fixed = TRUE)
dist_matrix_long$Timepiont <- as.numeric(dist_matrix_long$Timepiont)
dist_matrix_long$Group_combind <- factor(dist_matrix_long$Group_combind, 
                                         levels = c("UadOTU_Trans", "AdOTU_Trans"), 
                                         labels = c("Uncorrected","rrn_corrected"))

dist_matrix_long$value <- 1-dist_matrix_long$value
dist_matrix_long_1 <- subset(dist_matrix_long, dist_matrix_long$Treatment1 == "Control")

fig4b_All_box <- 
  ggplot(data = dist_matrix_long_1, aes(x = Group_combind, y = value)) + 
  geom_boxplot() + 
  labs(x = NULL, y = "Similarity to mRNA-seq") +
  geom_signif(comparisons = list(c("Uncorrected", "rrn_corrected")),
              textsize = 3, test = t.test, step_increase = 0.2, map_signif_level = T) +
  scale_y_continuous(limits = c(0.38, 0.90)) +
  facet_grid(. ~Timepiont) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 8, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
fig4b_All_box

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4b_pdf <- str_c("fig4b_", "All_box_", tm, ".pdf", sep = "")
fig4b_jpg <- str_c("fig4b_", "All_box_", tm, ".jpg", sep = "")

ggsave(fig4b_pdf, fig4b_All_box, width = 10.8, height = 3.07)
ggsave(fig4b_jpg, fig4b_All_box, width = 10.8, height = 3.07)


# fig4c-1 ------------------------
lmer_model <- lmer(value ~ Group_combind + (Timepiont|Var2_2), data = dist_matrix_long_1)
anova(lmer_model) %>% tidy()


fig4c_correct <- 
  ggplot(data = dist_matrix_long_1, aes(x = Group_combind, y = value)) + 
  geom_boxplot() + 
  labs(x = NULL, y = "Similarity to mRNA") +
  theme_bw() +
  annotate("text", x = 1.5, y = 0.85,
           label = expression("rRNA adjust:" ~~ F == "86.314," ~~ italic(p) == "9.64e-19"),
           size = 4, colour = "red") +
  theme(axis.title.x=element_blank(), 
        strip.text = element_text(size = 15, face = "bold"),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 8, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        aspect.ratio = 1)
fig4c_correct

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4c_pdf <- str_c("fig4c_", "correct_", tm, ".pdf", sep = "")
fig4c_jpg <- str_c("fig4c_", "correct_", tm, ".jpg", sep = "")

ggsave(fig4c_pdf, fig4c_correct, width = 4.35, height = 4.08)
ggsave(fig4c_jpg, fig4c_correct, width = 4.35, height = 4.08)


# fig4c-2 ---------------
dist_matrix_long_1_uad <- subset(dist_matrix_long_1, dist_matrix_long_1$Group_combind == "Uncorrected")
dist_matrix_long_1_uad <- subset(dist_matrix_long_1_uad, dist_matrix_long_1_uad$value < 0.7)
dist_matrix_long_1_0.25 <- as.data.frame(dist_matrix_long_1[dist_matrix_long_1$ID_OTU%in%dist_matrix_long_1_uad$ID_OTU, ])

lmer_model1 <- lmer(value ~ Group_combind + (Timepiont|Var2_2), data = dist_matrix_long_1_0.25)
anova(lmer_model1) %>% tidy()

fig4c_uncorrect <- 
  ggplot(data = dist_matrix_long_1_0.25, aes(x = Group_combind, y = value)) +
  geom_boxplot() + 
  labs(x = NULL, y = "Similarity to mRNA", title = "Uncorrected similarity < 0.7") +
  theme_bw() + 
  annotate("text", x = 1.5, y = 0.78,
           label = expression("rRNA adjust:" ~~ F == "223.17," ~~ italic(p) == "4.87e-32"),
           size = 4, colour = "red") +
  theme(axis.title.x=element_blank(),
        strip.text = element_text(size = 15, face="bold"),
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(colour = "black", size = 8, face = "bold"),
        legend.text = element_text(colour = "black", size = 8, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        legend.position = "bottom",
        aspect.ratio = 1)
fig4c_uncorrect

tm <- now() %>% str_split_i(pattern = " ", 1)
fig4c_pdf <- str_c("fig4c_", "uncorrect__", tm, ".pdf", sep = "")
fig4c_jpg <- str_c("fig4c_", "uncorrect__", tm, ".jpg", sep = "")

ggsave(fig4c_pdf, fig4c_uncorrect, width = 4.35, height = 4.35)
ggsave(fig4c_jpg, fig4c_uncorrect, width = 4.35, height = 4.35)

# done.

