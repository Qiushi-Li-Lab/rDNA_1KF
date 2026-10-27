



# add analyzed
# packages

library(tidyverse)
library(mgcv)
library(ggsci)

# remotes::install_version("mgcv", "1.9-1")

###### addded to XC_ROOT
source("0.scripts/fix_mgcv_194.R")
XC_ROOT_rDNAm1 <- readRDS("XC_ROOT_rDNAm1.RData")

XC_ROOT_gam <- XC_ROOT_rDNAm1


XC_ROOT_gam1 <- XC_ROOT_gam %>% 
  mutate(
    Site = as.factor(Site),
    Plantsp1 = as.factor(Plantsp1)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )


hist(XC_ROOT_gam1$rDNAm)
hist(XC_ROOT_gam1$rDNAm_log)

# Lat
XC_ROOT_gam_Lat <- gam(rDNAm_log ~ s(Latitude) + s(Site, bs = "re") + s(Plantsp1, bs = "re"), 
                       data = XC_ROOT_gam1, method = "REML")

XC_ROOT_gam1 %>% str()

gam.check(XC_ROOT_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


XC_ROOT_gam_Lat_sum <- summary(XC_ROOT_gam_Lat)
XC_ROOT_gam_Lat_sum


# plot
# predict(XC_ROOT_gam_Lat)

# predict model
XC_ROOT_Lat_model_data <- XC_ROOT_gam_Lat$model

# model plot data
XC_ROOT_Lat_plot <- data.frame(
  Latitude = seq(min(XC_ROOT_Lat_model_data$Latitude), 
                 max(XC_ROOT_Lat_model_data$Latitude),
                 length = 200),
  Site = XC_ROOT_Lat_model_data$Site[1],
  Plantsp1 = XC_ROOT_Lat_model_data$Plantsp1[1]
)


XC_ROOT_Lat_pred <- 
  predict(XC_ROOT_gam_Lat, 
          newdata = XC_ROOT_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)", "s(Plantsp1)"))


ROOT_Lat_plot <- data.frame(
  Latitude = XC_ROOT_Lat_plot$Latitude, 
  fit = XC_ROOT_Lat_pred$fit, 
  se = XC_ROOT_Lat_pred$se.fit
  )


anno_Lat <- data.frame(
  anno_x = range(ROOT_Lat_plot$Latitude)[1] + (range(ROOT_Lat_plot$Latitude)[2] - range(ROOT_Lat_plot$Latitude)[1]) / 2,
  anno_y = Inf
)


# conditional R2
XC_ROOT_gam_Lat_sum$r.sq
# 0.189

# Lat edf
XC_ROOT_gam_Lat_sum$s.table[1, 1]
# 1.955

# Lat F
XC_ROOT_gam_Lat_sum$s.table[1, 3]
# 6.729

# Lat p-value
XC_ROOT_gam_Lat_sum$s.table[1, 4]
# 0.00114



# fig3a, ROOT rDNAm ~ Lat --------------------------
fig_tmp_R_rDNAm_Lat <- 
  ggplot(XC_ROOT_Lat_model_data, aes(Latitude, rDNAm_log)) +
  geom_jitter(colour = "blue", size = 1.8, alpha = 0.15, width = 0.4) +
  # geom_smooth(aes(colour = group), method = "lm", alpha = 0.1) +
  geom_ribbon(data = ROOT_Lat_plot, aes(x = Latitude, y = fit,
                                        ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              inherit.aes = FALSE, alpha = 0.15, fill = "blue") +
  geom_line(data = ROOT_Lat_plot, aes(x = Latitude, y = fit),
            inherit.aes = FALSE, colour = "blue", linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  annotate(geom = "text", x = anno_Lat$anno_x, y = anno_Lat$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.189," ~~ "edf" == "1.955"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat$anno_x, y = anno_Lat$anno_y, 
           label = expression(italic(F) == "6.729," ~~ italic(p) == "0.00114"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 4) +
  labs(x = "Latitude (°)", y = "Community-weighted rDNA copy number\n(log-transformed)") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", face = "bold", size = 15),
    legend.title = element_text(face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.27, 0.93),
    aspect.ratio = 1
  )
fig_tmp_R_rDNAm_Lat

# tm <- now() %>% str_split_i(pattern = " ", 1)
# fig3a_pdf <- str_c("fig3a_", "RS_rDNAm_lat_", tm, ".pdf", sep = "")
# fig3a_jpg <- str_c("fig3a_", "RS_rDNAm_lat_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_3a_pdf <- str_c(fig_path, fig3a_pdf)
# fig_fullpath_3a_jpg <- str_c(fig_path, fig3a_jpg)
# 
# 
# ggsave(fig_fullpath_3a_pdf, fig3a_RS_rDNAm_lat, width = 6.1, height = 5.8)
# ggsave(fig_fullpath_3a_jpg, fig3a_RS_rDNAm_lat, width = 6.1, height = 5.8)




# MAT
XC_ROOT_gam_MAT <- gam(rDNAm_log ~ s(MAT) + s(Site, bs = "re") + s(Plantsp1, bs = "re"), 
                       data = XC_ROOT_gam1, method = "REML")

gam.check(XC_ROOT_gam_MAT)
# concurvity(XC_ROOT_gam_MAT, full = FALSE)


XC_ROOT_gam_MAT_sum <- summary(XC_ROOT_gam_MAT)
XC_ROOT_gam_MAT_sum

# Lat edf
XC_ROOT_gam_MAT_sum$s.table[1, 1]

# Lat F
XC_ROOT_gam_MAT_sum$s.table[1, 3]

# Lat p-value
XC_ROOT_gam_MAT_sum$s.table[1, 4]

# conditional R2
XC_ROOT_gam_MAT_sum$r.sq

# plot
# predict(XC_ROOT_gam_MAT)

# predict model
XC_ROOT_MAT_model_data <- XC_ROOT_gam_MAT$model

# model plot data
XC_ROOT_MAT_plot <- data.frame(
  MAT = seq(min(XC_ROOT_MAT_model_data$MAT), 
            max(XC_ROOT_MAT_model_data$MAT),
            length = 200),
  Site = XC_ROOT_MAT_model_data$Site[1],
  Plantsp1 = XC_ROOT_MAT_model_data$Plantsp1[1]
)


XC_ROOT_MAT_pred <- 
  predict(XC_ROOT_gam_MAT, 
          newdata = XC_ROOT_MAT_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)", "s(Plantsp1)"))


ROOT_MAT_plot <- data.frame(
  MAT = XC_ROOT_MAT_plot$MAT, 
  fit = XC_ROOT_MAT_pred$fit, 
  se = XC_ROOT_MAT_pred$se.fit
)


anno_MAT <- data.frame(
  anno_x = range(ROOT_MAT_plot$MAT)[1] + (range(ROOT_MAT_plot$MAT)[2] - range(ROOT_MAT_plot$MAT)[1]) / 2,
  anno_y = Inf
)


# conditional R2
XC_ROOT_gam_MAT_sum$r.sq
# conditional R2: 0.189

# Lat edf
XC_ROOT_gam_MAT_sum$s.table[1, 1]
# edf: 3.248

# Lat F
XC_ROOT_gam_MAT_sum$s.table[1, 3]
# F: 11.877

# Lat p-value
XC_ROOT_gam_MAT_sum$s.table[1, 4]
# p: <2e-16



# figS6a, ROOT rDNAm ~ MAT --------------------------
fig_S6a_new_XC_ROOT_rDNAm_MAT <- 
  ggplot(XC_ROOT_MAT_model_data, aes(MAT, rDNAm_log)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 0.4) +
  # geom_bin2d() +
  # geom_hex() +
  # geom_smooth(colour = "red", method = "lm") +
  geom_ribbon(data = ROOT_MAT_plot, aes(x = MAT, y = fit,
                                        ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              inherit.aes = FALSE, alpha = 0.15, fill = "red") +
  geom_line(data = ROOT_MAT_plot, aes(x = MAT, y = fit),
            inherit.aes = FALSE, colour = "red", linewidth = 1.1) +
  scale_y_continuous(
                     breaks = seq(3.5, 5, 0.5),
                     labels = floor(exp(seq(3.5, 5, 0.5)))) +
  annotate(geom = "text", x = anno_MAT$anno_x, y = anno_MAT$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.189," ~~ "edf" == "3.248"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_MAT$anno_x, y = anno_MAT$anno_y, 
           label = expression(italic(F) == "11.877," ~~ italic(p) == "< 2e-16"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 4) +
  labs(x = "Mean annual temperature", y = "Community-weighted rDNA copy number\n(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
fig_S6a_new_XC_ROOT_rDNAm_MAT

tm <- now() %>% str_split_i(pattern = " ", 1)
figS6a_new_pdf <- str_c("figS6a_new_", "XC_ROOT_MAT_", tm, ".pdf", sep = "")
figS6a_new_jpg <- str_c("figS6a_new_", "XC_ROOT_MAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S6a_new_pdf <- str_c(fig_path, figS6a_new_pdf)
fig_fullpath_S6a_new_jpg <- str_c(fig_path, figS6a_new_jpg)

ggsave(fig_fullpath_S6a_new_pdf, fig_S6a_new_XC_ROOT_rDNAm_MAT, width = 5.05, height = 4.97)
ggsave(fig_fullpath_S6a_new_jpg, fig_S6a_new_XC_ROOT_rDNAm_MAT, width = 5.05, height = 4.97)




# MAP
XC_ROOT_gam_MAP <- gam(rDNAm_log ~ s(MAP) + s(Site, bs = "re") + s(Plantsp1, bs = "re"), 
                       data = XC_ROOT_gam1, method = "REML")

gam.check(XC_ROOT_gam_MAP)
# concurvity(XC_ROOT_gam_MAT, full = FALSE)


XC_ROOT_gam_MAP_sum <- summary(XC_ROOT_gam_MAP)
XC_ROOT_gam_MAP_sum

# Lat edf
XC_ROOT_gam_MAP_sum$s.table[1, 1]

# Lat F
XC_ROOT_gam_MAP_sum$s.table[1, 3]

# Lat p-value
XC_ROOT_gam_MAP_sum$s.table[1, 4]

# conditional R2
XC_ROOT_gam_MAP_sum$r.sq

# plot
# predict(XC_ROOT_gam_MAT)

# predict model
XC_ROOT_MAP_model_data <- XC_ROOT_gam_MAP$model

# model plot data
XC_ROOT_MAP_plot <- data.frame(
  MAP = seq(min(XC_ROOT_MAP_model_data$MAP), 
            max(XC_ROOT_MAP_model_data$MAP),
            length = 200),
  Site = XC_ROOT_MAP_model_data$Site[1],
  Plantsp1 = XC_ROOT_MAP_model_data$Plantsp1[1]
)


XC_ROOT_MAP_pred <- 
  predict(XC_ROOT_gam_MAP, 
          newdata = XC_ROOT_MAP_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)", "s(Plantsp1)"))


ROOT_MAP_plot <- data.frame(
  MAP = XC_ROOT_MAP_plot$MAP, 
  fit = XC_ROOT_MAP_pred$fit, 
  se = XC_ROOT_MAP_pred$se.fit
)

anno_MAP <- data.frame(
  anno_x = range(ROOT_MAP_plot$MAP)[1] + (range(ROOT_MAP_plot$MAP)[2] - range(ROOT_MAP_plot$MAP)[1]) / 2,
  anno_y = Inf
  )


# conditional R2
XC_ROOT_gam_MAP_sum$r.sq
# conditional R2: 0.188

# Lat edf
XC_ROOT_gam_MAP_sum$s.table[1, 1]
# edf: 2.336

# Lat F
XC_ROOT_gam_MAP_sum$s.table[1, 3]
# F: 10.351

# Lat p-value
XC_ROOT_gam_MAP_sum$s.table[1, 4]
# p: 1.303e-5


# figS6b, RS rDNAm ~ MAP --------------------------
fig_S6b_new_XC_ROOT_rDNAm_MAP <- 
  ggplot(XC_ROOT_MAP_model_data, aes(MAP, rDNAm_log)) +
  geom_jitter(size = 2.8, alpha = 0.4, width = 0.4) +
  # geom_bin2d() +
  # geom_hex() +
  # geom_smooth(colour = "red", method = "lm") +
  geom_ribbon(data = ROOT_MAP_plot, aes(x = MAP, y = fit,
                                        ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              inherit.aes = FALSE, alpha = 0.15, fill = "red") +
  geom_line(data = ROOT_MAP_plot, aes(x = MAP, y = fit),
            inherit.aes = FALSE, colour = "red", linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  annotate(geom = "text", x = anno_MAP$anno_x, y = anno_MAP$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.188," ~~ "edf" == "2.336"),
                     parse = T,
                     colour = "blue",
                     size = 4.5,
                     vjust = 2) +
  annotate(geom = "text", x = anno_MAP$anno_x, y = anno_MAP$anno_y, 
           label = expression(italic(F) == "10.351," ~~ italic(p) == "1.303e-5"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 4) +
  labs(x = "Mean annual precipitation", y = "Community-weighted rDNA copy number\n(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        aspect.ratio = 1)
fig_S6b_new_XC_ROOT_rDNAm_MAP

tm <- now() %>% str_split_i(pattern = " ", 1)
figS6b_new_pdf <- str_c("figS6b_new_", "XC_ROOT_MAP_", tm, ".pdf", sep = "")
figS6b_new_jpg <- str_c("figS6b_new_", "XC_ROOT_MAP_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S6b_new_pdf <- str_c(fig_path, figS6b_new_pdf)
fig_fullpath_S6b_new_jpg <- str_c(fig_path, figS6b_new_jpg)

ggsave(fig_fullpath_S6b_new_pdf, fig_S6b_new_XC_ROOT_rDNAm_MAP, width = 5.05, height = 4.97)
ggsave(fig_fullpath_S6b_new_jpg, fig_S6b_new_XC_ROOT_rDNAm_MAP, width = 5.05, height = 4.97)


# figS6, XC ROOT MAT and MAP 
figS6_new_XC_ROOT_MATMAP <- fig_S6a_new_XC_ROOT_rDNAm_MAT + fig_S6b_new_XC_ROOT_rDNAm_MAP

tm <- now() %>% str_split_i(pattern = " ", 1)
figS6_new_pdf <- str_c("figS6_new_", "XC_ROOT_MATMAP_", tm, ".pdf", sep = "")
figS6_new_jpg <- str_c("figS6_new_", "XC_ROOT_MATMAP_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S6_new_pdf <- str_c(fig_path, figS6_new_pdf)
fig_fullpath_S6_new_jpg <- str_c(fig_path, figS6_new_jpg)

ggsave(fig_fullpath_S6_new_pdf, figS6_new_XC_ROOT_MATMAP, width = 9.47, height = 4.7)
ggsave(fig_fullpath_S6_new_jpg, figS6_new_XC_ROOT_MATMAP, width = 9.47, height = 4.7)



##### Zheng soil #####
Zsoil_cfb_rDNAm1 <- readRDS("Zsoil_cfb_rDNAm1.RData")
Zsoil_gam <- Zsoil_cfb_rDNAm1

# Zsoil_cfb_ECM_rDNAm1 %>% select(sample_id, Latitude1, Latitude) %>% view()

Zsoil_gam1 <- Zsoil_gam %>% 
  mutate(
    Site = as.factor(Site)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )

# hist(Zsoil_gam1$rDNAm)
# hist(Zsoil_gam1$rDNAm_log)

# Lat using GAMMs
Zsoil_gam_Lat <- gam(rDNAm_log ~ s(Latitude1) + s(Site, bs = "re"), 
                       data = Zsoil_gam1, method = "REML")

# gam.check(Zsoil_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_gam_Lat_sum <- summary(Zsoil_gam_Lat)
Zsoil_gam_Lat_sum
# edf = 1, the model needed to adjust


# Lat using LMM
Zsoil_gam_Lat <- gam(rDNAm_log ~ Latitude1 + s(Site, bs = "re"), 
                     data = Zsoil_gam1, method = "REML")

# gam.check(Zsoil_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_gam_Lat_sum <- summary(Zsoil_gam_Lat)
Zsoil_gam_Lat_sum


# plot
# predict(XC_ROOT_gam_Lat)

# predict model
Zsoil_gam_Lat_model_data <- Zsoil_gam_Lat$model

# model plot data
Zsoil_Lat_plot <- data.frame(
  Latitude1 = seq(min(Zsoil_gam_Lat_model_data$Latitude1), 
                 max(Zsoil_gam_Lat_model_data$Latitude1),
                 length = 200),
  Site = Zsoil_gam_Lat_model_data$Site[1]
)


Zsoil_Lat_pred <- 
  predict(Zsoil_gam_Lat, 
          newdata = Zsoil_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)"))


Soil_Lat_plot <- data.frame(
  Latitude1 = Zsoil_Lat_plot$Latitude1, 
  fit = Zsoil_Lat_pred$fit, 
  se = Zsoil_Lat_pred$se.fit
)


anno_Lat_Soil <- data.frame(
  anno_x = range(Soil_Lat_plot$Latitude1)[1] + (range(Soil_Lat_plot$Latitude1)[2] - range(Soil_Lat_plot$Latitude1)[1]) / 2,
  anno_y = Inf
)


# conditional R2
Zsoil_gam_Lat_sum$r.sq
# 0.0918

# Lat t
Zsoil_gam_Lat_sum$p.table[2, 3]
# 0.999

# Lat p-value
Zsoil_gam_Lat_sum$p.table[2, 4]
# 0.319



# fig3a
XC_ROOT_Lat_model_data
Zsoil_gam_Lat_model_data


ROOT_data <- XC_ROOT_Lat_model_data %>% select(rDNAm_log, Latitude) %>% 
  mutate(
    group = "Root"
  )

SOIL_data <- Zsoil_gam_Lat_model_data %>% select(rDNAm_log, Latitude1) %>% 
  mutate(
    group = "Soil"
  )
colnames(SOIL_data)[2] <- "Latitude"


RS_plot_data <- rbind(
  ROOT_data, SOIL_data
)

anno_Lat_RS <- data.frame(
  anno_x = range(RS_plot_data$Latitude)[1] + (range(RS_plot_data$Latitude)[2] - range(RS_plot_data$Latitude)[1]) / 2,
  anno_y = Inf
)

R_Lat_plot <- ROOT_Lat_plot %>% mutate(group = "Root")
S_Lat_plot <- Soil_Lat_plot %>% mutate(group = "Soil") 
colnames(S_Lat_plot)[1] <- "Latitude"


RS_Lat_plot <- rbind(
  R_Lat_plot, S_Lat_plot
)



# fig3a, Soil rDNAm ~ Lat --------------------------
fig_3a_RS_rDNAm_Lat <- 
  ggplot(RS_plot_data, aes(Latitude, rDNAm_log)) +
  geom_jitter(aes(colour = group), size = 1.8, alpha = 0.15, width = 0.4) +
  geom_ribbon(data = RS_Lat_plot, aes(x = Latitude, y = fit, fill = group,
                                        ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = RS_Lat_plot, aes(x = Latitude, y = fit, colour = group),
            inherit.aes = FALSE, linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  annotate(geom = "text", x = anno_Lat_Soil$anno_x, y = anno_Lat$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.0918"),
           parse = T,
           colour = "red",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat_Soil$anno_x, y = anno_Lat$anno_y, 
           label = expression(italic(t) == "0.999," ~~ italic(p) == "0.319"),
           parse = T,
           colour = "red",
           size = 4.5,
           vjust = 4) +
  annotate(geom = "text", x = anno_Lat$anno_x, y = anno_Lat$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.189," ~~ "edf" == "1.955"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 5) +
  annotate(geom = "text", x = anno_Lat$anno_x, y = anno_Lat$anno_y, 
           label = expression(italic(F) == "6.729," ~~ italic(p) == "0.00114"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 8) +
  guides(colour = guide_legend(title = NULL, nrow = 1, override.aes = list(size = 5)),
         fill = guide_legend(title = NULL, nrow = 1)) +
  scale_colour_manual(values = c("blue",
                                 "red")) +
  scale_fill_manual(values = c("blue",
                                 "red")) +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number\n(log-transformed)") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", face = "bold", size = 15),
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.background = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.1),
    aspect.ratio = 1
  )
fig_3a_RS_rDNAm_Lat


tm <- now() %>% str_split_i(pattern = " ", 1)
fig3a_new_pdf <- str_c("fig3a_new_", "RS_rDNAm_lat_", tm, ".pdf", sep = "")
fig3a_new_jpg <- str_c("fig3a_new_", "RS_rDNAm_lat_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_3a_new_pdf <- str_c(fig_path, fig3a_new_pdf)
fig_fullpath_3a_new_jpg <- str_c(fig_path, fig3a_new_jpg)


ggsave(fig_fullpath_3a_new_pdf, fig_3a_RS_rDNAm_Lat, width = 6.1, height = 5.8)
ggsave(fig_fullpath_3a_new_jpg, fig_3a_RS_rDNAm_Lat, width = 6.1, height = 5.8)



##### EMF ####
# XC
XC_ROOT_EMF_rDNAm1 <- readRDS("XC_ROOT_EMF_rDNAm1.RData")


XC_EMF_gam <- XC_ROOT_EMF_rDNAm1


XC_EMF_gam1 <- XC_EMF_gam %>% 
  mutate(
    Site = as.factor(Site),
    Plantsp1 = as.factor(Plantsp1),
    Family = as.factor(Family)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )


hist(XC_EMF_gam1$rDNAm)
hist(XC_EMF_gam1$rDNAm_log)

# EMF Lat
XC_EMF_gam_Lat <- gam(rDNAm_log ~ s(Latitude) + s(Site, bs = "re") + s(Family, bs = "re"), 
                       data = XC_EMF_gam1, method = "REML")

gam.check(XC_EMF_gam_Lat)
# concurvity(XC_EMF_gam_Lat, full = FALSE)


XC_EMF_gam_Lat_sum <- summary(XC_EMF_gam_Lat)
XC_EMF_gam_Lat_sum


# plot
# predict(XC_EMF_gam_Lat)

# predict model
XC_EMF_Lat_model_data <- XC_EMF_gam_Lat$model

# model plot data
XC_EMF_Lat_plot <- data.frame(
  Latitude = seq(min(XC_EMF_Lat_model_data$Latitude), 
                 max(XC_EMF_Lat_model_data$Latitude),
                 length = 200),
  Site = XC_EMF_Lat_model_data$Site[1],
  Family = XC_EMF_Lat_model_data$Family[1]
)


XC_EMF_Lat_pred <- 
  predict(XC_EMF_gam_Lat, 
          newdata = XC_EMF_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)", "s(Family)"))


EMF_Lat_plot <- data.frame(
  Latitude = XC_EMF_Lat_plot$Latitude, 
  fit = XC_EMF_Lat_pred$fit, 
  se = XC_EMF_Lat_pred$se.fit
)


anno_Lat_R_EMF <- data.frame(
  anno_x = range(EMF_Lat_plot$Latitude)[1] + (range(EMF_Lat_plot$Latitude)[2] - range(EMF_Lat_plot$Latitude)[1]) / 2,
  anno_y = Inf
)


# conditional R2
XC_EMF_gam_Lat_sum$r.sq
# 0.183

# Lat edf
XC_EMF_gam_Lat_sum$s.table[1, 1]
# 1.802

# Lat F
XC_EMF_gam_Lat_sum$s.table[1, 3]
# 5.351

# Lat p-value
XC_EMF_gam_Lat_sum$s.table[1, 4]
# 0.00382

# XC_EMF_Lat_model_data

fig_S7a_XC_ROOT_EMF_rDNAm_LAT_full <- 
  ggplot(XC_EMF_Lat_model_data, aes(Latitude, rDNAm_log, colour = Family)) +
  geom_jitter(size = 2.5, alpha = 0.5, width = 0.4) +
  geom_ribbon(data = EMF_Lat_plot, aes(x = Latitude, y = fit,
                                      ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              inherit.aes = FALSE, alpha = 0.15, fill = "red") +
  geom_line(data = EMF_Lat_plot, aes(x = Latitude, y = fit),
            inherit.aes = FALSE, linewidth = 1.1, colour = "red") +
  scale_y_continuous(limits = c(2.8, 5.4),
    breaks = seq(3, 5, 0.5),
    labels = floor(exp(seq(3, 5, 0.5)))) +
  annotate(geom = "text", x = anno_Lat_R_EMF$anno_x, y = anno_Lat_R_EMF$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.183," ~~ "edf" == "1.802"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat_R_EMF$anno_x, y = anno_Lat_R_EMF$anno_y, 
           label = expression(italic(F) == "5.351," ~~ italic(p) == "0.00382"),
           parse = T,
           colour = "blue",
           size = 4.5,
           vjust = 4) +
  scale_colour_nejm() +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal")) +
  labs(y = "Community-weighted rDNA copy number\nof ectomycorrhizal fungi(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold.italic"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.1))
fig_S7a_XC_ROOT_EMF_rDNAm_LAT_full

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7a_new_pdf <- str_c("figS7a_new_", "XCroot_EMF_LAT_", tm, ".pdf", sep = "")
figS7a_new_jpg <- str_c("figS7a_new_", "XCroot_EMF_LAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7a_new_pdf <- str_c(fig_path, figS7a_new_pdf)
fig_fullpath_S7a_new_jpg <- str_c(fig_path, figS7a_new_jpg)

ggsave(fig_fullpath_S7a_new_pdf, fig_S7a_XC_ROOT_EMF_rDNAm_LAT_full, width = 5.11, height = 4.67)
ggsave(fig_fullpath_S7a_new_jpg, fig_S7a_XC_ROOT_EMF_rDNAm_LAT_full, width = 5.11, height = 4.67)




# Zheng ECM Lat
Zsoil_cfb_ECM_rDNAm1 <- readRDS("Zsoil_cfb_ECM_rDNAm1.RData")
Zsoil_ECM_gam <- Zsoil_cfb_ECM_rDNAm1

Zsoil_ECM_gam1 <- Zsoil_ECM_gam %>% 
  mutate(
    Site = as.factor(Site)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )

# hist(Zsoil_ECM_gam1$rDNAm)
# hist(Zsoil_ECM_gam1$rDNAm_log)

# Lat
Zsoil_ECM_gam_Lat <- gam(rDNAm_log ~ s(Latitude1) + s(Site, bs = "re"), 
                     data = Zsoil_ECM_gam1, method = "REML")

# gam.check(Zsoil_ECM_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_ECM_gam_Lat_sum <- summary(Zsoil_ECM_gam_Lat)
Zsoil_ECM_gam_Lat_sum



# plot
# predict(XC_ROOT_gam_Lat)

# predict model
Zsoil_ECM_gam_Lat_model_data <- Zsoil_ECM_gam_Lat$model

# model plot data
Zsoil_ECM_Lat_plot <- data.frame(
  Latitude1 = seq(min(Zsoil_ECM_gam_Lat_model_data$Latitude1), 
                  max(Zsoil_ECM_gam_Lat_model_data$Latitude1),
                  length = 200),
  Site = Zsoil_ECM_gam_Lat_model_data$Site[1]
)


Zsoil_ECM_Lat_pred <- 
  predict(Zsoil_ECM_gam_Lat, 
          newdata = Zsoil_ECM_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)"))


ECM_Lat_plot_S <- data.frame(
  Latitude1 = Zsoil_ECM_Lat_plot$Latitude1, 
  fit = Zsoil_ECM_Lat_pred$fit, 
  se = Zsoil_ECM_Lat_pred$se.fit
)




# Zheng PP Lat
Zsoil_cfb_PP_rDNAm1 <- readRDS("Zsoil_cfb_PP_rDNAm1.RData")

Zsoil_PP_gam <- Zsoil_cfb_PP_rDNAm1


Zsoil_PP_gam1 <- Zsoil_PP_gam %>% 
  mutate(
    Site = as.factor(Site)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )

# hist(Zsoil_PP_gam1$rDNAm)
# hist(Zsoil_PP_gam1$rDNAm_log)

# Lat
Zsoil_PP_gam_Lat <- gam(rDNAm_log ~ s(Latitude1) + s(Site, bs = "re"), 
                         data = Zsoil_PP_gam1, method = "REML")

# gam.check(Zsoil_PP_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_PP_gam_Lat_sum <- summary(Zsoil_PP_gam_Lat)
Zsoil_PP_gam_Lat_sum
# edf = 1, should use LMM


Zsoil_PP_gam_Lat <- gam(rDNAm_log ~ Latitude1 + s(Site, bs = "re"), 
                        data = Zsoil_PP_gam1, method = "REML")

# gam.check(Zsoil_PP_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_PP_gam_Lat_sum <- summary(Zsoil_PP_gam_Lat)
Zsoil_PP_gam_Lat_sum

# plot
# predict(XC_ROOT_gam_Lat)

# predict model
Zsoil_PP_gam_Lat_model_data <- Zsoil_PP_gam_Lat$model

# model plot data
Zsoil_PP_Lat_plot <- data.frame(
  Latitude1 = seq(min(Zsoil_PP_gam_Lat_model_data$Latitude1), 
                  max(Zsoil_PP_gam_Lat_model_data$Latitude1),
                  length = 200),
  Site = Zsoil_PP_gam_Lat_model_data$Site[1]
)


Zsoil_PP_Lat_pred <- 
  predict(Zsoil_PP_gam_Lat, 
          newdata = Zsoil_PP_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)"))


PP_Lat_plot_S <- data.frame(
  Latitude1 = Zsoil_PP_Lat_plot$Latitude1, 
  fit = Zsoil_PP_Lat_pred$fit, 
  se = Zsoil_PP_Lat_pred$se.fit
)





# Sap Lat
Zsoil_cfb_Sap_rDNAm1 <- readRDS("Zsoil_cfb_Sap_rDNAm1.RData")

Zsoil_Sap_gam <- Zsoil_cfb_Sap_rDNAm1

Zsoil_Sap_gam1 <- Zsoil_Sap_gam %>% 
  mutate(
    Site = as.factor(Site)
  ) %>% 
  mutate(
    rDNAm_log = log(rDNAm)
  )

# hist(Zsoil_Sap_gam1$rDNAm)
# hist(Zsoil_Sap_gam1$rDNAm_log)

# Lat
Zsoil_Sap_gam_Lat <- gam(rDNAm_log ~ s(Latitude1) + s(Site, bs = "re"), 
                        data = Zsoil_Sap_gam1, method = "REML")

# gam.check(Zsoil_Sap_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_Sap_gam_Lat_sum <- summary(Zsoil_Sap_gam_Lat)
Zsoil_Sap_gam_Lat_sum
# edf = 1, use LMM


Zsoil_Sap_gam_Lat <- gam(rDNAm_log ~ Latitude1 + s(Site, bs = "re"), 
                         data = Zsoil_Sap_gam1, method = "REML")

# gam.check(Zsoil_Sap_gam_Lat)
# concurvity(XC_ROOT_gam_Lat, full = FALSE)


Zsoil_Sap_gam_Lat_sum <- summary(Zsoil_Sap_gam_Lat)
Zsoil_Sap_gam_Lat_sum

# plot
# predict(XC_ROOT_gam_Lat)

# predict model
Zsoil_Sap_gam_Lat_model_data <- Zsoil_Sap_gam_Lat$model

# model plot data
Zsoil_Sap_Lat_plot <- data.frame(
  Latitude1 = seq(min(Zsoil_Sap_gam_Lat_model_data$Latitude1), 
                  max(Zsoil_Sap_gam_Lat_model_data$Latitude1),
                  length = 200),
  Site = Zsoil_Sap_gam_Lat_model_data$Site[1]
)


Zsoil_Sap_Lat_pred <- 
  predict(Zsoil_Sap_gam_Lat, 
          newdata = Zsoil_Sap_Lat_plot, type = "link",
          se.fit = T,
          exclude = c("s(Site)"))


Sap_Lat_plot_S <- data.frame(
  Latitude1 = Zsoil_Sap_Lat_plot$Latitude1, 
  fit = Zsoil_Sap_Lat_pred$fit, 
  se = Zsoil_Sap_Lat_pred$se.fit
)


# Zheng Soil plot
Z_ECM_data <- Zsoil_ECM_gam_Lat_model_data %>% select(rDNAm_log, Latitude1) %>% 
  mutate(
    group = "Ectomycorrhizal fungi"
  )
colnames(Z_ECM_data)[2] <- "Latitude"

Z_PP_data <- Zsoil_PP_gam_Lat_model_data %>% select(rDNAm_log, Latitude1) %>% 
  mutate(
    group = "Plant pathogen"
  )
colnames(Z_PP_data)[2] <- "Latitude"

Z_Sap_data <- Zsoil_Sap_gam_Lat_model_data %>% select(rDNAm_log, Latitude1) %>% 
  mutate(
    group = "Saprotroph fungi"
  )
colnames(Z_Sap_data)[2] <- "Latitude"



EPS_plot_data <- rbind(
  Z_ECM_data, Z_PP_data, Z_Sap_data
)


anno_Lat_EPS <- data.frame(
  anno_x = range(EPS_plot_data$Latitude)[1] + (range(EPS_plot_data$Latitude)[2] - range(EPS_plot_data$Latitude)[1]) / 2,
  anno_y = Inf
)





EPS_Lat_plot <- rbind(
  ECM_Lat_plot_S %>% mutate(group = "Ectomycorrhizal fungi"),
  PP_Lat_plot_S %>% mutate(group = "Plant pathogen"),
  Sap_Lat_plot_S %>% mutate(group = "Saprotroph fungi")
)
colnames(EPS_Lat_plot)[1] <- "Latitude"


### ECM
# conditional R2
Zsoil_ECM_gam_Lat_sum$r.sq
# 0.222

# Lat edf
Zsoil_ECM_gam_Lat_sum$s.table[1, 1]
# 3.373

# Lat F
Zsoil_ECM_gam_Lat_sum$s.table[1, 3]
# 4.364

# Lat p-value
Zsoil_ECM_gam_Lat_sum$s.table[1, 4]
# 0.00252



### PP
# conditional R2
Zsoil_PP_gam_Lat_sum$r.sq
# 0.152

# Lat t
Zsoil_PP_gam_Lat_sum$p.table[2, 3]
# -1.850

# Lat p-value
Zsoil_PP_gam_Lat_sum$p.table[2, 4]
# 0.0662


### Sap
# conditional R2
Zsoil_Sap_gam_Lat_sum$r.sq
# 0.271

# Lat F
Zsoil_Sap_gam_Lat_sum$p.table[2, 3]
# -0.707

# Lat p-value
Zsoil_Sap_gam_Lat_sum$p.table[2, 4]
# 0.480


# plot
figS7b_EPS_rDNAm_Lat_E <- 
  ggplot(EPS_plot_data %>% filter(group == "Ectomycorrhizal fungi"), aes(Latitude, rDNAm_log)) +
  geom_jitter(colour = "#3B4992FF", size = 2.5, alpha = 0.4, width = 0.4) +
  geom_ribbon(data = EPS_Lat_plot %>% filter(group == "Ectomycorrhizal fungi"), aes(x = Latitude, y = fit,
                                      ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), fill = "#3B4992FF",
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = EPS_Lat_plot %>% filter(group == "Ectomycorrhizal fungi"), aes(x = Latitude, y = fit), colour = "#3B4992FF",
            inherit.aes = FALSE, linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
           label = expression("Conditional" ~~ italic(R^2) == "0.222," ~~ "edf" == "3.373"),
           parse = T,
           colour = "#3B4992FF",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
           label = expression(italic(F) == "4.364," ~~ italic(p) == "0.00252"),
           parse = T,
           colour = "#3B4992FF",
           size = 4.5,
           vjust = 4) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
  #          label = expression("Conditional" ~~ italic(R^2) == "0.152," ~~ "edf" == "1.000," ~~ italic(F) == "3.422," ~~ italic(p) == "0.0662"),
  #          parse = T,
  #          colour = "#EE0000FF",
  #          size = 3.5,
  #          vjust = 3.5) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
  #          label = expression("Conditional" ~~ italic(R^2) == "0.271," ~~ "edf" == "1.000," ~~ italic(F) == "0.500," ~~ italic(p) == "0.480"),
  #          parse = T,
  #          colour = "#008B45FF",
  #          size = 3.5,
  #          vjust = 5) +
  guides(colour = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal"),
         fill = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal")) +
  scale_colour_aaas() +
  scale_fill_aaas() +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number of\nectomycorrhizal fungi(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.1))
figS7b_EPS_rDNAm_Lat_E


tm <- now() %>% str_split_i(pattern = " ", 1)
figS7b_new_pdf <- str_c("figS7b_new_", "Zsoil_LAT_E_", tm, ".pdf", sep = "")
figS7b_new_jpg <- str_c("figS7b_new_", "Zsoil_LAT_E_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7b_new_pdf <- str_c(fig_path, figS7b_new_pdf)
fig_fullpath_S7b_new_jpg <- str_c(fig_path, figS7b_new_jpg)

ggsave(fig_fullpath_S7b_new_pdf, figS7b_EPS_rDNAm_Lat_E, width = 5.11, height = 4.67)
ggsave(fig_fullpath_S7b_new_jpg, figS7b_EPS_rDNAm_Lat_E, width = 5.11, height = 4.67)


figS7c_EPS_rDNAm_Lat_P <- 
  ggplot(EPS_plot_data %>% filter(group == "Plant pathogen"), aes(Latitude, rDNAm_log)) +
  geom_jitter(colour = "#EE0000FF", size = 2.5, alpha = 0.4, width = 0.4) +
  geom_ribbon(data = EPS_Lat_plot %>% filter(group == "Plant pathogen"), aes(x = Latitude, y = fit,
                                                                                    ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), fill = "#EE0000FF",
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = EPS_Lat_plot %>% filter(group == "Plant pathogen"), aes(x = Latitude, y = fit), colour = "#EE0000FF",
            inherit.aes = FALSE, linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
  #          label = expression("Conditional" ~~ italic(R^2) == "0.222," ~~ "edf" == "3.373," ~~ italic(F) == "4.364," ~~ italic(p) == "0.00252"),
  #          parse = T,
  #          colour = "#3B4992FF",
  #          size = 3.5,
  #          vjust = 2) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y,
           label = expression("Conditional" ~~ italic(R^2) == "0.152"),
           parse = T,
           colour = "#EE0000FF",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
           label = expression(italic(t) == "-1.850," ~~ italic(p) == "0.0662"),
           parse = T,
           colour = "#EE0000FF",
           size = 4.5,
           vjust = 4) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
  #          label = expression("Conditional" ~~ italic(R^2) == "0.271," ~~ "edf" == "1.000," ~~ italic(F) == "0.500," ~~ italic(p) == "0.480"),
  #          parse = T,
  #          colour = "#008B45FF",
  #          size = 3.5,
  #          vjust = 5) +
  guides(colour = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal"),
         fill = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                             nrow = 1,
                             direction = "horizontal")) +
  scale_colour_aaas() +
  scale_fill_aaas() +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number of\nPlant pathogen(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.1))
figS7c_EPS_rDNAm_Lat_P

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7c_new_pdf <- str_c("figS7c_new_", "Zsoil_LAT_P_", tm, ".pdf", sep = "")
figS7c_new_jpg <- str_c("figS7c_new_", "Zsoil_LAT_P_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7c_new_pdf <- str_c(fig_path, figS7c_new_pdf)
fig_fullpath_S7c_new_jpg <- str_c(fig_path, figS7c_new_jpg)

ggsave(fig_fullpath_S7c_new_pdf, figS7c_EPS_rDNAm_Lat_P, width = 5.11, height = 4.67)
ggsave(fig_fullpath_S7c_new_jpg, figS7c_EPS_rDNAm_Lat_P, width = 5.11, height = 4.67)


figS7d_EPS_rDNAm_Lat_S <- 
  ggplot(EPS_plot_data %>% filter(group == "Saprotroph fungi"), aes(Latitude, rDNAm_log)) +
  geom_jitter(colour = "#008B45FF", size = 2.5, alpha = 0.4, width = 0.4) +
  geom_ribbon(data = EPS_Lat_plot %>% filter(group == "Saprotroph fungi"), aes(x = Latitude, y = fit,
                                                                             ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), fill = "#008B45FF",
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = EPS_Lat_plot %>% filter(group == "Saprotroph fungi"), aes(x = Latitude, y = fit), colour = "#008B45FF",
            inherit.aes = FALSE, linewidth = 1.1) +
  scale_y_continuous(
    breaks = seq(3.5, 5, 0.5),
    labels = floor(exp(seq(3.5, 5, 0.5)))) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
  #          label = expression("Conditional" ~~ italic(R^2) == "0.222," ~~ "edf" == "3.373," ~~ italic(F) == "4.364," ~~ italic(p) == "0.00252"),
  #          parse = T,
  #          colour = "#3B4992FF",
  #          size = 3.5,
  #          vjust = 2) +
  # annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y,
  #          label = expression("Conditional" ~~ italic(R^2) == "0.152," ~~ italic(t) == "-1.850," ~~ italic(p) == "0.0662"),
  #          parse = T,
  #          colour = "#EE0000FF",
  #          size = 3.5,
  #          vjust = 2) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y,
           label = expression("Conditional" ~~ italic(R^2) == "0.271"),
           parse = T,
           colour = "#008B45FF",
           size = 4.5,
           vjust = 2) +
  annotate(geom = "text", x = anno_Lat_EPS$anno_x, y = anno_Lat_EPS$anno_y, 
           label = expression(italic(t) == "-0.707," ~~ italic(p) == "0.480"),
           parse = T,
           colour = "#008B45FF",
           size = 4.5,
           vjust = 4) +
  guides(colour = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                               nrow = 1,
                               direction = "horizontal"),
         fill = guide_legend(NULL, override.aes = list(size = 5, alpha = 1),
                             nrow = 1,
                             direction = "horizontal")) +
  scale_colour_aaas() +
  scale_fill_aaas() +
  labs(x = "Latitude", y = "Community-weighted rDNA copy number of\nSaprotroph fungi(log-transformed)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", colour = "black", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold"),
        # legend.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.1))
figS7d_EPS_rDNAm_Lat_S

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7d_new_pdf <- str_c("figS7d_new_", "Zsoil_LAT_S_", tm, ".pdf", sep = "")
figS7d_new_jpg <- str_c("figS7d_new_", "Zsoil_LAT_S_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7d_new_pdf <- str_c(fig_path, figS7d_new_pdf)
fig_fullpath_S7d_new_jpg <- str_c(fig_path, figS7d_new_jpg)

ggsave(fig_fullpath_S7d_new_pdf, figS7d_EPS_rDNAm_Lat_S, width = 5.11, height = 4.67)
ggsave(fig_fullpath_S7d_new_jpg, figS7d_EPS_rDNAm_Lat_S, width = 5.11, height = 4.67)

# tm <- now() %>% str_split_i(pattern = " ", 1)
# figS7b_new_pdf <- str_c("figS7b_new_", "Zsoil_EPS_LAT_", tm, ".pdf", sep = "")
# figS7b_new_jpg <- str_c("figS7b_new_", "Zsoil_EPS_LAT_", tm, ".jpg", sep = "")
# 
# # figs were saved in 4.figs
# fig_path <- "./4.figs/"
# 
# fig_fullpath_S7b_new_pdf <- str_c(fig_path, figS7b_new_pdf)
# fig_fullpath_S7b_new_jpg <- str_c(fig_path, figS7b_new_jpg)
# 
# ggsave(fig_fullpath_S7b_new_pdf, fig_tmp_EPS_rDNAm_Lat, width = 5.36, height = 5.05)
# ggsave(fig_fullpath_S7b_new_jpg, fig_tmp_EPS_rDNAm_Lat, width = 5.36, height = 5.05)




# figS7
figS7_new <- fig_S7a_XC_ROOT_EMF_rDNAm_LAT_full + figS7b_EPS_rDNAm_Lat_E +
  figS7c_EPS_rDNAm_Lat_P + figS7d_EPS_rDNAm_Lat_S + plot_layout(ncol =  2)

figS7_new

tm <- now() %>% str_split_i(pattern = " ", 1)
figS7_new_pdf <- str_c("figS7_new_", "XC_ROOT_Zsoil_LAT_", tm, ".pdf", sep = "")
figS7_new_jpg <- str_c("figS7_new_", "XC_ROOT_Zsoil_LAT_", tm, ".jpg", sep = "")

# figs were saved in 4.figs
fig_path <- "./4.figs/"

fig_fullpath_S7_new_pdf <- str_c(fig_path, figS7_new_pdf)
fig_fullpath_S7_new_jpg <- str_c(fig_path, figS7_new_jpg)

ggsave(fig_fullpath_S7_new_pdf, figS7_new, width = 11, height = 10.06)
ggsave(fig_fullpath_S7_new_jpg, figS7_new, width = 11, height = 10.06)
# Saving 10.4 x 5.06 in image



