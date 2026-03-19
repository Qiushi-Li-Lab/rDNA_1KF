

# Fungal Names process
# By Qiushi-Li, IM-CAS, 2024.12.24

# packages
library(tidyverse)

library(readxl)
library(writexl)

# read
Fungal_names_raw <- read_delim("./database/Fungal_names_all_v2024_10_22.txt", delim = "\t")
# write_xlsx(Fungal_names_raw, "./database/Fungal_names_all_v2024_10_22.xlsx")

Fungal_names_sub <- Fungal_names_raw %>% filter(`Name status` == "Current Name")
colnames(Fungal_names_sub) <- str_replace_all(colnames(Fungal_names_sub), pattern = " ", replacement = "_")

Fungal_names_sub

# write_xlsx(Fungal_names_sub, "./database/Fungal_names_all_v2024_10_22_currentname.xlsx")

# annotate taxa level -----
# data
rCNV_rlt <- read_excel("./output/rCNV_rlt_20250223.xlsx", sheet = 1)
rCNV_rlt_Q0 <- rCNV_rlt %>% filter(Quality == "Q0")
rCNV_rlt_Q0$project %>% length()

rCNV_rlt_Q0$project %>% unique() %>% length()

rCNV_rlt_sub <- rCNV_rlt_Q0 %>% filter(diff < 150)
# next path file
# rCNV_rlt_sub

rCNV_rlt_sub %>% filter(str_detect(project, "AMF"))

# read Fungal name
FungalName_current <- read_excel(
  "./database/Fungal_names_all_v2024_10_22_currentname.xlsx",
  guess_max = 1000000,
  sheet = 1
)
FungalName_curr_gen <- FungalName_current %>% filter(Rank == "gen.") %>% select(Fungal_name, Classification)


# read JGI F project
JGI_F_proj <- read_excel("./JGI_F_proj/JGI_totalstatus_20250223.xlsx", sheet = 1)
JGI_F_proj <- JGI_F_proj %>% select(DataSet, Order, Name, Assembly_Length, Genes, Project_ID) %>% drop_na(Project_ID)

# JGI_F_proj %>% filter(DataSet == "AMF_re")



rCNV_rlt_total_tmp <-
  rCNV_rlt_sub %>% left_join(JGI_F_proj, by = c("project" = "Project_ID"))

# rCNV_rlt_total_tmp %>% view()
# rCNV_rlt_total_tmp$Name %>% is.na() %>% table()


rCNV_rlt_taxa <- rCNV_rlt_total_tmp %>%
  mutate(Gen = str_split_i(
    Name,
    pattern = " ",
    i = 1
  )) %>%
  mutate(Gen = str_c("g_", Gen)) %>%
  left_join(
    FungalName_curr_gen %>% mutate(Gen = str_c("g_", Fungal_name)),
    by = "Gen"
  ) %>% distinct(project, .keep_all = T)

dim(rCNV_rlt_sub)
dim(rCNV_rlt_total_tmp)
dim(rCNV_rlt_taxa)


# rCNV_rlt_taxa %>% filter(is.na(Fungal_name)) %>% view()
# maybe needed re-check

# re-check and annoteted with Index Fungrum --------------
rCNV_rlt_taxa_spl_with_na <- 
  rCNV_rlt_taxa %>% separate_wider_delim(cols = Classification,
                                         delim = "|",
                                         names = c("kin", "phy", "cla", "ord", "fam"),
                                         too_few = "debug") %>% select(-Assembly_Length, -Genes)

#write_xlsx(rCNV_rlt_taxa_spl_with_na, "./output/rCNV_rlt_taxa_20250223.xlsx")

rCNV_rlt_taxa_spl_only_na <- rCNV_rlt_taxa_spl_with_na %>% filter(is.na(Fungal_name)) %>% arrange(Name)
rCNV_rlt_taxa_spl_only_na

write_xlsx(rCNV_rlt_taxa_spl_only_na, "./output/rCNV_rlt_taxa_spl_only_na_20250223.xlsx")

# fixed by Index Fungrum --------
rCNV_rlt_taxa_na_fix <- read_excel("./output/rCNV_rlt_taxa_spl_only_na_fix_20250223.xlsx", sheet = 1)
# view(rCNV_rlt_taxa_na_fix)



rCNV_rlt_taxa_spl0 <- rbind(
  rCNV_rlt_taxa_spl_with_na %>% filter(!is.na(Fungal_name)),
  rCNV_rlt_taxa_na_fix
)
dim(rCNV_rlt_taxa_spl0)

# restricted
restricted_list <- read_excel("./output/restricted_proj_list.xlsx", sheet = 1)
restricted_list %>% view()

rCNV_rlt_taxa_spl0 <- rCNV_rlt_taxa_spl0 %>% 
  filter(!project %in% restricted_list$project_id) %>% 
  filter(Classification != "no_fungi")

# dim(rCNV_rlt_taxa_spl)

rCNV_rlt_taxa_spl0


# Q20 -----------------

# data
rCNV_rlt_Q20 <- read_excel("./output/rCNV_rlt_20250223_adj.xlsx", sheet = 1)

rCNV_rlt_Q20 <- rCNV_rlt_Q20 %>% filter(diff < 150)

rCNV_rlt_Q20$project %>% length()
rCNV_rlt_Q20$project %>% unique() %>% length()

rCNV_rlt_Q20 %>% filter(str_detect(project, "AMF"))

rCNV_rlt_total_tmp_Q20 <-
  rCNV_rlt_Q20 %>% left_join(JGI_F_proj, by = c("project" = "Project_ID"))

# rCNV_rlt_total_tmp %>% view()
# rCNV_rlt_total_tmp$Name %>% is.na() %>% table()


rCNV_rlt_taxa_Q20 <- rCNV_rlt_total_tmp_Q20 %>%
  mutate(Gen = str_split_i(
    Name,
    pattern = " ",
    i = 1
  )) %>%
  mutate(Gen = str_c("g_", Gen)) %>%
  left_join(
    FungalName_curr_gen %>% mutate(Gen = str_c("g_", Fungal_name)),
    by = "Gen"
  ) %>% distinct(project, .keep_all = T)

dim(rCNV_rlt_Q20)
dim(rCNV_rlt_total_tmp_Q20)
dim(rCNV_rlt_taxa_Q20)


# re-check and annoteted with Index Fungrum --------------
rCNV_rlt_taxa_spl_with_na_Q20 <- 
  rCNV_rlt_taxa_Q20 %>% separate_wider_delim(cols = Classification,
                                         delim = "|",
                                         names = c("kin", "phy", "cla", "ord", "fam"),
                                         too_few = "debug") %>% select(-Assembly_Length, -Genes)

rCNV_rlt_taxa_spl_only_na_Q20 <- rCNV_rlt_taxa_spl_with_na_Q20 %>% filter(is.na(Fungal_name)) %>% arrange(Name)
rCNV_rlt_taxa_spl_only_na_Q20

write_xlsx(rCNV_rlt_taxa_spl_only_na_Q20, "./output/rCNV_rlt_taxa_spl_only_na_Q20_20250223.xlsx")

# fixed by Index Fungrum --------
rCNV_rlt_taxa_na_fix_Q20 <- read_excel("./output/rCNV_rlt_taxa_spl_only_na_Q20_fix_20250223.xlsx", sheet = 1)
# view(rCNV_rlt_taxa_na_fix)



rCNV_rlt_taxa_spl_Q20 <- rbind(
  rCNV_rlt_taxa_spl_with_na_Q20 %>% filter(!is.na(Fungal_name)),
  rCNV_rlt_taxa_na_fix_Q20
)
dim(rCNV_rlt_taxa_spl_Q20)

rCNV_rlt_taxa_spl_Q20 <- 
  rCNV_rlt_taxa_spl_Q20 %>% 
  filter(!project %in% restricted_list$project_id) %>% 
  filter(Classification != "no_fungi")

rCNV_rlt_taxa_spl_Q20

# switch ------------

#rCNV_rlt_taxa_spl <- rCNV_rlt_taxa_spl0
rCNV_rlt_taxa_spl <- rCNV_rlt_taxa_spl_Q20

# done.


