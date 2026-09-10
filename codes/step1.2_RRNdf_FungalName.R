

# Fungal Names process
# By Qiushi-Li, IM-CAS, 2024.12.24

# packages
library(tidyverse)

library(readxl)
library(writexl)

# read
Fungal_names_raw <- read_delim("./2.database/Fungal_names_all_v2024_10_22.txt", delim = "\t")
# write_xlsx(Fungal_names_raw, "./2.database/Fungal_names_all_v2024_10_22.xlsx")

Fungal_names_sub <- Fungal_names_raw %>% filter(`Name status` == "Current Name")
colnames(Fungal_names_sub) <- str_replace_all(colnames(Fungal_names_sub), pattern = " ", replacement = "_")

Fungal_names_sub

# write_xlsx(Fungal_names_sub, "./database/Fungal_names_all_v2024_10_22_currentname.xlsx")

# annotate taxa level -----
# data
FRRN_rlt <- read_excel("./3.tables/FRRN_rlt_20260707.xlsx", sheet = 1)

FRRN_rlt_F <- FRRN_rlt %>% filter(diff < 150)

FRRN_rlt_F


# restricted
restricted_list <- read_excel("./2.database/restricted_proj_list.xlsx", sheet = 1)
# restricted_list

FRRN_rlt_F1 <- 
  FRRN_rlt_F %>% 
  filter(!project %in% restricted_list$project_id)


# check AMF
# FRRN_rlt_F1 %>% filter(str_detect(project, "AMF")) %>% view()



# read Fungal name
FungalName_current <- read_excel(
  "./2.database/Fungal_names_all_v2024_10_22_currentname.xlsx",
  guess_max = 1000000,
  sheet = 1
)
FungalName_curr_gen <- FungalName_current %>% filter(Rank == "gen.") %>% select(Fungal_name, Classification)


# read JGI F project
FRRN_proj <- read_excel("./2.database/FRRN_project_totalstatus_20260223.xlsx", sheet = 1)
FRRN_proj <- FRRN_proj %>% select(DataSet, Order, Name, Assembly_Length, Genes, Project_ID) %>% drop_na(Project_ID)
FRRN_proj


FRRN_rlt_F2 <-
  FRRN_rlt_F1 %>% left_join(FRRN_proj, by = c("project" = "Project_ID"))



FRRN_rlt_taxa <- FRRN_rlt_F2 %>%
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


dim(FRRN_rlt_F1)
dim(FRRN_rlt_F2)
dim(FRRN_rlt_taxa)

# FRRN_rlt_taxa %>% filter(is.na(Fungal_name)) %>% view()
# maybe needed re-check
# FRRN_rlt_taxa

# FRRN_rlt_taxa %>% filter(str_detect(project, "GaoLab")) %>% view()


# re-check and annotated with Index Fungorum --------------

FRRN_rlt_taxa_spl_with_na <- 
  FRRN_rlt_taxa %>% separate_wider_delim(cols = Classification,
                                         delim = "|",
                                         names = c("kin", "phy", "cla", "ord", "fam"),
                                         too_few = "debug") %>% select(-Assembly_Length, -Genes)


FRRN_rlt_taxa_spl_only_na <- FRRN_rlt_taxa_spl_with_na %>% filter(is.na(Fungal_name)) %>% arrange(Name)
FRRN_rlt_taxa_spl_only_na
# write_xlsx(FRRN_rlt_taxa_spl_only_na, "./3.tables/FRRN_rlt_taxa_spl_na2.xlsx")

# fixed by Index Fungorum --------
FRRN_rlt_taxa_na_fix <- read_excel("./3.tables/FRRN_rlt_taxa_spl_na2_fix.xlsx", sheet = 1)
# view(rCNV_rlt_taxa_na_fix)


dim(FRRN_rlt_taxa_spl_with_na) # 1167
dim(FRRN_rlt_taxa_spl_only_na) # 87
dim(FRRN_rlt_taxa_na_fix) # 87
dim(FRRN_rlt_taxa_spl_with_na %>% filter(!is.na(Fungal_name))) # 1080

FRRN_rlt_taxa_spl <- rbind(
  FRRN_rlt_taxa_spl_with_na %>% filter(!is.na(Fungal_name)),
  FRRN_rlt_taxa_na_fix
) %>% filter(Classification != "no_fungi")


dim(FRRN_rlt_taxa_spl)

# done.


