
# this script could process RCNV_pipeline.sh results into a table-like sheets
# by Qiushi-Li, IMCAS
# 2024.10.11

# Rdata tmp


# the package we needed to loading
library(tidyverse) # data science

# read & writexl
library(readxl)
library(writexl)

# multi session calcu
library(furrr)
library(future)
plan(multisession)


# FRRN fungi path & rCNV project name
FRRN_data_path <- dir("./1.data/Fungal_RRN/", full.names = T)

str_detect(FRRN_data_path, "none", negate = F) %>% table()
str_detect(FRRN_data_path, "sobig", negate = F) %>% table()
str_detect(FRRN_data_path, "noITS", negate = F) %>% table()


FRRN_data_path_F <- FRRN_data_path[str_detect(FRRN_data_path, "none", negate = T)]
FRRN_data_path_F <- FRRN_data_path_F[str_detect(FRRN_data_path_F, "sobig", negate = T)]
FRRN_data_path_F <- FRRN_data_path_F[str_detect(FRRN_data_path_F, "noITS", negate = T)]

# FRRN_rlt_path

FRRN_proj_name <- str_split(FRRN_data_path_F, pattern = "\\.", n = 4, simplify = T)[, 4]

# view(FRRN_proj_name)
# view(FRRN_rlt_path)
# view(FRRN_proj_name)


# clean & extract
# read the data
# tmp path

# list.files(FRRN_rlt_path[1])

FRRN_rlt_process <- function(project_path, project_name) {
  
  # filtered rows
  need_row <- c("SCG depth average=",
                "estimated rDNA CN (ITS only)=",
                "estimated rDNA CN (ITS / LSU)=",
                "% diff ITS / LSU=")
  
  # project file path
  project_rlt_path <- str_c(project_path, "/CN_rlt")
  project_files <- list.files(project_rlt_path)
  
  # project file rlt nums
  project_files_num <- length(project_files)
  
  if(project_files_num > 1) {
    
    # read_Q0...
    Q0_path <- str_c(project_rlt_path, "/", project_files[1])
    project_tmp_Q0 <- read_delim(Q0_path, col_names = F, delim = "\t")
    project_tmp_Q0_rlt <- project_tmp_Q0 %>% mutate(Quality = "Q0", Confidence = "Possible")
    
    # read_Q20
    Q20_path <- str_c(project_rlt_path, "/", project_files[2])
    project_tmp_Q20 <- read_delim(Q20_path, col_names = F, delim = "\t")
    project_tmp_Q20_rlt <- project_tmp_Q20 %>% mutate(Quality = "Q20", Confidence = "Probable")
    
    project_rlt_tmp <- rbind(project_tmp_Q0_rlt, project_tmp_Q20_rlt)
    
    
  } else {
    
    # read_Q0...
    Q0_path <- str_c(project_rlt_path, "/", project_files[1])
    project_tmp_Q0 <- read_delim(Q0_path, col_names = F, delim = "\t")
    project_tmp_Q0_rlt <- project_tmp_Q0 %>% mutate(Quality = "Q0", Confidence = "Possible")
    
    project_rlt_tmp <- rbind(project_tmp_Q0_rlt)
    
  }
  
  # data process
  
  det_sub <- c("ITS=", "LSU=")
  
  project_det <-
    project_rlt_tmp %>%
    group_by(Quality) %>%
    filter(!(X1 %in% need_row | X1 %in% det_sub)) %>%
    mutate(det_tmp = str_c(X1, X2, sep = " ")) %>%
    select(Quality, det_tmp) %>%
    mutate(det_tmp = str_c(det_tmp, collapse = ", ")) %>%
    distinct()
  
  project_det_sub <-
    project_rlt_tmp %>%
    group_by(Quality) %>%
    filter(X1 %in% det_sub) %>%
    mutate(det_tmp_sub = str_c(X1, X2, sep = " ")) %>%
    select(Quality, det_tmp_sub) %>%
    mutate(det_tmp_sub = str_c(det_tmp_sub, collapse = ", ")) %>%
    distinct()
  
  project_FRRN_rlt <-
    project_rlt_tmp %>% filter(X1 %in% need_row) %>%
    mutate(project = project_name, .before = X1) %>%
    pivot_wider(names_from = X1, values_from = X2) %>%
    left_join(project_det, by = "Quality") %>%
    left_join(project_det_sub, by = "Quality")
  
  colnames(project_FRRN_rlt) <- c("project", "Quality", "Accurary", "depth_avg",
                                  "ITS_only", "Both_ITS_LSU", "diff", "det_single", "det_multi")
  
  project_FRRN_rlt <-
    project_FRRN_rlt %>% select(project, ITS_only, Both_ITS_LSU, depth_avg, diff,
                                everything())
  
  # results
  return(project_FRRN_rlt)
  
}

FRRN_data_path_F %>% length()
FRRN_proj_name %>% length()

# future so fast !!!
FRRN_rlt_summary <- future_map2_dfr(FRRN_data_path_F, FRRN_proj_name, FRRN_rlt_process, .progress = T)
FRRN_rlt_summary


# use Q20 or only Q0 -----------
project_Q20_tab <- FRRN_rlt_summary %>% filter(Quality == "Q20")
project_Q0_tab <- FRRN_rlt_summary %>% filter(Quality == "Q0")

Q20_lst <- project_Q20_tab %>% select(project) %>% distinct() %>% pull()
project_Q0_tab_f <- project_Q0_tab %>% filter(!project %in% Q20_lst)

FRRN_rlt_summary1 <- rbind(project_Q0_tab_f, project_Q20_tab)

#FRRN_rlt_summary1 %>% filter(project == "Helsp1")
#FRRN_rlt_summary1 %>% filter(project == "LecAK0013_1")

# write_xlsx(FRRN_rlt_summary1, "./3.tables/FRRN_rlt_20260629.xlsx")
write_xlsx(FRRN_rlt_summary1, "./3.tables/FRRN_rlt_20260707.xlsx")
# results should be carefully checked #


# done...



