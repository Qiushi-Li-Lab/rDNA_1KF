
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


# rCNV fungi path & rCNV project name
rCNV_data_path <- dir("./JGI_rlt/", full.names = T)

# needed checked
str_detect(rCNV_data_path, "none", negate = F) %>% table()
str_detect(rCNV_data_path, "sobig", negate = F) %>% table()
str_detect(rCNV_data_path, "noITS", negate = F) %>% table()


rCNV_data_path_F <- rCNV_data_path[str_detect(rCNV_data_path, "none", negate = T)]
rCNV_data_path_F <- rCNV_data_path_F[str_detect(rCNV_data_path_F, "sobig", negate = T)]
rCNV_data_path_F <- rCNV_data_path_F[str_detect(rCNV_data_path_F, "noITS", negate = T)]

# rCNV_rlt_path

rCNV_proj_name <- str_split(rCNV_data_path_F, pattern = "\\.", n = 3, simplify = T)[, 3]

### view(rCNV_proj_name)
### view(rCNV_rlt_path)

view(rCNV_proj_name)


# clean & extract
# read the data
# tmp path

# list.files(rCNV_rlt_path[1])

rCNV_rlt_process <- function(project_path, project_name) {
  
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
  
  project_rCNV_rlt <-
    project_rlt_tmp %>% filter(X1 %in% need_row) %>%
    mutate(project = project_name, .before = X1) %>%
    pivot_wider(names_from = X1, values_from = X2) %>%
    left_join(project_det, by = "Quality") %>%
    left_join(project_det_sub, by = "Quality")
  
  colnames(project_rCNV_rlt) <- c("project", "Quality", "Accurary", "depth_avg",
                                  "ITS_only", "Both_ITS_LSU", "diff", "det_single", "det_multi")
  
  project_rCNV_rlt <-
    project_rCNV_rlt %>% select(project, ITS_only, Both_ITS_LSU, depth_avg, diff,
                                everything())
  
  # results
  return(project_rCNV_rlt)
  
}

rCNV_data_path_F %>% length()
rCNV_proj_name %>% length()

# future so fast !!!
rCNV_rlt_summary <- future_map2_dfr(rCNV_data_path_F, rCNV_proj_name, rCNV_rlt_process, .progress = T)
rCNV_rlt_summary

# view(rCNV_rlt_summary)
write_xlsx(rCNV_rlt_summary, "./output/rCNV_rlt_20250223.xlsx")

# done.

