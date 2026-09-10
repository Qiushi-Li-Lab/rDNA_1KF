

# this R-script could aggregate busco results into a table-like sheets
# by Qiushi-Li, IMCAS
# 2026.02.11



# packages
library(tidyverse)

# read & writexl
library(readxl)
library(writexl)

# multi session calcu
library(furrr)
library(future)
plan(multisession)

# busco result fullpath
busco_rlt_fullpath <- list.files("./1.data/busco_summaries/", full.names = T)

# processing function
busco_rlt_process <- function(full_file_path) {
  
  # result file info
  proj_id <- str_split_i(full_file_path, pattern = "\\/", 4) %>% str_split_i(pattern = "_busco.txt", 1)
  
  # read logs
  busco_rlt <- read_log(full_file_path)
  
  # extract rlt element
  comp_tmp <- busco_rlt$X1[9]
  
  # split element
  comp_tmp1 <- str_split_i(comp_tmp, pattern = "\\[", 1) %>% 
    str_split_i(pattern = ":", 2) %>% parse_number()
  
  # trans into numeric
  comp_rlt <- comp_tmp1 / 100
  
  # results data.frame
  rlt_df <- data.frame(
    project_id = proj_id,
    completeness = comp_rlt
  )
  
  return(rlt_df)
  
}


busco_summary <- future_map_dfr(busco_rlt_fullpath, busco_rlt_process, .progress = T)

# view(busco_summary)

write_xlsx(busco_summary, "./3.tables/busco_summary_20260707.xlsx")

# checked
FRRN_rlt_taxa_spl1 <- FRRN_rlt_taxa_spl %>% left_join(busco_summary, by = c("project" = "project_id"))

# done.


