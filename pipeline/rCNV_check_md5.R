

### check md5
# by Qiushi-Li
# IM-CAS, 2024.08.14

# packages
library(tidyverse)

# read the data
raw_data_md5 <- read.table("raw_data_md5.txt", header = T)
download_md5 <- read.table("md5summary.txt")

colnames(download_md5) <- c("md5", "fil")

# check md5
raw_data_md5_check <- raw_data_md5 %>% 
  left_join(download_md5, by = "fil") %>% 
  mutate(check_md5 = if_else(md5.x == md5.y, "Yes !", "Oh, No !"))

# result
write.table(raw_data_md5_check, "md5_check.txt", quote = F, row.names = F, sep = "\t")


# done.