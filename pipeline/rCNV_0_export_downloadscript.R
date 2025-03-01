
# export download JGI Raw-data scripts
# by Qiushi-Li, IM-CAS
# 2024.08.14

# packages
library(xml2)
library(tidyverse)

# read xml file
xml_file <- commandArgs(trailingOnly = T)
project_xml <- read_xml(xml_file)
# project_xml <- read_xml("Tercla1.xml")

# get proj nodes
proj_name <- xml_find_all(project_xml, xpath = "./@name") %>% xml_text()

# get rawdata nodes
raw_data_url <- xml_find_all(project_xml, xpath = "//folder[@name='Raw Data']") %>% 
  xml_find_all(xpath = ".//@url") %>% xml_text()
raw_data_md5 <- xml_find_all(project_xml, xpath = "//folder[@name='Raw Data']") %>% 
  xml_find_all(xpath = ".//@md5") %>% xml_text()


# http
raw_data_http <- str_c("https://genome-downloads.jgi.doe.gov" , raw_data_url)
raw_data_filename <- str_split(raw_data_http, pattern = "/", simplify = T)[, 12]

downloaded_http <- str_c("curl '",
                        raw_data_http,
                        "'",
                        " -b cookies > ",
                        raw_data_filename)

downloaded_http1 <- downloaded_http[str_detect(downloaded_http, "fastq")]
downloaded_http2 <- downloaded_http1[!str_detect(downloaded_http1, "filter")]

# sff
downloaded_http3 <- c(downloaded_http2, downloaded_http[str_detect(downloaded_http, "sff")])

proj_message <- str_c("###", proj_name, sep = " ")


### scripts
download_scripts <- c(
  "#! /usr/bin/bash",
  proj_message,
  "### login",
  "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=your_login' --data-urlencode 'password=your_passwd' -c cookies > /dev/null",
  "### download",
  downloaded_http3,
  "### md5 summary",
  "md5sum * > md5summary.txt",
  "Rscript ../../0.scripts/rCNV_check_md5.R",
  "echo 'done!' "
)

### file name
script_name <- str_c("JGI_download", "_", proj_name, ".sh")

### export
write_lines(download_scripts, script_name)

### md5
raw_data_md5 <- 
  data.frame(
    fil = raw_data_filename,
    md5 = raw_data_md5
  )
raw_data_md5_f <- raw_data_md5 %>% 
  filter(str_detect(fil, "fastq")) %>% 
  filter(!str_detect(fil, "filter"))

write.table(raw_data_md5_f, "raw_data_md5.txt", quote = F, row.names = F, sep = "\t")

# done.

