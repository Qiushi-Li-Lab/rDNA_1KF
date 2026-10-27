
##### adj AMP #####
# by Qiushi Li, IM-CAS, 2026.10.17


### opt ###
library(optparse)

# settings
options(rlang_backtrace_on_error = "none")

opt_list <- list(
  make_option(
    c("-r", "--ref_db"),
    type = "character",
    default = NULL,
    help = "RRNdf table, obtained from https://nmdc.cn/rrndf/"
  ),
  make_option(
    c("-i", "--table"),
    type = "character",
    default = NULL,
    help = "Un-rarefied OTU/ASV table, name of first column must be <OTU_ID>"
  ),
  make_option(
    c("-t", "--taxa_metadata"),
    type = "character",
    default = NULL,
    help = "must contain <OTU_ID> and may include some of <Phylum>, <Class>, <Order>, <Family>, <Genus>, <Species>"
  ),
  make_option(
    c("-m", "--adj_mode"),
    type = "character",
    default = "balance",
    help = "Adjustment mode: <strict> (un-targeted OTU/ASV set to -1), <balance> (only OTU/ASV matched to RRNdf are adjusted, others kept as original counts), <broad> (all OTU/ASV adjusted, un-targeted filled with the sample mean RRN)"
  ),
  make_option(
    c("-p", "--Phylum"),
    type = "character",
    default = "all",
    help = "Taxa to adjust: <all> or <Glomeromycota>"
  ),
  make_option(
    c("-o", "--out_put"),
    type = "character",
    default = "adj_table.tsv",
    help = "Adjusted OTU/ASV table, saved in TSV format."
  )
)


opt_parser <- 
  OptionParser(
    usage = "usage %prog [options]",
    option_list = opt_list,
    add_help_option = T,
    prog = NULL,
    description = "An R script for adjusting ASV/OTU table with RRNdf (obtained from https://nmdc.cn/rrndf/)."
  )

opt <- parse_args(opt_parser)

# check the opt
for (o in c("ref_db", "table", "taxa_metadata")) {
  if (is.null(opt[[o]])) {
    stop("Missing required argument: --", o, call. = FALSE)
  }
}

##### main #####
#### packages

library(tidyverse) # R4DS


#### input
input_table <- read_tsv(opt$table)
input_tax <- read_tsv(opt$taxa_metadata)
input_db <- read_tsv(opt$ref_db)
input_adj_mode <- opt$adj_mode
input_adj_Phylum <- opt$Phylum

#### check the input
taxa_cols <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# input_table：OTU_ID must be contained
if (!"OTU_ID" %in% names(input_table))
  stop("--table must contain a column named 'OTU_ID'", call. = FALSE)

# input_tax：OTU_ID must be contained
if (!"OTU_ID" %in% names(input_tax))
  stop("--taxa_metadata must contain a column named 'OTU_ID'", call. = FALSE)

taxa_ranks <- taxa_cols[taxa_cols %in% names(input_tax)]
if (length(taxa_ranks) == 0)
  stop("--taxa_metadata must contain at least one taxonomy column (",
       str_c(taxa_cols, collapse = "/"), ")", call. = FALSE)

# input_db：check the RRN col
if (!"RRN" %in% names(input_db))
  stop("--ref_db must contain a column named 'RRN'", call. = FALSE)

missing_db <- setdiff(taxa_ranks, names(input_db))
if (length(missing_db) > 0)
  stop("--ref_db is missing taxonomy column(s): ",
       str_c(missing_db, collapse = ", "), call. = FALSE)

# check the OTU table and taxa table
# if all OTUs is unique
stopifnot(!anyDuplicated(input_table$OTU_ID))
stopifnot(!anyDuplicated(input_tax$OTU_ID))

# is numeric ?
is_num <- vapply(input_table[-1], is.numeric, logical(1))
if (any(!is_num)) stop("Non-numeric sample columns: ", ...)
if (!is.numeric(input_db$RRN)) stop(...)


### RRNdf parse
input_db1 <- input_db %>% mutate(
  Phylum = str_sub(Phylum, start = 3),
  Class = str_sub(Class, start = 3),
  Order = str_sub(Order, start = 3),
  Family = str_sub(Family, start = 3),
  Genus = str_sub(Genus, start = 3)
)


# adj Phylum
adj_Phylum <- rlang::arg_match(input_adj_Phylum, c("all", "Glomeromycota"))

if(adj_Phylum != "all") {
  
  input_db1 <- input_db1 %>% filter(
    Phylum %in% adj_Phylum
  )
  
}


# check the Phylum
if (nrow(input_db1) == 0)
  stop("No records left in ref_db after filtering by Phylum = ",
       adj_Phylum, call. = FALSE)


# Species tab
RRNdf_spc <- input_db1 %>% mutate(
  Species = str_c(str_split_i(Species, pattern = " ", 1), str_split_i(Species, pattern = " ", 2), sep = " ")
) %>% group_by(Species) %>% summarise(RRN_m = mean(RRN))

# Genus tab
RRNdf_gen <- input_db1 %>% select(Genus, RRN) %>%
  group_by(Genus) %>% summarise(RRN_m = mean(RRN))

# Family tab
RRNdf_fam <- input_db1 %>% select(Family, RRN) %>%
  group_by(Family) %>% summarise(RRN_m = mean(RRN))

# Order tab
RRNdf_ord <- input_db1 %>% select(Order, RRN) %>%
  group_by(Order) %>% summarise(RRN_m = mean(RRN))

# Class tab
RRNdf_cla <- input_db1 %>% select(Class, RRN) %>%
  group_by(Class) %>% summarise(RRN_m = mean(RRN))

# Phylum tab
RRNdf_phy <- input_db1 %>% select(Phylum, RRN) %>%
  group_by(Phylum) %>% summarise(RRN_m = mean(RRN))
### done.


# Taxa RRN
# Species RRN
Table_RRN_spc <- input_tax %>% left_join(RRNdf_spc, by = "Species") %>% 
  filter(!is.na(RRN_m))

Table_RRN_spc_na <- input_tax %>% left_join(RRNdf_spc, by = "Species") %>% 
  filter(is.na(RRN_m))

# Genus RRN
Table_RRN_gen <- Table_RRN_spc_na %>% select(-RRN_m) %>% left_join(RRNdf_gen, by = "Genus") %>% 
  filter(!is.na(RRN_m))

Table_RRN_gen_na <- Table_RRN_spc_na %>% select(-RRN_m) %>% left_join(RRNdf_gen, by = "Genus") %>% 
  filter(is.na(RRN_m))

# Family RRN
Table_RRN_fam <- Table_RRN_gen_na %>% select(-RRN_m) %>% left_join(RRNdf_fam, by = "Family") %>% 
  filter(!is.na(RRN_m))

Table_RRN_fam_na <- Table_RRN_gen_na %>% select(-RRN_m) %>% left_join(RRNdf_fam, by = "Family") %>% 
  filter(is.na(RRN_m))

# Order RRN
Table_RRN_ord <- Table_RRN_fam_na %>% select(-RRN_m) %>% left_join(RRNdf_ord, by = "Order") %>% 
  filter(!is.na(RRN_m))

Table_RRN_ord_na <- Table_RRN_fam_na %>% select(-RRN_m) %>% left_join(RRNdf_ord, by = "Order") %>% 
  filter(is.na(RRN_m))

# Class RRN
Table_RRN_cla <- Table_RRN_ord_na %>% select(-RRN_m) %>% left_join(RRNdf_cla, by = "Class") %>% 
  filter(!is.na(RRN_m))

Table_RRN_cla_na <- Table_RRN_ord_na %>% select(-RRN_m) %>% left_join(RRNdf_cla, by = "Class") %>% 
  filter(is.na(RRN_m))

# Phylum RRN
Table_RRN_phy <- Table_RRN_cla_na %>% select(-RRN_m) %>% left_join(RRNdf_phy, by = "Phylum") %>% 
  filter(!is.na(RRN_m))

Table_RRN_phy_na <- Table_RRN_cla_na %>% select(-RRN_m) %>% left_join(RRNdf_phy, by = "Phylum") %>% 
  filter(is.na(RRN_m))



# parsed RRNdf table --------------------------
Table_RRN_full <- rbind(Table_RRN_phy, Table_RRN_cla, Table_RRN_ord,
                        Table_RRN_fam, Table_RRN_gen, Table_RRN_spc)

# write_tsv(Table_RRN_full, "adj_RRN_db.tsv")

#### table parse
### sample info
sample_id <- colnames(input_table)[-1]

input_table_RRN <- 
  input_table %>% 
  left_join(Table_RRN_full %>% select(OTU_ID, RRN_m), by = "OTU_ID")



RRN_adj <- function(OTU_TABLE_RRN, MODE = input_adj_mode, SCALE = 1000) {
  
  MODE <- rlang::arg_match(MODE, c("strict", "balance", "board"))
  
  # info
  na_ids <- OTU_TABLE_RRN$OTU_ID[is.na(OTU_TABLE_RRN$RRN_m)]
  n_otu <- nrow(OTU_TABLE_RRN)
  n_sample <- ncol(OTU_TABLE_RRN) - 2
  
  # fill
  filled_RRN <- OTU_TABLE_RRN %>% 
    summarise(
      across(
        -c(OTU_ID, RRN_m),
        # ~ mean(.x / .data[["RRN_m"]], na.rm = T),
        ~ sum(.x * RRN_m, na.rm = TRUE) / sum(.x[!is.na(RRN_m)], na.rm = TRUE)
      )
    )
  
  
  # OTU/ASV table with RRN_m
  adj_table <- OTU_TABLE_RRN %>% 
    mutate(
      across(
        .cols = -c(OTU_ID, RRN_m),
        .fns = ~ case_when(
          is.na(.data[["RRN_m"]]) & MODE == "strict" ~ -1,
          is.na(.data[["RRN_m"]]) & MODE == "balance" ~ round(.x * SCALE, 0),
          is.na(.data[["RRN_m"]]) & MODE == "board" ~ round(.x / filled_RRN[[cur_column()]] * SCALE, 0),
          TRUE ~ round(.x / .data[["RRN_m"]] * SCALE, 0)
          
        )
      )
    ) %>% 
    select(-RRN_m)
  
  # summary
  cli::cli_inform(c(
    "v" = "Adj with RRNdf finished",
    "*" = "Adj mode: {MODE}",
    "*" = "Adj Phylum: {adj_Phylum}",
    "*" = "SCALE: {format(SCALE, big.mark = ',', scientific = FALSE)}",
    "*" = "All OTU: {n_otu} OTU × {n_sample} sample",
    "*" = "Un-adjusted OTU: {length(na_ids)} OTU {cli::col_grey(if (length(na_ids) > 0) paste0(' (handled by ', MODE, ')') else '')}"
  ))
  
  return(adj_table)
  
}

# adj
adj_table <- RRN_adj(input_table_RRN)

# export
write_tsv(adj_table, opt$out_put)

# done...


