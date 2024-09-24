library(tidyverse)
library(archive)
library(parallel)

n_cores = 16

# Inferred origin
files_origins_xz <- list.files(here::here("data/inferred_origins/"), pattern = "_origin.txt.xz", full.names = TRUE)
if(length(files_origins_xz)<1){ # have no origin data
  # source(here::here("scripts/data_processing/Prepare_origin_data.r"))
  warning("No origin data available, please run the shell script `scripts/data_processing/Infer_origins` first.")
  stop()
} else {
  df_origin <- mclapply(files_origins_xz, function(file_i){
    file_variant <- NA
    if(grepl("BA1_", file_i)){
      file_variant <- "BA.1"
    } else if(grepl("BA2_", file_i)){
      file_variant <- "BA.2"
    } else if(grepl("Delta_", file_i)){
      file_variant <- "Delta"
    } else if(grepl("Others_", file_i)){
      file_variant <- "Others"
    } else {
      stop("Unknown variant in ", file_i)
    }
    df <- read_delim(file_i, delim = " ")
    df$variant <- file_variant
    df
  }, mc.cores = n_cores) %>% bind_rows()
}

names(df_origin)[names(df_origin)=="virus"] <- "Accession ID"
names(df_origin)[names(df_origin)=="country_infer"] <- "parent_country"
names(df_origin)[names(df_origin)=="code"] <- "sequencing_country"
df_origin <- df_origin %>% select("Accession ID", "date", "parent_country", "sequencing_country", "variant")
df_origin %>% filter(date>=ymd("2021-09-15")) # 5,296,695 rows
write_csv(df_origin %>% select(`Accession ID`), here::here("data/gisaid_data/df_origin_accession_id_all.csv"))

# GISAID data
tar_xz_file <- here::here("data/gisaid_data/metadata_tsv_2023_10_16.tar.xz")
df_meta <- read_tsv(archive_read(tar_xz_file, file = "metadata.tsv"))

# select columns to use
df_meta$label <- gsub(" ","", df_meta$`Virus name`)
df_meta <- df_meta %>% select(`Accession ID`, `Collection date`, Location, `Clade`, `Pango lineage`, label)

# extract continent, country, region
df_meta$continent <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][1], mc.cores = n_cores))
df_meta$country <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][2], mc.cores = n_cores))
df_meta$region <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][3], mc.cores = n_cores))
df_meta <- df_meta %>% select(-Location)

# filter collection dates
df_meta <- df_meta %>% mutate(date_decimal=decimal_date(lubridate::ymd(`Collection date`))) 
df_meta <- df_meta %>% arrange(`Collection date`, `Accession ID`)
df_meta <- df_meta %>% filter(`Collection date`>="2019-12-01")

# remove rows with NA in `Collection date` or `Pango lineage`
check_date <- is.na(df_meta$date_decimal)
check_lineage <- is.na(df_meta$`Pango lineage`) | (df_meta$`Pango lineage`=="Unassigned")

df_meta <- df_meta[!(check_lineage|check_date),]
names(df_meta)[names(df_meta)=="Pango lineage"] <- "lineage"

# merge with origin data
df_meta <- left_join(df_meta, df_origin, by="Accession ID", multiple="first")

# fix country names
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)
df_countries_gisaid <- df_meta %>% select(continent, country) %>% unique() 
df_countries_gisaid %>% filter(country %in% country[duplicated(country)]) # some countries have multiple continents, e.g. Turkey (Europe and Asia) and Haiti (North America and Caribbean)

df_countries_gisaid$code <- countrycode::countrycode(df_countries_gisaid$country, origin = 'country.name', destination = 'iso2c') # Belgium /Luxembourg, Bonaire, Canary Islands, Crimea, Kosovo, Micronesia, Saint Martin, Sint Eustatius not matched
df_countries_gisaid <- df_countries_gisaid %>% filter(code %in% cross_check_table$code) # keep those in studying coutries
cross_check_table$loc_name[!cross_check_table$code %in% df_countries_gisaid$code] # "San Marino" "Yemen" not in GISAID data
df_countries_gisaid$code3 <- countrycode::countrycode(df_countries_gisaid$code, origin = 'iso2c', destination = 'iso3c')

df_meta <- left_join(df_meta, df_countries_gisaid, c("continent", "country"))
df_meta <- df_meta %>% filter(!is.na(code))

df_meta <- df_meta %>% filter(!is.na(parent_country))
df_meta$parent_country <- factor(df_meta$parent_country, levels=df_countries_gisaid$code3, labels=df_countries_gisaid$code)
df_meta$sequencing_country <- factor(df_meta$sequencing_country, levels=df_countries_gisaid$code3, labels=df_countries_gisaid$code)

df_meta$parent_country <- as.character(df_meta$parent_country)
df_meta$sequencing_country <- as.character(df_meta$sequencing_country)

# group by country and date variants, sum the counts
df_meta_lin_sum <- df_meta %>% group_by(`Collection date`, date_decimal, continent, sequencing_country, parent_country, variant) %>% summarise(N=n()) %>% ungroup()

# fix continent
unique(df_meta_lin_sum$continent)
df_meta_lin_sum$continent[df_meta_lin_sum$continent=="Caribbean"] <- "North America"
df_meta_lin_sum$continent[df_meta_lin_sum$sequencing_country=="TR"] <- "Asia"

# output
# df_meta_lin_sum$variant[df_meta_lin_sum$variant=="BA_1"] <- "BA.1"
# df_meta_lin_sum$variant[df_meta_lin_sum$variant=="BA_2"] <- "BA.2"
write_csv(df_meta_lin_sum, here::here("data/gisaid_data/df_meta_lin_sum.csv"))
