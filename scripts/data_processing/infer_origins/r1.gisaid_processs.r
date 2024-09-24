library(tidyverse)
library(parallel)
library(here) 
n_cores = 16
setwd("global_genomic_surveillance_spatpomp/")
# GISAID data
df_meta0 <- read_tsv("data/gisaid_data/metadata.tsv")
df_meta<-df_meta0
# select columns to use
df_meta <- df_meta %>% select(`Virus name`,`Accession ID`, `Collection date`, Host,Location, `Clade`, `Pango lineage`)

# extract continent, country, region
df_meta$continent <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][1], mc.cores = n_cores))
df_meta$country <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][2], mc.cores = n_cores))
df_meta$region <- unlist(mclapply(df_meta$Location, function(x) strsplit(x, " / ")[[1]][3], mc.cores = n_cores))
df_meta <- df_meta %>% select(-Location)

# filter collection dates
df_meta <- df_meta %>% mutate(date_decimal=decimal_date(lubridate::ymd(`Collection date`))) 
df_meta <- df_meta %>% arrange(`Collection date`, `Accession ID`)
df_meta <- df_meta %>% filter(`Collection date`>="2019-12-01")
df_meta <- df_meta %>% filter(`Collection date`<="2022-03-15")
df_meta <- df_meta %>% filter(`Host`=="Human")

# remove rows with NA in `Collection date` or `Pango lineage`
check_date <- is.na(df_meta$date_decimal)
sum(check_date)/nrow(df_meta)*100
check_lineage <- is.na(df_meta$`Pango lineage`) | (df_meta$`Pango lineage`=="Unassigned")
sum(check_lineage)/nrow(df_meta)*100
sum(check_lineage|check_date)/nrow(df_meta)*100
df_meta <- df_meta[!(check_lineage|check_date),]

names(df_meta)[names(df_meta)=="Pango lineage"] <- "lineage"
# fix country names
cross_check_table <- readxl::read_excel("data/our_airports_data/cross_check_table_ihme_input_completed.xlsx") %>% filter(level==3)
cross_check_table$code<-cross_check_table$country_iso3c
df_countries_gisaid <- df_meta %>% select(continent, country) %>% unique() 
df_countries_gisaid %>% filter(country %in% country[duplicated(country)]) # some countries have multiple continents, e.g. Turkey (Europe and Asia) and Haiti (North America and Caribbean)

df_countries_gisaid$code <- countrycode::countrycode(df_countries_gisaid$country, origin = 'country.name', destination = 'iso3c') # Belgium /Luxembourg, Bonaire, Canary Islands, Crimea, Kosovo, Micronesia, Saint Martin, Sint Eustatius not matched
df_countries_gisaid <- df_countries_gisaid %>% filter(code %in% cross_check_table$code) # keep those in studying coutries
cross_check_table$loc_name[!cross_check_table$code %in% df_countries_gisaid$code] # "San Marino" "Yemen" not in GISAID data

df_meta <- left_join(df_meta, df_countries_gisaid, c("continent", "country"))
df_meta <- df_meta %>% filter(!is.na(code))
write_csv(df_meta, here::here("data/gisaid_data/df_meta_clean.csv"))

# output
### different VOCs
BA1<-df_meta[grepl("^BA.1", df_meta$lineage),"Virus name"]
BA2<-df_meta[grepl("^BA.2", df_meta$lineage),"Virus name"] 
Delta1<-df_meta[grepl("^B.1.617", df_meta$lineage),"Virus name"] 
Delta2<-df_meta[grepl("^AY", df_meta$lineage),"Virus name"] 
Delta<-rbind(Delta1,Delta2)
others<-df_meta[!grepl("^BA.1|^BA.2|^B.1.617|^AY", df_meta$lineage),"Virus name"]
write_csv(BA1, "results/model_data/origin_infer/BA1.csv")
write_csv(BA2, "results/model_data/origin_infer/BA2.csv")
write_csv(Delta, "results/model_data/origin_infer/Delta.csv")
write_csv(others, "results/model_data/origin_infer/others.csv")

