#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
setwd(here::here())
base_model_name <- "Omicron20"
model_name <- "M4"
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")
dir.create(dir_model_data, showWarnings = FALSE, recursive = TRUE)
check_redo <- TRUE
if(check_redo){
  file.remove(list.files(dir_model_data, pattern = "df_all_values_M4_\\d+.rds$", full.names = TRUE))
}

# Use the same fitting data as the M0 model
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

# Construct the parameter set with different combination of change_IDR, change_DSR, and sequencing_propensity.
df_all_values <- tibble(
  change_IDR = list(NA),
  change_DSR = list(NA),
  sequencing_propensity = list(NA),
  id_unit_intro = list(which(data_fitting$country_under_investigation=="ZA")-1)
)
df_all_values$scenario <- "base"

values_sequencing_propensity <- c(0, 1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 0.9, 0.99, 1)
source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
values_top_n <- c(3, 5, 7, 9, 11, 13, 15, 17, NA)
values_id_unit_intro <- seq(0, 28, 1)

df_all_values_TH_SP_ID <- bind_rows(
  lapply(values_sequencing_propensity, function(sequencing_propensity_i){
    lapply(values_top_n, function(top_n){
      df_tmp <- get_params_for_travel_hubs(top_n = top_n, base_model_name = base_model_name, data_fitting = data_fitting)
      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR),
          change_DSR = list(df_tmp$change_DSR),
          sequencing_propensity = list(sequencing_propensity_i),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("TH_SP_ID_", top_n, "_", sequencing_propensity_i, "_", id_unit_intro)
        ))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) 
)
df_all_values <- bind_rows(df_all_values, df_all_values_TH_SP_ID)

## save the values to N files with 32 rows each for parallel computing, with more rows if total files >50
if(nrow(df_all_values) <= 50*32){
  n_files <- ceiling(nrow(df_all_values)/32)
  lapply(seq_len(n_files), function(i){
    df_all_values_i <- df_all_values[(i*32-31):min(i*32, nrow(df_all_values)), ]
    saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", i, ".rds"))
    return(NULL)
  }) %>% invisible()
} else{
  n_files <- 50
  nrows_per_file <- ceiling(nrow(df_all_values)/n_files)
  ## cut the df_all_values into 50 files
  lapply(seq_len(n_files), function(i){
    df_all_values_i <- df_all_values[(i*nrows_per_file-nrows_per_file+1):min(i*nrows_per_file, nrow(df_all_values)), ]
    saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", i, ".rds"))
    return(NULL)
  }) %>% invisible()
  
}

## generate a bash command for submitting i slurm jobs
bash_command <- c(
  "#!/bin/bash",
  "",
  paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1 ", n_files, " amd", " 32 ", model_name)
)
writeLines(bash_command, paste0("scripts/model_simulation/Simulate_", model_name, "_2_HPC_command.sh"))
system(paste0("chmod +x scripts/model_simulation/Simulate_", model_name, "_2_HPC_command.sh"))
