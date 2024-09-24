#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
setwd(here::here())
base_model_name <- "Omicron20"
model_name <- "M3"
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")
dir.create(dir_model_data, showWarnings = FALSE, recursive = TRUE)
check_redo <- TRUE
if(check_redo){
  file.remove(list.files(dir_model_data, pattern = "df_all_values_M3_\\d+.rds$", full.names = TRUE))
}

# Use the same fitting data as the M0 model
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

# Construct the parameter set with different combination of change_IDR, change_DSR, and sequencing_propensity.
values_change_IDR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_change_DSR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_sequencing_propensity <- c(0, 1/1000, 1/100, 1/10, 0.3, 0.5, 0.7, 0.9, 0.99, 1)
values_top_n <- c(3, 5, 7, 9, 11, 13, 15, 17, NA)

df_all_values <- bind_rows(
  expand.grid(
    change_IDR = values_change_IDR,
    change_DSR = NA,
    sequencing_propensity = NA
  ),
  expand.grid(
    change_IDR = NA,
    change_DSR = values_change_DSR,
    sequencing_propensity = NA
  ),
  expand.grid(
    change_IDR = NA,
    change_DSR = NA,
    sequencing_propensity = values_sequencing_propensity
  )
)
df_all_values <- as_tibble(unique(df_all_values))
df_all_values$scenario <- NA

df_all_values$scenario[is.na(df_all_values$change_DSR) & is.na(df_all_values$sequencing_propensity)] <- paste0("detection", "_", df_all_values$change_IDR[is.na(df_all_values$change_DSR) & is.na(df_all_values$sequencing_propensity)])
df_all_values$scenario[is.na(df_all_values$change_IDR) & is.na(df_all_values$sequencing_propensity)] <- paste0("sequencing", "_", df_all_values$change_DSR[is.na(df_all_values$change_IDR) & is.na(df_all_values$sequencing_propensity)])
df_all_values$scenario[is.na(df_all_values$change_IDR) & is.na(df_all_values$change_DSR)] <- paste0("SP", "_", df_all_values$sequencing_propensity[is.na(df_all_values$change_IDR) & is.na(df_all_values$change_DSR)])
df_all_values <- df_all_values %>% group_by(row_number()) %>% mutate(change_IDR=list(change_IDR), change_DSR=list(change_DSR), sequencing_propensity=list(sequencing_propensity)) %>% ungroup()

df_all_values$fig <- "fig_3a"

## prepare for the travel hubs scenario. For this scenario, we will try to keep the total detection and sequencing amount the same but increase the IDR and DSR among travel hubs (top 3, 5, 7, 9, 11, 13, 15, 17), and perform no sequencing at all for other units. Top travel hubs have the highest incoming traffic flow.
source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
df_all_values_travel_hubs <- bind_rows(
  lapply(values_top_n, function(top_n){
    df_tmp <- get_params_for_travel_hubs(top_n = top_n, base_model_name = base_model_name, data_fitting = data_fitting)
    return(tibble(
      change_IDR = list(df_tmp$change_IDR),
      change_DSR = list(df_tmp$change_DSR),
      sequencing_propensity = list(NA),
      scenario = paste0("TH_", top_n)
    ))
  })
)
df_all_values_travel_hubs$fig <- "fig_3a"
df_all_values <- bind_rows(df_all_values, df_all_values_travel_hubs)

## Other potentially interesting scenarios (interacting between two parameters): (0) interaction between IDR and DSR; (1) sequencing propensity with IDR; (2) sequencing propensity with DSR; (3) travel hubs with IDR; (4) travel hubs with DSR; (5) travel hubs with sequencing capacity. (use stat_density_2d(aes(fill = ..density..))

## interaction between IDR and DSR
df_all_values_IDR_DSR <- bind_rows(
  lapply(values_change_IDR, function(change_IDR_i){
    lapply(values_change_DSR, function(change_DSR_i){
      return(tibble(
        change_IDR = list(change_IDR_i),
        change_DSR = list(change_DSR_i),
        sequencing_propensity = list(NA),
        scenario = paste0("detection_sequencing_", change_IDR_i, "_", change_DSR_i)
      ))
    })
  })
)
df_all_values_IDR_DSR$fig <- "fig_3b"
df_all_values <- bind_rows(df_all_values, df_all_values_IDR_DSR)

## sequencing propensity with IDR
df_all_values_SP_detection <- bind_rows(
  lapply(values_change_IDR, function(change_IDR_i){
    lapply(values_sequencing_propensity, function(sequencing_propensity_i){
      return(tibble(
        change_IDR = list(change_IDR_i),
        change_DSR = list(NA),
        sequencing_propensity = list(sequencing_propensity_i),
        scenario = paste0("SP_detection_", sequencing_propensity_i, "_", change_IDR_i)
      ))
    })
  })
)
df_all_values_SP_detection$fig <- "fig_3c"
df_all_values <- bind_rows(df_all_values, df_all_values_SP_detection)

## sequencing propensity with DSR
df_all_values_SP_sequencing <- bind_rows(
  lapply(values_change_DSR, function(change_DSR_i){
    lapply(values_sequencing_propensity, function(sequencing_propensity_i){
      return(tibble(
        change_IDR = list(NA),
        change_DSR = list(change_DSR_i),
        sequencing_propensity = list(sequencing_propensity_i),
        scenario = paste0("SP_sequencing_", sequencing_propensity_i, "_", change_DSR_i)
      ))
    })
  })
)
df_all_values_SP_sequencing$fig <- "fig_3d"
df_all_values <- bind_rows(df_all_values, df_all_values_SP_sequencing)

## travel hubs with IDR
df_all_values_TH_IDR <- bind_rows(
  lapply(values_change_IDR, function(change_IDR_i){
    lapply(values_top_n, function(top_n){
      df_tmp <- get_params_for_travel_hubs(top_n = top_n, base_model_name = base_model_name, data_fitting = data_fitting)
      return(tibble(
        change_IDR = list(df_tmp$change_IDR*change_IDR_i),
        change_DSR = list(df_tmp$change_DSR),
        sequencing_propensity = list(NA),
        scenario = paste0("TH_detection_", top_n, "_", change_IDR_i)
      ))
    })
  })
)
df_all_values_TH_IDR$fig <- "fig_3e"
df_all_values <- bind_rows(df_all_values, df_all_values_TH_IDR)

## travel hubs with DSR
df_all_values_TH_DSR <- bind_rows(
  lapply(values_change_DSR, function(change_DSR_i){
    lapply(values_top_n, function(top_n){
      df_tmp <- get_params_for_travel_hubs(top_n = top_n, base_model_name = base_model_name, data_fitting = data_fitting)
      return(tibble(
        change_IDR = list(df_tmp$change_IDR),
        change_DSR = list(df_tmp$change_DSR*change_DSR_i),
        sequencing_propensity = list(NA),
        scenario = paste0("TH_sequencing_", top_n, "_", change_DSR_i)
      ))
    })
  })
)
df_all_values_TH_DSR$fig <- "fig_3f"
df_all_values <- bind_rows(df_all_values, df_all_values_TH_DSR)

## travel hubs with sequencing capacity
df_all_values_TH_SP <- bind_rows(
  lapply(values_sequencing_propensity, function(sequencing_propensity_i){
    lapply(values_top_n, function(top_n){
      df_tmp <- get_params_for_travel_hubs(top_n = top_n, base_model_name = base_model_name, data_fitting = data_fitting)
      return(tibble(
        change_IDR = list(df_tmp$change_IDR),
        change_DSR = list(df_tmp$change_DSR),
        sequencing_propensity = list(sequencing_propensity_i),
        scenario = paste0("TH_SP_", top_n, "_", sequencing_propensity_i)
      ))
    })
  })
)
df_all_values_TH_SP$fig <- "fig_3g"
df_all_values <- bind_rows(df_all_values, df_all_values_TH_SP)

## save the values to N files with 32 rows each for parallel computing
n_files <- ceiling(nrow(df_all_values)/32)
lapply(seq_len(n_files), function(i){
  df_all_values_i <- df_all_values[(i*32-31):min(i*32, nrow(df_all_values)), ]
  saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_M3_", i, ".rds"))
  return(NULL)
}) %>% invisible()

## generate a bash command for submitting i slurm jobs
bash_command <- c(
  "#!/bin/bash",
  "",
  paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1 ", n_files, " amd", " 64 ", model_name)
)
writeLines(bash_command, "scripts/model_simulation/Simulate_M3_2_HPC_command.sh")
system("chmod +x scripts/model_simulation/Simulate_M3_2_HPC_command.sh")
system("chmod +x scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh")
