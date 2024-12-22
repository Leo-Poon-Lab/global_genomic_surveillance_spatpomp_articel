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

# Construct the parameter set with different combination of change_IDR, change_DSR, and traveler_weight.
values_change_IDR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_change_DSR <- c(1/1000, 1/100, 1/10, 0.5, 1, 1.5, 2)
values_traveler_weight <- c(1/1000, 1/100, seq(0.1, 0.9, 0.2), 1)
values_strategies <- c(expand.grid(c("T3", "T3a", "T7", "T7a", "T11", "T11a"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))

df_all_values <- bind_rows(
  expand.grid(
    change_IDR = values_change_IDR,
    change_DSR = NA,
    traveler_weight = NA
  ),
  expand.grid(
    change_IDR = NA,
    change_DSR = values_change_DSR,
    traveler_weight = NA
  ),
  expand.grid(
    change_IDR = NA,
    change_DSR = NA,
    traveler_weight = values_traveler_weight
  )
)
df_all_values <- as_tibble(unique(df_all_values))
df_all_values$scenario <- NA

df_all_values$scenario[is.na(df_all_values$change_DSR) & is.na(df_all_values$traveler_weight)] <- paste0("detection", "_", df_all_values$change_IDR[is.na(df_all_values$change_DSR) & is.na(df_all_values$traveler_weight)])
df_all_values$scenario[is.na(df_all_values$change_IDR) & is.na(df_all_values$traveler_weight)] <- paste0("sequencing", "_", df_all_values$change_DSR[is.na(df_all_values$change_IDR) & is.na(df_all_values$traveler_weight)])
df_all_values$scenario[is.na(df_all_values$change_IDR) & is.na(df_all_values$change_DSR)] <- paste0("TW", "_", df_all_values$traveler_weight[is.na(df_all_values$change_IDR) & is.na(df_all_values$change_DSR)])
df_all_values <- df_all_values %>% group_by(row_number()) %>% mutate(change_IDR=list(change_IDR), change_DSR=list(change_DSR), traveler_weight=list(traveler_weight)) %>% ungroup()

df_all_values$fig <- "fig_3a"

## Other potentially interesting scenarios (interacting between two parameters): (0) interaction between IDR and DSR; (1) traveler weight with IDR; (2) traveler weight with DSR; (3) traveler weight and different strategies.

## interaction between IDR and DSR
df_all_values_IDR_DSR <- bind_rows(
  lapply(values_change_IDR, function(change_IDR_i){
    lapply(values_change_DSR, function(change_DSR_i){
      return(tibble(
        change_IDR = list(change_IDR_i),
        change_DSR = list(change_DSR_i),
        traveler_weight = list(NA),
        scenario = paste0("detection_sequencing_", change_IDR_i, "_", change_DSR_i)
      ))
    })
  })
)
df_all_values_IDR_DSR$fig <- "fig_3b"
df_all_values <- bind_rows(df_all_values, df_all_values_IDR_DSR)

## traveler weight with IDR
df_all_values_TW_detection <- bind_rows(
  lapply(values_change_IDR, function(change_IDR_i){
    lapply(values_traveler_weight, function(traveler_weight_i){
      return(tibble(
        change_IDR = list(change_IDR_i),
        change_DSR = list(NA),
        traveler_weight = list(traveler_weight_i),
        scenario = paste0("TW_detection_", traveler_weight_i, "_", change_IDR_i)
      ))
    })
  })
)
df_all_values_TW_detection$fig <- "fig_3c"
df_all_values <- bind_rows(df_all_values, df_all_values_TW_detection)

## traveler weight with DSR
df_all_values_TW_sequencing <- bind_rows(
  lapply(values_change_DSR, function(change_DSR_i){
    lapply(values_traveler_weight, function(traveler_weight_i){
      return(tibble(
        change_IDR = list(NA),
        change_DSR = list(change_DSR_i),
        traveler_weight = list(traveler_weight_i),
        scenario = paste0("TW_sequencing_", traveler_weight_i, "_", change_DSR_i)
      ))
    })
  })
)
df_all_values_TW_sequencing$fig <- "fig_3d"
df_all_values <- bind_rows(df_all_values, df_all_values_TW_sequencing)

## traveler weight and strategies (per capita)
source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
df_all_values_STp_TW <- bind_rows(
  lapply(values_traveler_weight, function(traveler_weight_i){
    # traveler_weight_i=values_traveler_weight[5]
    lapply(values_strategies, function(strategy_i){
      # strategy_i = values_strategies[8]
      df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = "per_capita")
      return(tibble(
        change_IDR = list(df_tmp$change_IDR),
        change_DSR = list(df_tmp$change_DSR),
        traveler_weight = list(df_tmp$traveler_weight),
        scenario = paste0("STp_TW_", strategy_i, "_", traveler_weight_i),
        travel_hubs_selected = list(df_tmp$travel_hubs_selected[[1]])
      ))
    })
  })
)
df_all_values_STp_TW$fig <- "fig_3e"
df_all_values <- bind_rows(df_all_values, df_all_values_STp_TW)

## traveler weight and strategies (total)
df_all_values_STt_TW <- bind_rows(
  lapply(values_traveler_weight, function(traveler_weight_i){
    lapply(values_strategies, function(strategy_i){
      df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = "total")
      return(tibble(
        change_IDR = list(df_tmp$change_IDR),
        change_DSR = list(df_tmp$change_DSR),
        traveler_weight = list(df_tmp$traveler_weight),
        scenario = paste0("STt_TW_", strategy_i, "_", traveler_weight_i),
        travel_hubs_selected = list(df_tmp$travel_hubs_selected[[1]])
      ))
    })
  })
)
df_all_values_STt_TW$fig <- "fig_3f"
df_all_values <- bind_rows(df_all_values, df_all_values_STt_TW)

## save the values to 40 files
n_files <- 35

## cut the df_all_values into n_files files
split_indices <- cut(seq(1, nrow(df_all_values)), breaks = n_files, labels = FALSE)
for (i in 1:n_files) {
  df_all_values_i <- df_all_values[split_indices == i, ]
  saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", i, ".rds"))
}


## generate a bash command for submitting i slurm jobs
bash_command <- c(
  "#!/bin/bash",
  "",
  paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1 ", n_files, " amd", " 64 ", model_name)
)
writeLines(bash_command, "scripts/model_simulation/Simulate_M3_2_HPC_command.sh")
system("chmod +x scripts/model_simulation/Simulate_M3_2_HPC_command.sh")
system("chmod +x scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh")
