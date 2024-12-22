#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
setwd(here::here())
base_model_name <- "Omicron20"
model_name <- "M4"
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")
dir.create(dir_model_data, showWarnings = FALSE, recursive = TRUE)
check_redo <- TRUE
if(check_redo){
  file.remove(list.files(dir_model_data, pattern = "df_all_values_M4_\\d{1,2}.rds$", full.names = TRUE))
}

# Use the same fitting data as the M0 model
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

# Construct the parameter set with different combination of change_IDR, change_DSR, and traveler_weight.
# df_all_values <- tibble(
#   change_IDR = list(NA),
#   change_DSR = list(NA),
#   traveler_weight = list(NA),
#   id_unit_intro = list(which(data_fitting$country_under_investigation=="ZA")-1)
# )
# df_all_values$scenario <- "base"

df_all_values <- tibble()

values_traveler_weight <- c(1/1000, 1/100, seq(0.1, 0.9, 0.2), 1)
values_strategies <- c(expand.grid(c("T2", "T2a", "T3", "T3a", "T7", "T7a", "T11"), c("R", "H", "M")) %>% apply(1, paste, collapse="-"))
values_id_unit_intro <- seq(0, 28, 1)

source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
df_all_values_STp_TW_ID <- bind_rows(
  lapply(values_traveler_weight, function(traveler_weight_i){
    lapply(values_strategies, function(strategy_i){
      df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = "per_capita", travel_data_period = "normal")

      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR),
          change_DSR = list(df_tmp$change_DSR),
          traveler_weight = list(df_tmp$traveler_weight),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("STp_TW_ID_", strategy_i, "_", traveler_weight_i, "_", id_unit_intro)
        ))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) 
)
df_all_values <- bind_rows(df_all_values, df_all_values_STp_TW_ID)

df_all_values_STt_TW_ID <- bind_rows(
  lapply(values_traveler_weight, function(traveler_weight_i){
    lapply(values_strategies, function(strategy_i){
      df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = "total", travel_data_period = "normal")

      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR),
          change_DSR = list(df_tmp$change_DSR),
          traveler_weight = list(df_tmp$traveler_weight),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("STt_TW_ID_", strategy_i, "_", traveler_weight_i, "_", id_unit_intro)
        ))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) 
)
df_all_values <- bind_rows(df_all_values, df_all_values_STt_TW_ID)

## a reference scenario
df_all_values_STt_TW_ID_add <- lapply(values_id_unit_intro, function(id_unit_intro){
  traveler_weight_i <- 1/100
  return(tibble(
    change_IDR = list(1),
    change_DSR = list(1),
    traveler_weight = list(traveler_weight_i),
    id_unit_intro = list(id_unit_intro),
    scenario = paste0("STt_TW_ID_", "Max", "_", traveler_weight_i, "_", id_unit_intro)
  ))
}) %>% bind_rows()

df_all_values <- bind_rows(df_all_values, df_all_values_STt_TW_ID_add)

## adding a testing senario
# adding the second to seventh, rather than the first (T3 strategy), and test different combinations
values_strategies_add_one <- c(expand.grid(paste0("T2+", c(0, 2:21)), "H") %>% apply(1, paste, collapse="-")) 
df_all_values_STp_TW_ID_add_one <- bind_rows(
  lapply(values_traveler_weight, function(traveler_weight_i){
    # traveler_weight_i = values_traveler_weight[3]
    lapply(values_strategies_add_one, function(strategy_i){
      # strategy_i = values_strategies_add_one[11]
      df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = "per_capita", travel_data_period = "normal")

      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR),
          change_DSR = list(df_tmp$change_DSR),
          traveler_weight = list(df_tmp$traveler_weight),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("STp_TW_ID_", strategy_i, "_", traveler_weight_i, "_", id_unit_intro)
        ))
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) 
)
df_all_values <- bind_rows(df_all_values, df_all_values_STp_TW_ID_add_one) 

## save the values to N files
n_files <- 50

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
  paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1 ", n_files, " amd", " 32 ", model_name)
)
writeLines(bash_command, paste0("scripts/model_simulation/Simulate_", model_name, "_2_HPC_command.sh"))
system(paste0("chmod +x scripts/model_simulation/Simulate_", model_name, "_2_HPC_command.sh"))
