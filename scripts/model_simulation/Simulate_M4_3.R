#!/usr/bin/Rscript

pacman::p_load(tidyverse, ggplot2, patchwork, parallel, doParallel)

## Run on local machine
model_name <- "M4"
base_model_name <- "Omicron20"
dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")

n_cores <- 64

## Load data
values_sequencing_propensity <- c(0, 1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 0.9, 0.99, 1)
values_top_n <- c(3, 5, 7, 9, 11, 13, 15, 17, NA)
values_id_unit_intro <- seq(0, 28, 1)
df_all_values <- list.files(dir_model_data, "df_all_values_M4_\\d{1,2}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

df_all_values$change_DSR[df_all_values$scenario == "TH_SP_ID_9_0.01_0"]
df_all_values$sequencing_propensity[df_all_values$scenario == "TH_SP_ID_9_0.01_0"]

data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

## Below is for visualizing the results
source("scripts/model_simulation/helper/plot_fig_4abcd.R")

list_p_fig_4abcd <- plot_fig_4abcd(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values = df_all_values,
  data_fitting = data_fitting,
  bootstrap_each_param_pair=1024*10,
  n_cores = n_cores
)

## simulate the data for Fig 4c and 4d
list_p_fig_4abcd <- readRDS(paste0(dir_rst, "list_p_fig_4abcd.rds"))
source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
source("scripts/model_simulation/helper/prepare_simulation_fig_4ef.R")

prepare_simulation_fig_4ef(
  best_param_set = list_p_fig_4abcd$best_param_set %>% slice_sample(n=1),
  base_model_name = base_model_name,
  values_change_IDR = c(1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 1),
  values_change_DSR = c(1/1000, 1/100, 1/10, seq(0.2, 0.8, 0.2), 1),
  values_id_unit_intro = seq(0, 28, 1),
  data_fitting = data_fitting
)
