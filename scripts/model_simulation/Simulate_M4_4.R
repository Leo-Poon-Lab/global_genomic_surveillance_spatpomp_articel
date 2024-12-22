#!/usr/bin/Rscript

pacman::p_load(tidyverse, ggplot2, patchwork, parallel, doParallel)

## Run on local machine
model_name <- "M4"
base_model_name <- "Omicron20"
n_cores <- 64
dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")

# ## Load data
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

df_all_values <- list.files(dir_model_data, "df_all_values_M4_1\\d{3}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

## load Fig 4abcd
list_p_fig_4abcd <- readRDS(paste0(dir_rst, "list_p_fig_4abcd.rds"))

source("scripts/model_simulation/helper/plot_fig_5abc.R")
## Plotting Fig. 5a and 5bc
list_fig_5abc <- plot_fig_5abc(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values = df_all_values,
  data_fitting = data_fitting,
  bootstrap_each_param=1024*10,
  n_cores = n_cores
)

## prepare the data for the next step - round 3 simulations
### round three testing changes in viral transmissibility and population immunity

list_fig_5abc <- readRDS(paste0(dir_rst, "list_p_fig_5abc.rds"))
best_param_set <- list_p_fig_4abcd$df_EDT_fig_4bcd_percentile %>%
  filter((parval_2=="P2-H" & parval_1=="0.1%") | (parval_2=="T2-R" & parval_1=="70%")) %>% 
  mutate(parval_1 = as.numeric(gsub("%","",as.character(parval_1)))/100) 

write_csv(best_param_set, paste0(dir_rst, "best_param_set.csv"))

source("scripts/model_simulation/helper/prepare_simulation_fig_5d.R")

df_all_values_new <- prepare_simulation_fig_5d(
  best_param_set = best_param_set,
  base_model_name = base_model_name,
  values_id_unit_intro = seq(0, 28, 1),
  values_change_Reff = c(1.5, 3, 5),
  values_change_Immunity_setting = c("Few", "Moderate", "More"),
  change_IDR = list_fig_5abc$change_IDR_min_better,
  change_DSR = list_fig_5abc$change_DSR_min_better,
  values_Scenarios = c("fewer_diagnostics", "fewer_sequencing", "half_resources"),
  data_fitting = data_fitting,
  check_redo = TRUE
)
