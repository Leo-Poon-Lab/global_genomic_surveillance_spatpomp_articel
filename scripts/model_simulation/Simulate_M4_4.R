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

df_all_values <- list.files(dir_model_data, "df_all_values_M4_\\d{4}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

## testing
print(paste0("DSR 0.4: ", df_all_values$change_DSR[df_all_values$scenario == "TH_SP_ID_DSR_9_1_28_0.4"]))
print(paste0("DSR 0.6: ", df_all_values$change_DSR[df_all_values$scenario == "TH_SP_ID_DSR_9_1_28_0.6"]))
print(paste0("DSR 0.8: ", df_all_values$change_DSR[df_all_values$scenario == "TH_SP_ID_DSR_9_1_28_0.8"]))
print(paste0("DSR 1: ", df_all_values$change_DSR[df_all_values$scenario == "TH_SP_ID_DSR_9_1_28_1"]))

print(paste0("IDR 0.4: ", df_all_values$change_IDR[df_all_values$scenario == "TH_SP_ID_IDR_9_1_28_0.4"]))
print(paste0("IDR 0.6: ", df_all_values$change_IDR[df_all_values$scenario == "TH_SP_ID_IDR_9_1_28_0.6"]))
print(paste0("IDR 0.8: ", df_all_values$change_IDR[df_all_values$scenario == "TH_SP_ID_IDR_9_1_28_0.8"]))
print(paste0("IDR 1: ", df_all_values$change_IDR[df_all_values$scenario == "TH_SP_ID_IDR_9_1_28_1"]))

## load Fig 4abcd
list_p_fig_4abcd <- readRDS(paste0(dir_rst, "list_p_fig_4abcd.rds"))

source("scripts/model_simulation/helper/plot_fig_4efg.R")
## Plotting Fig. 4e and 4f
list_fig_4efg <- plot_fig_4efg(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values = df_all_values,
  data_fitting = data_fitting,
  bootstrap_each_param=1024*50,
  n_cores = n_cores
)

## combine all fig 4 sub plots
p_fig_4 <- list_p_fig_4abcd$p_fig_4abcd / list_fig_4efg$p_fig_4efg
ggsave(paste0(dir_rst, "fig_4.jpg"), p_fig_4, width = 16, height = 20, dpi = 300)
ggsave(paste0(dir_rst, "fig_4.pdf"), p_fig_4, width = 16, height = 20)
