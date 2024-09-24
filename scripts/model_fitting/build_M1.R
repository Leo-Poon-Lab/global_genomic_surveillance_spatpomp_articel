#!/usr/bin/Rscript

#### Note: This script is used to fit the M_1 model ####
library(tidyverse)
library(slider)
library(doParallel)
source(here::here("scripts/data_processing/install_prerequisite.R"))
setwd(here::here())
check_redo <- TRUE
# check_redo <- FALSE
model_name <- "M1"
base_model_name <- "Omicron20"
variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others")
# colors_lineage <- MetBrewer::met.brewer("Degas", 4, "continuous")
colors_lineage <- RColorBrewer::brewer.pal(4, "Set1")
names(colors_lineage) <- variants
n_cores = 48 # use 48 or 72 when running on the "low" or "mid" search profile
# n_cores = 96 # use 96 when running on the "high" search profile

# modeling the Omicron variant emergence and circulation in the first 4 months #
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

# Use the same fitting data as the M0 model
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

source("scripts/model_fitting/helper/fixed_parameters.R")

## This step do Model building ##
source(here::here(paste0("scripts/model_building/spatPomp_", model_name, ".R")))
model_M1 <- build_spatPomp_M1(
  global_data=data_fitting,
  dt_years = 1 / 365, 
  variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others"),
  fixed_shared_params_values = fixed_shared_params_values,
  fixed_unit_specific_params_values = fixed_unit_specific_params_values,
  N_introdution = 2000, # the number of Omicron BA.1 cases introduced into E and I compartments in South Africa, 1000 correspond to 3% of the infected population on day 0
  for_IBPF = FALSE
)

# params for testing
if(any(coef(model_M1)==0)){
  coef(model_M1)[coef(model_M1)==0] <- 0.5
}

# # For testing
time_used <- system.time(print(paste("bpfilter logLik for", model_name, "model:", logLik(bpfilter(model_M1, Np=1, block_size=1)))))
print(paste0("Time used for logLik bpfilter:"))
print(time_used)

saveRDS(model_M1, paste0("results/model_data/model_", model_name, ".rds"))
