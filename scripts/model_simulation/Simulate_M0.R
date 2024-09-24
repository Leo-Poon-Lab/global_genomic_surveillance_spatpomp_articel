#!/usr/bin/Rscript

## The script is used for simulation of the M0 model. We plot Figs of arrival time, detection lags, importation and exportation dynamics in M0.

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(spatPomp))
suppressPackageStartupMessages(library(doParallel))
ncores = 64
simulation_replicates = 100
setwd(here::here())
model_name <- "Omicron20"
sim_model_name <- paste0(model_name, "_simulation")
variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others")
num_of_strains = length(variants)
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("input dir must be supplied.n", call.=FALSE)
} else {
  input_dir <- args[1]
}
# input_dir = "results/model_data/model_simulation/Omicron20"
if(input_dir == "results/model_data/model_simulation/Omicron20"){
  simulation_replicates <- 1024
  print("simulation_replicates: 1024")
}
# select the highest loglik file in the input_dir
files_loglik <- list.files(input_dir, pattern = "^logLik_.+\\.txt$", full.names = TRUE)
logliks <- as.numeric(unlist(lapply(files_loglik, readLines)))
input_mle <- gsub("logLik_", "", files_loglik[which.max(logliks)])
# input_mle <- gsub("logLik_", "", files_loglik[2])
input_mle <- gsub(".txt$", ".rds", input_mle)

stopifnot(file.exists(input_mle))
print(paste0("input_mle: ", input_mle))
params_mle <- readRDS(input_mle)
if(grepl("profiling", input_mle)){
  dir_rst <- paste0(gsub("model_data/", "figs/", dirname(input_mle), fixed=T), "/simulation/")
  dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
  writeLines(input_mle, paste0(dir_rst, "input_mle_Fig1.txt"))
} else {
  dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
}

params_mle_original <- params_mle
if(!any(grepl("beta_naught",names(params_mle)))){
  params_mle <- unlist(params_mle$params[which.min(params_mle$logLiks$logLik), ])
} else {
  params_mle <- unlist(params_mle)
}
names(params_mle) <- names(params_mle_original)

# # rebuild the model to make sure it is the most recent version
# suppressWarnings(suppressMessages(source(paste0("scripts/model_fitting/build_", model_name, ".R"), echo = FALSE)))
# suppressWarnings(suppressMessages(source(paste0("scripts/model_fitting/build_", model_name, ".R"), echo = FALSE))) # rebuild the model to make sure it is the most recent version
model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))

# plug in MLE here
source("scripts/model_simulation/helper/put_params_into_model.R")
model_spatpomp <- put_params_into_model(model_spatpomp, params_mle)
source("scripts/model_fitting/helper/write_parameters_table.R")
write_parameter_tables(
  model_name,
  output_path = dir_rst,
  mod_params_fitted=params_mle,
  country_under_investigation=model_spatpomp@unit_names,
  n_of_units=length(unit_names(model_spatpomp)),
  cross_check_table=cross_check_table
)

# Get quantiles of simulations from the calibrated model.
set.seed(2023)
source("scripts/model_fitting/helper/sim_cases.R")
system.time(model_sims <- simulate_step_mc(model=model_spatpomp, model_nsim=simulation_replicates, n_cores=ncores, if_aggregate_sequenced_cases=FALSE))

model_quants <- get_quants_from_sim(
  model_sims, 
  events=c("daily_cases", "C_1_all", "C_2_all", "C_3_all", "C_4_all", "daily_deaths"),
  quants=c(0.05, 0.5, 0.95)
) %>% mutate(
  date = lubridate::date_decimal(date_decimal)
) %>% mutate(
  date = as.Date(lubridate::round_date(date, unit = 'day'))
)

suppressMessages(suppressWarnings(source("scripts/model_simulation/helper/plot_fig_1.R")))
if(grepl("profiling_[asu]", input_mle)){
  suppressMessages(suppressWarnings(plot_fig_1c(model_name, model_sims, model_quants, dir_rst)))
} else {
  suppressMessages(suppressWarnings(plot_fig_1(model_name, model_sims, model_quants, dir_rst)))
  source("scripts/model_simulation/helper/plot_fig_2.R")
  suppressMessages(suppressWarnings(plot_fig_2(model_sims, dir_rst)))
  saveRDS(model_sims %>% filter(`.id` %in% 1:100), "results/figs/model_simulation/Omicron20/model_sims.rds")
}
