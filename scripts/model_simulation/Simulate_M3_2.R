#!/usr/bin/Rscript

## Run on HPC
setwd(here::here())
path_mle <- "results/model_data/profiling_HKU/Omicron20/profiling_shared_11/mlev1_high.rds"
params_mle <- readRDS(path_mle)

simulation_replicates <- 512
check_redo_model <- FALSE
check_redo_simulation <- FALSE
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("input df_all_values rds, and ncores must be supplied.n", call.=FALSE)
} else {
  input_rds <- args[1] # for df_all_values_M3_i.rds
  stopifnot(file.exists(input_rds))
  
  n_cores <- as.numeric(args[2])
  stopifnot(!is.na(n_cores))
  stopifnot(n_cores > 0)
}
# input_rds <- "results/model_data/model_simulation/M3/df_all_values_M3_49.rds"
# n_cores <- 64

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doParallel))

base_model_name <- "Omicron20"
model_name <- "M3"
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")
dir.create(dir_model_data, showWarnings = FALSE, recursive = TRUE)

data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

df_all_values <- readRDS(input_rds)

# load introduction dates
lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
df_vline <- tibble(lag=lags_intro)
df_vline$date_decimal <- floor(df_vline$lag)/365.25 + data_fitting$dates$date_start_d
df_vline$date <- date(date_decimal(df_vline$date_decimal))
df_vline$variant <- c("Omicron BA.1", "Omicron BA.2")

# simulate the models
source(here::here(paste0("scripts/model_building/spatPomp_", model_name, ".R")))
source("scripts/model_fitting/helper/sim_cases.R")
source("scripts/model_simulation/helper/put_params_into_model.R")
source("scripts/model_fitting/helper/fixed_parameters.R")

set.seed(2023)

## first generate and save the models in parallel
print("Generating the models in parallel")
print(date())
registerDoParallel(n_cores)

foreach(i = 1:nrow(df_all_values)) %dopar% {
  # i=1
  change_IDR_i <- unlist(df_all_values$change_IDR[i])
  change_DSR_i <- unlist(df_all_values$change_DSR[i])
  traveler_weight_i <- unlist(df_all_values$traveler_weight[i])
  model_name_i <- paste0("M3_", df_all_values$scenario[i])

  model_out_file_i <- paste0(dir_model_data, "Spatpomp_", model_name_i, ".rds")
  if(check_redo_model){
    check_redo_i <- TRUE
  } else{
    if(file.exists(model_out_file_i)){
      check_redo_i <- FALSE
      print(paste("Model", model_name_i, "already exists. Skip."))
    } else{
      check_redo_i <- TRUE
    }
  }

  if(check_redo_i){
    ## build model
    print(paste("Building the", i, "th model:", model_name_i))
    model_M3 <- suppressWarnings(suppressMessages(
      build_spatPomp_M3(
        global_data=data_fitting,
        dt_years = 1 / 365, 
        variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others"),
        fixed_shared_params_values = fixed_shared_params_values,
        fixed_unit_specific_params_values = fixed_unit_specific_params_values,
        N_introdution = 2000, # the number of Omicron BA.1 cases introduced into E and I compartments in South Africa, 1000 correspond to 3% of the infected population on day 0
        for_IBPF = FALSE,
        change_IDR = change_IDR_i,
        change_DSR = change_DSR_i,
        traveler_weight = traveler_weight_i
      )
    ))

    ## plug in MLE here
    model_M3 <- put_params_into_model(model_M3, params_mle)
    
    ## testing
    # system.time(print(paste("bpfilter logLik for", model_name_i, "model:", logLik(bpfilter(model_M3, Np=1, block_size=1)))))

    saveRDS(model_M3, model_out_file_i)
    rm(model_M3)
  }
} %>% invisible()
stopImplicitCluster()
print(date())
print("Finished generating the models in parallel")

## secondly simulate the models in parallel
print("Simulating the models in parallel")
print(date())
set.seed(NULL)
for(i in sample(1:nrow(df_all_values))){
  change_IDR_i <- unlist(df_all_values$change_IDR[i])
  change_DSR_i <- unlist(df_all_values$change_DSR[i])
  traveler_weight_i <- unlist(df_all_values$traveler_weight[i])
  model_name_i <- paste0("M3_", df_all_values$scenario[i])

  model_out_file_i <- paste0(dir_model_data, "Spatpomp_", model_name_i, ".rds")
  df_sims_out_file_i <- paste0(dir_model_data, "df_sims_", model_name_i, ".rds")

  if(check_redo_simulation){
    check_redo_i <- TRUE
  } else{
    if(file.exists(df_sims_out_file_i)){
      check_redo_i <- FALSE
      print(paste("df_sims for Model", model_name_i, "already exists. Skip."))
    } else{
      check_redo_i <- TRUE
    }
  }

  if(check_redo_i){
    ## simulation process
    model_M3 <- readRDS(model_out_file_i)

    set.seed(2023)
    print(paste("Simulating the", i, "th model:", model_name_i))
    print(date())
    df_sims <- simulate_step_mc(model=model_M3, model_nsim=simulation_replicates, n_cores=n_cores, if_aggregate_sequenced_cases=TRUE, if_return_EDT=TRUE, df_vline = df_vline)
    saveRDS(df_sims, df_sims_out_file_i)
    rm(df_sims, model_M3)
    print(date())
    print(paste("Finished simulating the", i, "th model:", model_name_i))
    gc()
  }
}

print(date())
print("Finished simulating the models in parallel")
