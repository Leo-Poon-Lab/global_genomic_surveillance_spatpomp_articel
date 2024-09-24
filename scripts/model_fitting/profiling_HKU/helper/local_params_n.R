#!/usr/bin/Rscript

#### Note: This script is used to fit the Omicron model basing on previous profiling results ####

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=9) {
  stop("input mle, local_param_type, output dir, search_level, partition, ncores, task_id, and optimizer must be supplied.", call.=FALSE)
} else {
  input_mle <- args[1]
  stopifnot(file.exists(input_mle))

  local_param_type <- args[2]
  output_dir <- args[3]
  search_level <- args[4]
  partition <- args[5]

  n_cores <- as.numeric(args[6])
  stopifnot(!is.na(n_cores))
  stopifnot(n_cores > 0)

  task_id <- as.numeric(args[7])

  optimizer <- args[8]

  model_name <- args[9]

}

if(partition == "amd"){
  random_seed <- task_id + n_cores
} else if(partition == "intel"){
  random_seed <- task_id + 1000 + n_cores
} else if(partition == "condo_amd"){
  random_seed <- task_id + 2000 + n_cores
} else if(partition == "self"){
  random_seed <- task_id + 3000 + n_cores
} else if(partition == "hugemem"){
  random_seed <- task_id + 4000 + n_cores
} else {
  stop("partition must be one of 'amd', 'intel', 'self', 'hugemem', or 'condo_amd'.")
}
print(paste0("input_mle: ", input_mle, ", local_param_type: ", local_param_type, ", output_dir: ", output_dir, ", search_level: ", search_level, ", partition: ", partition, ", n_cores: ", n_cores, ", task_id: ", task_id, ", optimizer: ", optimizer))
# input_mle <- "results/model_data/profiling_HKU/profiling_shared_1_Omicron20/mle_mid.rds"
# output_dir <- "tmp_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(spatPomp))
setwd(here::here())
variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others")
num_of_strains <- length(variants)

source("scripts/model_fitting/helper/fixed_parameters.R")
source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")

check_redo <- FALSE
name_fit <- paste0(output_dir, "/", search_level, "_", partition, "_", n_cores, "_", task_id, ".rds")
if(check_redo){file.remove(name_fit)}

if(search_level == "low"){
  search_profiles <- list(
    # for iubf
    Nubf = 4,
    Nrep_per_param = 32,
    Nparam = 32,
    nbhd = fn_nbhd,
    prop = 0.5,
    cooling.type = "geometric",
    cooling.fraction.50 = 0.5,
    # for ibpf
    NREPS = 32,
    NP_EVAL = 500,
    NREPS_EVAL = 3,
    NBPF = 8,
    NP = 500,
    COOLING = 0.5,
    KEEP_TRACES = FALSE,
    KEEP_LIK_MAT = FALSE,
    SPAT_REGRESSION = 0.05,
    # for both
    scale_RW_SD_unit_specific = 1,
    scale_RW_SD_shared = 0.75,
    RW_SD = 0.02
  )
} else if(search_level == "mid"){
  search_profiles <- list(
    # for iubf
    Nubf = 6,
    Nrep_per_param = 64,
    Nparam = 32,
    nbhd = fn_nbhd,
    prop = 0.5,
    cooling.type = "geometric",
    cooling.fraction.50 = 0.5,
    # for ibpf
    NREPS = 32,
    NP_EVAL = 100,
    NREPS_EVAL = 4,
    NBPF = 8, 
    NP = 1000,
    COOLING = 0.5,
    KEEP_TRACES = FALSE,
    KEEP_LIK_MAT = FALSE,
    SPAT_REGRESSION = 0.05,
    # for both
    scale_RW_SD_unit_specific = 1,
    scale_RW_SD_shared = 0.75,
    RW_SD = 0.015
  )
} else if(search_level == "high"){
  search_profiles <- list(
    # for iubf
    Nubf = 1,
    Nrep_per_param = 128,
    Nparam = 32,
    nbhd = fn_nbhd,
    prop = 0.5,
    cooling.type = "geometric",
    cooling.fraction.50 = 0.5,
    # for ibpf
    NREPS = 32,
    NP_EVAL = 1500,
    NREPS_EVAL = 6,
    NBPF = 10, 
    NP = 1500,
    COOLING = 0.5,
    KEEP_TRACES = FALSE,
    KEEP_LIK_MAT = FALSE,
    SPAT_REGRESSION = 0.05,
    # for both
    scale_RW_SD_unit_specific = 0.75,
    scale_RW_SD_shared = 1,
    RW_SD = 0.015
  )
} else {
  stop("search_level must be one of 'low', 'mid', or 'high'.")
}

model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
mod_params <- readRDS(input_mle)
source("scripts/model_simulation/helper/put_params_into_model.R")
model_spatpomp <- put_params_into_model(model_spatpomp, mod_params)

## model fitting ##
source("scripts/model_fitting/helper/parameter_names_est.R")
profile_local <- function(
  model_spat = model_spatpomp,
  ncores,
  search,
  local_para_type,
  optimizer
){
  stopifnot(local_para_type %in% c("unit", "shared", "all"))
  ## Create a list to save all of the results.
  results <- list()
  n_of_units <- length(unit_names(model_spat))
  # Add unit numbers to each parameter
  est_param_names <- c(shared_param_names_est, unit_specific_names_est)
  est_param_names_expanded <- paste0(rep(est_param_names, each = n_of_units), seq_len(n_of_units))
  stopifnot(all(sort(est_param_names_expanded)%in%sort(names(coef(model_spat)))))

  # create rw.sd for search 
  rw.sd_s <- rep(list(search$RW_SD),times=length(est_param_names_expanded))
  names(rw.sd_s) <- est_param_names_expanded

  ## Whether to adjust the scale of rw.sd for the unit-specific / shared parameters
  rw.sd_s[grepl("^beta_naught",names(rw.sd_s))] <- search$RW_SD * search$scale_RW_SD_unit_specific # unit-specific
  rw.sd_s[!grepl("^beta_naught",names(rw.sd_s))] <- search$RW_SD * search$scale_RW_SD_shared # shared

  if(local_para_type == "unit"){
    ## If only fitting beta_naught, set the rw.sd for the other parameters to 0
    rw.sd_s[!grepl("beta_naught", names(rw.sd_s))] <- 0
  }
  if(local_para_type == "shared"){
    ## If only fitting shared parameters, set the rw.sd for the beta_naught parameters to 0
    rw.sd_s[grepl("beta_naught", names(rw.sd_s))] <- 0
  }

  rw.sd_s <- do.call(rw_sd, rw.sd_s)

  doParallel::registerDoParallel(ncores)
  doRNG::registerDoRNG(random_seed)

  if(optimizer == "iubf"){
    source("scripts/model_fitting/profiling_HKU/helper/fn_iubf.R", local = TRUE)
  } else if(optimizer == "ibpf"){
    source("scripts/model_fitting/profiling_HKU/helper/fn_ibpf.R", local = TRUE)
  } else {
    stop("optimizer must be one of 'iubf' or 'ibpf'.")
  }

  return(search_results)
}

set.seed(random_seed)
Omicron_profile_local <- bake(
  file = name_fit, {
    profile_local(
      model_spat = model_spatpomp,
      ncores = n_cores,
      search = search_profiles,
      local_para_type = local_param_type,
      optimizer = optimizer
    )
  }
)

