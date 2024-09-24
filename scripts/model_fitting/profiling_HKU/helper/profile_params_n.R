#!/usr/bin/Rscript
# source(here::here("scripts/data_processing/install_prerequisite.R"))
suppressPackageStartupMessages(library(spatPomp))
suppressPackageStartupMessages(library(doParallel))

task_id <- as.numeric(Sys.getenv("TASK_ID"))
check_redo <- FALSE

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop("input_mle, param_type, round_n, n_cores, search_level, estimator, and model_name must be supplied.", call.=FALSE)
} else {
  input_mle <- args[1]
  if(input_mle != "NA"){
    stopifnot(file.exists(input_mle))
  }
  
  param_type <- args[2]
  
  round_n <- as.numeric(args[3])
  stopifnot(!is.na(round_n))
  stopifnot(round(round_n) > 0)
  
  n_cores <- as.numeric(args[4])
  stopifnot(!is.na(n_cores))
  stopifnot(n_cores > 0)

  search_level <- args[5]

  estimator <- args[6]

  model_name <- args[7]
}
print(paste0("task_id: ", task_id, "; input_mle: ", input_mle, "; param_type: ", param_type, "; round_n: ", round_n, "; n_cores: ", n_cores, "; search_level: ", search_level, "; estimator: ", estimator))

source("scripts/model_fitting/helper/parameter_names_est.R") # load shared_param_names_est, unit_specific_names_est
if(param_type == "unit"){
  profile_para <- c(paste0(rep(unit_specific_names_est, each = 29), seq_len(29)))[task_id]
} else if(param_type == "shared"){
  profile_para <- shared_param_names_est[task_id]
} else if(param_type == "all"){
  profile_para <- c(paste0(rep(unit_specific_names_est, each = 29), seq_len(29)), shared_param_names_est)[task_id]
}

stopifnot(!is.na(profile_para))

output_dir <- paste0("results/model_data/profiling_HKU/", model_name, "/profiling_", param_type, "_", round_n, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read data
model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
source("scripts/model_simulation/helper/put_params_into_model.R")
source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")

# set up the search profiles
if(search_level == "low"){
  search_profiles <- list(
    # for ubf
    NREPS = 32,
    # for bpf
    NP_EVAL = 500,
    N_EVAL = 3,
    # for all
    N_INTERVAL = 10
  )
} else if(search_level == "mid"){
  search_profiles <- list(
    # for ubf
    NREPS = 64,
    # for bpf
    NP_EVAL = 800,
    N_EVAL = 4,
    # for all
    N_INTERVAL = 16
  )
} else if(search_level == "high"){
  search_profiles <- list(
    # for ubf
    NREPS = 128,
    # for bpf
    NP_EVAL = 1200,
    N_EVAL = 6,
    # for all
    N_INTERVAL = 32
  )
} else {
  stop("search_level must be low, mid or high.", call.=FALSE)
}

# Other settings
if(input_mle == "NA"){
  mle_local <- NA
} else {
  mle_local <- readRDS(input_mle)
}

# main function to profile the parameter
profile_Omicron20 <- function(
  model_spat = model_spatpomp,
  ncores,
  search,
  profile_para,
  fn_nbhd,
  estimator
){
  n_of_units <- length(unit_names(model_spat))

  ## Create a list to save all of the results.
  results <- list()

  # Add unit numbers to each parameter
  est_param_names <- c(shared_param_names_est, unit_specific_names_est)
  est_param_names_expanded <- paste0(rep(est_param_names, each = n_of_units), seq_len(n_of_units))
  stopifnot(all(sort(est_param_names_expanded)%in%sort(names(coef(model_spat)))))

  # Get lower bound for unit parameters 
  unit_lb <- rep(rep(NA, each=length(unit_specific_names_est)), each = n_of_units)
  names(unit_lb) <- paste0(rep(unit_specific_names_est, each = n_of_units), seq_len(n_of_units))
  unit_lb <- coef(model_spat)[names(unit_lb)]

  # Get upper bound for unit parameters 
  unit_ub <- rep(rep(NA, each=length(unit_specific_names_est)), each = n_of_units)
  names(unit_ub) <- paste0(rep(unit_specific_names_est, each = n_of_units), seq_len(n_of_units))
  unit_ub <- coef(model_spat)[names(unit_ub)]

  # Get lower bound for shared parameters 
  shared_lb <- coef(model_spat)[paste0(shared_param_names_est, 1)]
  names(shared_lb) <- shared_param_names_est

  # Get upper bound for shared parameters 
  shared_ub <- coef(model_spat)[paste0(shared_param_names_est, 1)]
  names(shared_ub) <- shared_param_names_est

  # Create data.frame with random unit parameters
  guesses_unit <- pomp::runif_design(
    lower = unit_lb,
    upper = unit_ub,
    nseq = search$N_INTERVAL
  )
  # guesses_unit[names(unit_ub)] <- mle[names(unit_ub)]
  if(profile_para %in% names(guesses_unit)){
    guesses_unit[profile_para] <- seq(from = 0.01, to = 0.99, length.out = search$N_INTERVAL)
  }
  
  # Create data.frame with random shared parameters
  guesses_shared <- pomp::runif_design(
    lower = shared_lb,
    upper = shared_ub,
    nseq = search$N_INTERVAL
  )
  # guesses_shared[names(shared_ub)] <- mle[names(shared_ub)]
  if(profile_para %in% names(guesses_shared)){
    guesses_shared[profile_para] <- seq(from = 0.01, to = 0.99, length.out = search$N_INTERVAL)
  }

  guesses_shared <- guesses_shared[, rep(1:ncol(guesses_shared), each = n_of_units)]
  names(guesses_shared) <- paste0(rep(shared_param_names_est, each = n_of_units), seq_len(n_of_units))
  
  # Combine the unit and shared parameters
  guesses <- cbind(guesses_shared, guesses_unit)

  # We need to add fixed parameters
  all_params <- coef(model_spat)
  fixed_params <- all_params[!names(all_params) %in% colnames(guesses)]
  fixed_mat <- matrix(
    rep(fixed_params, search$N_INTERVAL),
    byrow = TRUE, nrow = search$N_INTERVAL
  )
  colnames(fixed_mat) <- names(all_params[!names(all_params) %in% colnames(guesses)])

  # Combine estimated and fixed parameters, and reorder based on original order.
  guesses_all <- cbind(guesses, fixed_mat)[, names(coef(model_spat))]
  stopifnot(all(sort(names(guesses_all)) == sort(names(coef(model_spat)))))

  # Memory clean-up
  rm(guesses_unit, guesses_shared, fixed_mat, shared_lb,
     shared_ub, unit_lb, unit_ub,  all_params, fixed_params)
  gc()

  if(estimator == "ubf"){
    source("scripts/model_fitting/profiling_HKU/helper/fn_ubf.R", local = TRUE)
  } else if(estimator == "bpf"){
    source("scripts/model_fitting/profiling_HKU/helper/fn_bpf.R", local = TRUE)
  } else {
    stop("estimator must be one of 'ubf' or 'bpf'.")
  }

  return(search_results)
}

model_spatpomp <- put_params_into_model(model_spatpomp, mle_local)

# To be run on multiple nodes with HPC
set.seed(2023)
name_fit <- paste0(output_dir, "/profile_", profile_para, "_", search_level, ".rds")
if(check_redo){file.remove(name_fit)}
Omicron_profile <- pomp::bake(
  file = name_fit, {
    profile_Omicron20(
      model_spat = model_spatpomp,
      ncores = n_cores,
      search = search_profiles,
      profile_para = profile_para,
      fn_nbhd = fn_nbhd,
      estimator = estimator
    )
  }
)
