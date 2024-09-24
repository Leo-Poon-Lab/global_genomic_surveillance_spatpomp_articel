#!/usr/bin/Rscript

library(spatPomp)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("input_mle and model_name must be supplied.n", call.=FALSE)
} else {
  input_dir <- args[1]
  output_dir <- input_dir

  model_name <- args[2]
}
# input_dir <- "results/model_data/profiling_HKU/profiling_shared_1_local_4/"
# output_dir <- input_dir

# local search rds
files_rds <- c(
  list.files(input_dir, pattern = "^low_\\D+_\\d+_\\d+.rds$", full.names = TRUE),
  list.files(input_dir, pattern = "^mid_\\D+_\\d+_\\d+.rds$", full.names = TRUE),
  list.files(input_dir, pattern = "^high_\\D+_\\d+_\\d+.rds$", full.names = TRUE)
  )

if(all(grepl("^low_\\D+_\\d+_\\d+.rds$", basename(files_rds)))){
  search_level <- "low"
  Nreps <- 32
} else if(all(grepl("^mid_\\D+_\\d+_\\d+.rds$", basename(files_rds)))){
  search_level <- "mid"
  Nreps <- 64
} else if(all(grepl("^high_\\D+_\\d+_\\d+.rds$", basename(files_rds)))){
  search_level <- "high"
  Nreps <- 128
} else{
  stop("The files in the input directory should be all low, mid or high.")
}

source("scripts/model_fitting/helper/parameter_names_est.R")
files_local_prams <- list.files(input_dir, pattern = "^local_.+.rds$", full.names = TRUE)

if(all(grepl("local_iubf", basename(files_local_prams)))){
  # iubf
  mod_params_all <- lapply(files_rds, function(input_file){
    # input_file <- files_rds[1]
    mod_params_all <- readRDS(input_file)
    mod_params <- as.data.frame(t(mod_params_all$params))
    rownames(mod_params) <- NULL
    mod_params <- cbind(logLik=mod_params_all$logLiks, mod_params)
  })
  mod_params_all <- do.call(rbind, mod_params_all)
  mod_params_best <- mod_params_all[which.max(mod_params_all$logLik),]

  for (p in shared_param_names_est) {
    for(i in 1:29){
      # use the first unit's parameters for all units
      mod_params_best[paste0(p, i)] <- unlist(mod_params_best[paste0(p, 1)])
    }
  }
  saveRDS(mod_params_best, paste0(output_dir, "/mle_", search_level, ".rds"))

  if(grepl("results/model_data/profiling_HKU/profiling_shared_1_local_4", input_dir)){
    # special case for this run, needs to estimate the loglik again using ubf
    model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
    source("scripts/model_simulation/helper/put_params_into_model.R")
    model_spatpomp <- put_params_into_model(model_spatpomp, as.data.frame(mod_params_best))
    source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")
    
    set.seed(2023)
    doParallel::registerDoParallel(64)
    doRNG::registerDoRNG(20231118)
    log_lik_mle <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
    doParallel::stopImplicitCluster()

    print(paste0("logLog: ", log_lik_mle))
    writeLines(as.character(log_lik_mle), paste0(output_dir, "/logLik_mle_", search_level, ".txt"))
  } else {
    print(paste0("logLik: ", mod_params_best$logLik))
    writeLines(as.character(mod_params_best$logLik), paste0(output_dir, "/logLik_mle_", search_level, ".txt"))
  }

} else if(all(grepl("local_ibpf", basename(files_local_prams)))){
  # ibpf
  mod_params_all <- lapply(files_rds, function(input_file){
    mod_params_all <- readRDS(input_file)
    mod_params <- t(as.data.frame(mod_params_all$params[which.max(mod_params_all$logLik$logLik),]))
    rownames(mod_params) <- NULL
    mod_params <- cbind(mod_params_all$logLiks[which.max(mod_params_all$logLik$logLik),],mod_params)
  })
  mod_params_all <- do.call(rbind, mod_params_all)
  mod_params_best <- mod_params_all[which.max(mod_params_all$logLik),]

  saveRDS(mod_params_best, paste0(output_dir, "/mle_", search_level, ".rds"))
  writeLines(as.character(mod_params_best$logLik), paste0(output_dir, "/ibpf_logLik_mle_", search_level, ".txt"))

  # estimate the loglik using ubf
  model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
  source("scripts/model_simulation/helper/put_params_into_model.R")
  model_spatpomp <- put_params_into_model(model_spatpomp, as.data.frame(mod_params_best))

  source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")
  set.seed(2023)
  doParallel::registerDoParallel(64)
  doRNG::registerDoRNG(20231118)
  log_lik_mle <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
  doParallel::stopImplicitCluster()

  print(paste0("logLog: ", log_lik_mle))
  writeLines(as.character(log_lik_mle), paste0(output_dir, "/logLik_mle_", search_level, ".txt"))
} else {
  stop("The files in the input directory should be either local_iubf or local_ibpf.")
}
