#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("input/output dir, and model_name must be supplied.", call.=FALSE)
} else {
  input_dir <- args[1]
  output_dir <- input_dir
  
  model_name <- args[2]
}
# input_dir <- "results/model_data/profiling_HKU/Omicron20/profiling_all_1"
# input_dir <- "results/model_data/profiling_HKU/Omicron20/profiling_shared_10/"
# output_dir <- input_dir

check_mle_in_CI <- TRUE # in the latter part of the fitting, we switch this to TRUE

files_rds <- list.files(input_dir, full.names = TRUE, pattern = "profile_.+\\.rds$")

# ggsave may not run on servers lacking X11 support
library(tidyverse)
library(pomp)
source("scripts/model_fitting/helper/parameter_names_est.R")

params_all <- c(shared_param_names_est, paste0(rep(unit_specific_names_est, each = 29), seq_len(29)))

mle <- as.data.frame(matrix(NA, nrow = 1, ncol = length(params_all))) # the peak of the curve obtained using mcap
names(mle) <- params_all
mlev1_max <- mle
loglik_mle <- mle
loglik_mlev1_max <- mle
mlev2 <- as.data.frame(matrix(NA, nrow = 0, ncol = 757)) # the one with the highest log-likelihood, may not be the peak
check_CI_valid <- rep(NA, length(files_rds))
df_mcap <- tibble()

if(all(grepl("_low.rds$", files_rds))){
  search_level <- "low"
  Nreps <- 32
} else if(all(grepl("_mid.rds$", files_rds))){
  search_level <- "mid"
  Nreps <- 64
} else if(all(grepl("_high.rds$", files_rds))){
  search_level <- "high"
  Nreps <- 128
} else{
  stop("The files in the input directory should be all low, mid or high.")
}

for (i in seq_along(files_rds)) {
  file_rds_i <- files_rds[i]
  profile_para <- gsub(paste0("_", search_level, ".rds"), "", gsub("profile_", "", basename(file_rds_i)), fixed=T)

  profile_data <- readRDS(file_rds_i)
  if("bpf_time" %in% names(profile_data)){
    estimator <- "bpf"
  } else if ("ubf_time" %in% names(profile_data)){
    estimator <- "ubf"
  } else{
    stop("The estimator is not recognized.")
  }
  profile_data <- cbind(profile_data$logLiks,profile_data$params)
  if(any(names(profile_data)=="ll")){
    names(profile_data)[names(profile_data)=="ll"] <- "logLik"
  }
  names(mlev2) <- names(profile_data)

  if(profile_para %in% shared_param_names_est){
    profile_para1 <- paste0(profile_para,"1")
  }
  if(profile_para %in% paste0(rep(unit_specific_names_est, each = 29), seq_len(29))){
    profile_para1 <- profile_para
  }
  prof_res <- profile_data %>%
    group_by(.data[[profile_para1]]) %>%
    summarize(logLik = max(logLik))
  mcap_results <- pomp::mcap(prof_res$logLik, prof_res[[profile_para1]])
  if(is.na(mcap_results$se_stat)){
    check_CI_valid[i] <- FALSE
  } else{
    check_CI_valid[i] <- TRUE
  }
  names(check_CI_valid)[i] <- profile_para

  df_mcap <- rbind(df_mcap, tibble(profile_para = profile_para, mle=mcap_results$mle, CI_low=min(mcap_results$ci), CI_high=max(mcap_results$ci), delta=mcap_results$delta, se_stat=mcap_results$se_stat, se_mc=mcap_results$se_mc, se=mcap_results$se))

  if(check_mle_in_CI){ # constrain the MLEv2 to be within the CI
    id_filtered <- which(profile_data[[profile_para1]]>=sort(mcap_results$ci)[1] & profile_data[[profile_para1]]<=sort(mcap_results$ci)[2])
    if(length(id_filtered)<1){ # if the CI is too narrow, choose the closest point
      id_filtered <- which.min(abs(profile_data[[profile_para1]]-mcap_results$mle))
    }
    profile_data_filtered <- profile_data[id_filtered,]
    ## MLEv2
    mlev2 <- rbind(mlev2, profile_data_filtered[which.max(profile_data_filtered$logLik),])
    ## MLEv1
    mle[profile_para] <- profile_data_filtered[which.max(profile_data_filtered$logLik),][[profile_para1]] # we use this value for MLEv1, on and after profiling_unit_8
    loglik_mle[profile_para] <- max(profile_data_filtered$logLik)
  } else{
    ## MLEv2
    mlev2 <- rbind(mlev2, profile_data[which.max(profile_data$logLik),])
    ## MLEv1
    mle[profile_para] <- profile_data[which.max(profile_data$logLik),][[profile_para1]] # we use this value for MLEv1, on and after profiling_unit_8
    loglik_mle[profile_para] <- max(profile_data$logLik)
  }

  ## MLEv1_max
  mlev1_max[profile_para] <- profile_data[which.max(profile_data$logLik),][[profile_para1]]
  loglik_mlev1_max[profile_para] <- max(profile_data$logLik)
}

df_mcap <- df_mcap[match(params_all, df_mcap$profile_para),]
write_csv(df_mcap, paste0(output_dir, "/df_mcap_", search_level, ".csv"))

source("scripts/model_fitting/helper/tranform_log_params.R")
df_mcap_transformed <- df_mcap
df_mcap_transformed$mle_transformed <- NA
df_mcap_transformed$CI_low_transformed <- NA
df_mcap_transformed$CI_high_transformed <- NA
lapply(seq_len(nrow(df_mcap_transformed)), function(i){
  profile_para_i <- df_mcap_transformed$profile_para[i]
  if(grepl("beta_naught", profile_para_i)){
    profile_para_i <- "beta_naught"
  }
  df_mcap_transformed$mle_transformed[i] <<- k_values[profile_para_i] * df_mcap_transformed$mle[i] + b_values[profile_para_i]
  df_mcap_transformed$CI_low_transformed[i] <<- k_values[profile_para_i] * df_mcap_transformed$CI_low[i] + b_values[profile_para_i]
  df_mcap_transformed$CI_high_transformed[i] <<- k_values[profile_para_i] * df_mcap_transformed$CI_high[i] + b_values[profile_para_i]
  NULL
})
df_mcap_transformed <- df_mcap_transformed %>% select(profile_para, mle_transformed, CI_low_transformed, CI_high_transformed, delta, se_stat, se_mc, se, mle, CI_low, CI_high)
write_csv(df_mcap_transformed, paste0(output_dir, "/df_mcap_transformed_", search_level, ".csv"))

## MLEv2
mlev2 <- mlev2[which.max(mlev2$logLik),]
log_lik_mlev2 <- mlev2$logLik
mlev2 <- mlev2[,c(paste0(shared_param_names_est, 1), paste0(rep(unit_specific_names_est, each = 29), seq_len(29)))]
names(mlev2)[names(mlev2) %in% paste0(shared_param_names_est,1)] <- shared_param_names_est
saveRDS(mlev2, paste0(output_dir, "/mlev2_", search_level, ".rds"))

if(estimator == "ubf"){
  print(paste0("MLEv2 Loglik: ", log_lik_mlev2))
  writeLines(as.character(log_lik_mlev2), paste0(output_dir, paste0("/logLik_mlev2_", search_level, ".txt")))
} else if(estimator == "bpf"){
  writeLines(as.character(log_lik_mlev2), paste0(output_dir, paste0("/bpf_logLik_mlev2_", search_level, ".txt")))

  # needs to estimate the log-likelihood using ubf
  model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
  source("scripts/model_simulation/helper/put_params_into_model.R")
  model_spatpomp <- put_params_into_model(model_spatpomp, mlev2)
  source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")

  set.seed(2023)
  doParallel::registerDoParallel(64)
  doRNG::registerDoRNG(20231118)
  log_lik_mlev2 <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
  doParallel::stopImplicitCluster()
  
  print(paste0("MLEv2 Loglik: ", log_lik_mlev2))
  writeLines(as.character(log_lik_mlev2), paste0(output_dir, paste0("/logLik_mlev2_", search_level, ".txt")))
}

## MLEv1
names_mle_na <- names(mle)[is.na(mle)]
if(any(names_mle_na %in% names(mlev2))){
  names_direct_sub <- names_mle_na[names_mle_na %in% names(mlev2)]
  mle[names_direct_sub] <- mlev2[names_direct_sub]
}
names_mle_na <- names(mle)[is.na(mle)]
if(any(paste0(names_mle_na, 1) %in% names(mlev2))){
  names_sub <- names_mle_na[paste0(names_mle_na, 1) %in% names(mlev2)]
  mle[names_sub] <- mlev2[paste0(names_sub, 1)]
}
saveRDS(mle, paste0(output_dir, "/mlev1_", search_level, ".rds"))

model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
source("scripts/model_simulation/helper/put_params_into_model.R")
model_spatpomp <- put_params_into_model(model_spatpomp, mle)
source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")

set.seed(2023)
doParallel::registerDoParallel(64)
doRNG::registerDoRNG(20231118)
log_lik_mle <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
doParallel::stopImplicitCluster()

print(paste0("MLEv1 Loglik: ", log_lik_mle))
writeLines(as.character(log_lik_mle), paste0(output_dir, paste0("/logLik_mlev1_", search_level, ".txt")))

## MLEv3
### replace the MLEv1 with the MLEv2 if the MLEv1 is outside the CI
mlev3 <- mle
replacement_1 <- mlev2[names(check_CI_valid)[!check_CI_valid]]
if(length(replacement_1)>0){
  mlev3[names(replacement_1)] <- replacement_1
}
### replace the MLEv1 with the MLEv2 if the MLEv1 is too close to the boundary
replacement_2 <- mlev2[names(mle)[mle<1/64 | mle>(1-1/64)]]
if(length(replacement_2)>0){
  mlev3[names(replacement_2)] <- replacement_2
}
saveRDS(mlev3, paste0(output_dir, "/mlev3_", search_level, ".rds"))


if(all(mlev3==mle)){
  print("MLEv3 is the same as MLEv1.")
} else{
  print("MLEv3 is different from MLEv1.")
  model_spatpomp <- put_params_into_model(model_spatpomp, mlev3)
  set.seed(2023)
  doParallel::registerDoParallel(64)
  doRNG::registerDoRNG(20231118)
  log_lik_mlev3 <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
  doParallel::stopImplicitCluster()
  print(paste0("MLEv3 Loglik: ", log_lik_mlev3))
  writeLines(as.character(log_lik_mlev3), paste0(output_dir, paste0("/logLik_mlev3_", search_level, ".txt")))
}

mlev2_mod <- mlev2
for (i in seq_along(mlev1_max)){ # in profiling_unit_11
  para_name_i <- names(unlist(mlev1_max[i]))
  loglik_increase_i <- (loglik_mlev1_max[[para_name_i]] - loglik_mle[[para_name_i]])/abs(loglik_mle[[para_name_i]])

  if(!is.na(mlev1_max[i])){
    if(loglik_increase_i>0.01){
      print(paste0("The log-likelihood of ", para_name_i, " increased by ", loglik_increase_i, "."))
      mlev2_mod[[names(mlev1_max[i])]] <- unlist(mlev1_max[i])
    }
  }
}
# mlev2_mod$beta_naught7 <- 0.3125 # for profiling_unit_11

model_spatpomp <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
source("scripts/model_simulation/helper/put_params_into_model.R")
model_spatpomp <- put_params_into_model(model_spatpomp, mlev2_mod)
source("scripts/model_fitting/profiling_HKU/helper/fn_nbhd.R")

set.seed(2023)
doParallel::registerDoParallel(64)
doRNG::registerDoRNG(20231118)
log_lik_mlev2_mod <- logLik(abf(model_spatpomp, Np=1, Nrep=Nreps, nbhd=fn_nbhd)) # ubf is a special case when abf with Np=1
doParallel::stopImplicitCluster()

print(paste0("MLEv2_mod Loglik: ", log_lik_mlev2_mod))
writeLines(as.character(log_lik_mlev2_mod), paste0(output_dir, paste0("/logLik_mlev2_", search_level, "_mod.txt")))
saveRDS(mlev2_mod, paste0(output_dir, "/mlev2_", search_level, "_mod.rds"))


