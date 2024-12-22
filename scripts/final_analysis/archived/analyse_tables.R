library(tidyverse)

# 95% CI for the fitted parameters
final_local_fit_path <- "results/model_data/profiling_HKU/Omicron20/profiling_unit_11_local_1"
files_local_fits <- list.files(final_local_fit_path, pattern = "local_iubf_[aic]", full.names = T)

source("scripts/model_fitting/helper/parameter_names_est.R") # load shared_param_names_est, unit_specific_names_est
para_names_fitting <- c(paste0(rep(unit_specific_names_est, each = 29), seq_len(29)), paste0(shared_param_names_est, 1))

df_paramMatrix_all <- lapply(files_local_fits, function(x){
  # x=files_local_fits[1]
  local_iubf <- readRDS(x)
  iubf_paramMatrix <- local_iubf@paramMatrix
  as_tibble(t(iubf_paramMatrix)) %>% select(all_of(para_names_fitting))
})
df_paramMatrix_all <- bind_rows(df_paramMatrix_all)

df_paramMatrix_CI <- df_paramMatrix_all %>% apply(2, quantile, c(0.025, 0.975)) %>% as_tibble() 
df_out <- t(df_paramMatrix_CI) %>% as_tibble()
names(df_out) <- c("2.5%", "97.5%")
df_out$parameter <- colnames(df_paramMatrix_CI)
writexl::write_xlsx(df_out, paste0(final_local_fit_path, "/local_iubf_CI.xlsx"))

source("scripts/model_fitting/helper/tranform_log_params.R")
df_paramMatrix_CI_transformed <- df_paramMatrix_CI

lapply(seq_len(ncol(df_paramMatrix_CI)), function(i){
  para_i <- names(df_paramMatrix_CI)[i]
  if(grepl("beta_naught", para_i)){
    para_i <- "beta_naught"
  } else{
    para_i <- gsub("\\d$", "", para_i)
  }

  df_paramMatrix_CI_transformed[[i]] <<- k_values[para_i] * df_paramMatrix_CI[[i]] + b_values[para_i]
  NULL
})

df_out_transformed <- t(df_paramMatrix_CI_transformed) %>% as_tibble()
names(df_out_transformed) <- c("2.5%", "97.5%")
df_out_transformed$parameter <- colnames(df_paramMatrix_CI_transformed)
writexl::write_xlsx(df_out_transformed, paste0(final_local_fit_path, "/local_iubf_CI_transformed.xlsx"))
