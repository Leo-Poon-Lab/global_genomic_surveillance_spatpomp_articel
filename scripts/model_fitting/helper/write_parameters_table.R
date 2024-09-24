library(tidyverse)

write_parameter_tables <- function(
  model_name,
  output_path,
  mod_params_fitted,
  country_under_investigation,
  n_of_units,
  cross_check_table
){
  # Get the fitted parameters
  source("scripts/model_fitting/helper/parameter_names_est.R")
  mod_params_u <- matrix(nrow = length(unit_specific_names_est), ncol = n_of_units)
  rownames(mod_params_u) <- unit_specific_names_est
  mod_params_s <- numeric(length(shared_param_names_est))
  names(mod_params_s) <- shared_param_names_est
  for (p in shared_param_names_est) {
    if(p %in% names(mod_params_fitted)){
      mod_params_s[p] <- unlist(mod_params_fitted[p])
    } else {
      mod_params_s[p] <- unlist(mod_params_fitted[paste0(p, 1)])
    }
  }
  for (p in unit_specific_names_est) {
    pat <- paste0("^", p, "[[:digit:]]+$")
    mod_params_u[p, ] <- unlist(mod_params_fitted[grepl(pat, names(mod_params_fitted))])
  }

  df_mod_params <- bind_cols(as_tibble(t(mod_params_u)), as_tibble(t(mod_params_s)))
  df_mod_params$code <- country_under_investigation
  df_mod_params <- left_join(df_mod_params, cross_check_table %>% select(code, loc_name), "code")
  writexl::write_xlsx(df_mod_params, paste0(output_path, "/params_", model_name, ".xlsx"))

  source("scripts/model_fitting/helper/tranform_log_params.R")
  df_mod_params_transformed <- transform_log_params(df_mod_params, params_to_est, k_values, b_values)
  writexl::write_xlsx(df_mod_params_transformed, paste0(output_path, "/params_transformed_", model_name, ".xlsx"))
}
