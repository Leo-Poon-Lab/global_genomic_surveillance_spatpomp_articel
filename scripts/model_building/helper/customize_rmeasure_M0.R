get_code_rmeasure_M0 <- function(U, num_of_strains=4){
  code_rmeasure_M0 <- paste0(readLines("scripts/model_building/helper/code_rmeasure_M0.C"), collapse = "\n")

  code_rmeasure_C_k_i_o <- apply(expand.grid(seq_len(num_of_strains), seq_len(U)-1), 1, function(x){
    name_x <- paste0("C_", x[1], "_i_sequenced_origin_", x[2], "_new[u]")
    name_y <- paste0("C_", x[1], "_i_o_", x[2], "_[u]")

    paste0("m = ", name_x, ";
      if (m > 0.0) {
        v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
        ", name_y, "= nearbyint(rnorm(m, sqrt(v)));
      } else {
        ", name_y, " = 0.0;
      }")
    }) %>% paste0(collapse="\n")

  code_rmeasure_M0 <- gsub("code_rmeasure_C_k_i_o_TO_BE_REPLACED", code_rmeasure_C_k_i_o, code_rmeasure_M0, fixed=T)

  # # Deal with parameters transformation
  # source("scripts/model_fitting/helper/tranform_log_params.R")
  # id_params_to_est_measure <- grepl("tau_", params_to_est)
  # code_initialize_transformed_params <- paste0(paste0("double ", params_to_est[id_params_to_est_measure], "_u=0;"), collapse=" ")
  # code_rmeasure_M0 <- gsub("initialize_transformed_params_TO_BE_REPLACED", code_initialize_transformed_params, code_rmeasure_M0, fixed=T)

  # code_transform_parameters <- paste0(paste0(params_to_est[id_params_to_est_measure], "_u = ", round(k_values[id_params_to_est_measure],3), " * (", params_to_est[id_params_to_est_measure], "[u*", params_to_est[id_params_to_est_measure], "_unit", "]) + ", round(b_values[id_params_to_est_measure],3), ";"), collapse=" ")
  # code_rmeasure_M0 <- gsub("code_transform_parameters_TO_BE_REPLACED", code_transform_parameters, code_rmeasure_M0, fixed=T)

  return(code_rmeasure_M0)
}
