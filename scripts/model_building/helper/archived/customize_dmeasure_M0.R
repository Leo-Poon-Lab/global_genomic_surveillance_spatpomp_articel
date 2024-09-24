get_code_dmeasure_M0 <- function(U, num_of_strains=4, num_of_units=29){
  code_dmeasure_M0 <- paste0(readLines("scripts/model_building/helper/code_dmeasure_M0.C"), collapse = "\n")

  code_dmeasure_C_k_i_o <- apply(expand.grid(seq_len(num_of_strains), seq_len(U)-1), 1, function(x){
    name_x <- paste0("C_", x[1], "_i_sequenced_origin_", x[2], "_new[u]")
    name_y <- paste0("C_", x[1], "_i_o_", x[2], "_[u]")

    paste0("m = ", name_x, " > 0 ? ", name_x, " : 0;
    v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
    if (", name_y, " > tol) {
      lik += (log(pnorm(", name_y, "+0.5,m,sqrt(v),1,0)-pnorm(", name_y, "-0.5,m,sqrt(v),1,0) + tol))/N_units_TO_BE_REPLACED;
    }"
    )
  }) %>% paste0(collapse="\n")

  code_dmeasure_M0 <- gsub("code_dmeasure_C_k_i_o_TO_BE_REPLACED", code_dmeasure_C_k_i_o, code_dmeasure_M0, fixed=T)

  code_dmeasure_M0 <- gsub("N_units_TO_BE_REPLACED", num_of_units, code_dmeasure_M0, fixed=T)

  return(code_dmeasure_M0)
}
