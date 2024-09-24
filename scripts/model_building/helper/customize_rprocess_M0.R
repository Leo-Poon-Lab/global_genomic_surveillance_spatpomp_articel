customize_rprocess_M0 <- function(code_rprocess, U, units_all, num_of_strains = 4, N_introdution){

  # Update the aggregation of E_k_i and I_k_i compartments
  code_update_E_k_i <- lapply(seq_len(num_of_strains), function(k){
    paste0("E_", k, "_i[u] = ", paste0("E_", k, "_i_direct_in_origin_", seq_len(U)-1, "_[u]", collapse = " + "), " + ", paste0("E_", k, "_i_transit_first_origin_", seq_len(U)-1, "_[u]", collapse = " + "), ";")
  }) %>% paste0(collapse = "\n")
  code_rprocess <- gsub("code_update_E_k_i_TO_BE_REPLACED", code_update_E_k_i, code_rprocess, fixed=T)

  code_update_I_k_i <- lapply(seq_len(num_of_strains), function(k){
    paste0("I_", k, "_i[u] = ", paste0("I_", k, "_i_direct_in_origin_", seq_len(U)-1, "_[u]", collapse = " + "), " + ", paste0("I_", k, "_i_transit_first_origin_", seq_len(U)-1, "_[u]", collapse = " + "), ";")
  }) %>% paste0(collapse = "\n")
  code_rprocess <- gsub("code_update_I_k_i_TO_BE_REPLACED", code_update_I_k_i, code_rprocess, fixed=T)

  # The outbreak is introduced in ZA at first (South Africa)
  code_rprocess <- gsub("u==ZA_code", paste0("u==", which(units_all=="ZA")-1), code_rprocess, fixed=T)
  code_rprocess <- gsub("N_BA_one_introduced", N_introdution, code_rprocess, fixed=T)

  # Deal with parameters transformation
  source("scripts/model_fitting/helper/tranform_log_params.R")
  code_initialize_transformed_params <- paste0(paste0("double ", params_to_est, "_u=0;"), collapse=" ")
  code_rprocess <- gsub("initialize_transformed_params_TO_BE_REPLACED", code_initialize_transformed_params, code_rprocess, fixed=T)

  code_transform_parameters <- paste0(paste0(params_to_est, "_u = ", round(k_values,3), " * (", params_to_est, "[u*", params_to_est, "_unit", "]) + ", round(b_values,3), ";"), collapse=" ")
  code_rprocess <- gsub("code_transform_parameters_TO_BE_REPLACED", code_transform_parameters, code_rprocess, fixed=T)

  # Deal with transitions between compartments

  # Community related transitions, illustrated as below:
  #     > V        > T      > T > D
  #    /   \      /        /   /
  #   /     \>   /        /   /
  #  S  -->  E_k_c  -->  I_k_c  -->  R

  # Travel related transitions, illustrated as below:
  travel_levels <- c(0,1)
  # 1. Direct travel #
  # T                             T              > D
  #  \(T_Eki0)                     \(T_Iki0)    /(Iki0_D)
  #   \>             (Eki0_Iki0)    \>         /     (Iki0_R)
  #  E_k_i_direct_in      -->      I_k_i_direct_in      -->      R
  #     \                             \
  #      \(Eki0_P)                     \(Iki0_P)
  #       \>P                           \>P
  # 2. Transit first #
  # P                                  P                 > D
  #  \(P_Eki1)                          \(P_Iki1)       /(Iki1_D)
  #   \>                  (Eki1_Iki1)    \>            /      (Iki1_R)
  #   E_k_i_transit_first      -->      I_k_i_transit_first      -->      R

  # travel level 0
  routes_E_direct_in <- c("T_Eki0")
  routes_E_direct_out <- c("Eki0_P", "Eki0_Iki0")
  routes_I_direct_in <- c("T_Iki0")
  routes_I_direct_out <- c("Iki0_P", "Iki0_D", "Iki0_R")
  # travel level 1
  routes_E_transit_first_out <- c("Eki1_P", "Eki1_Iki1")
  routes_I_transit_first_out <- c("Iki1_P", "Iki1_D", "Iki1_R")

  enumerate_all_routes <- function(routes, num_of_strains, U){
    # enumerate for strains
    routes_strains <- lapply(seq_len(num_of_strains), function(k) {
      gsub("k", k, routes)
    }) %>% unlist()
    # enumerate for units
    paste0(paste0(routes_strains, "_origin_"), rep(seq_len(U)-1, each=num_of_strains*length(routes)))
  }

  transitions_M0 <- c(
    ## Before exposed
    "S_V", "S_E1c", "S_E2c", "S_E3c", "S_E4c",
    "V_E1c", "V_E2c", "V_E3c", "V_E4c",
    ## community related E to I, and community related E to T (direct travel)
    "E1c_T", "E1c_I1c",
    "E2c_T", "E2c_I2c",
    "E3c_T", "E3c_I3c",
    "E4c_T", "E4c_I4c",
    ## community related I to R/D, and community related I to T (direct travel)
    "I1c_T", "I1c_D", "I1c_R",
    "I2c_T", "I2c_D", "I2c_R",
    "I3c_T", "I3c_D", "I3c_R",
    "I4c_T", "I4c_D", "I4c_R",
    ## E_direct in
    enumerate_all_routes(routes_E_direct_in, num_of_strains, U),
    ## E_direct out
    enumerate_all_routes(routes_E_direct_out, num_of_strains, U),
    ## I_direct in
    enumerate_all_routes(routes_I_direct_in, num_of_strains, U),
    ## I_direct out
    enumerate_all_routes(routes_I_direct_out, num_of_strains, U),
    ## E_transit_first out
    enumerate_all_routes(routes_E_transit_first_out, num_of_strains, U),
    ## I_transit_first out
    enumerate_all_routes(routes_I_transit_first_out, num_of_strains, U),
    NULL
  )

  num_of_transition <- length(transitions_M0)
  code_rprocess <- gsub("num_of_transition_TO_BE_REPLACED", num_of_transition, code_rprocess, fixed=T)

  code_initialize_transition_order <- paste0(paste0("int ", (transitions_M0), "=", seq_along(transitions_M0)-1, ";"), collapse = " ")
  code_rprocess <- gsub("initialize_transition_order_TO_BE_REPLACED", code_initialize_transition_order, code_rprocess, fixed=T)

  code_reset_rate_and_dN <- paste0("for(int i=0; i<", num_of_transition, "; i++){rate[i]=0; dN[i]=0;}", collapse = " ")
  code_rprocess <- gsub("reset_rate_and_dN_TO_BE_REPLACED", code_reset_rate_and_dN, code_rprocess, fixed=T)

  # A function for replacing v and origin in the below example codes
  replace_origin_to_v <- function(code, v){
    code <- gsub("v==0", paste0("v==", v), code)
    code <- gsub("origin_0", paste0("origin_", v), code)
    code
  }
  # codes for different levels of traveling events for origin v (based on the example of `v==0`)
  code_direct_travelers_v0 <- paste0(readLines("scripts/model_building/helper/code_rprocess_main_M0_direct_travelers_v0.C"), collapse = "\n")
  code_travel_related_transitions_direct <- sapply(seq_len(U)-1, function(v){
    code_v <- replace_origin_to_v(code_direct_travelers_v0, v)
    code_v
  }) %>% paste0(collapse = "\n")
  code_rprocess <- gsub("code_travel_related_transitions_direct_TO_BE_REPLACED", code_travel_related_transitions_direct, code_rprocess, fixed=T)

  code_transit_first_v0 <- paste0(readLines("scripts/model_building/helper/code_rprocess_main_M0_transit_first_v0.C"), collapse = "\n")
  code_travel_related_transitions_transit_first <- sapply(seq_len(U)-1, function(v){
    code_v <- replace_origin_to_v(code_transit_first_v0, v)
    code_v
  }) %>% paste0(collapse = "\n")
  code_rprocess <- gsub("code_travel_related_transitions_transit_first_TO_BE_REPLACED", code_travel_related_transitions_transit_first, code_rprocess, fixed=T)

  # Transitions in R and D compartments
  code_dN_Ikin_R_o <- paste0("dN[I", seq_len(num_of_strains), "i", travel_levels, "_R_origin_", rep(seq_len(U)-1, each=num_of_strains*length(travel_levels)), "]") %>% paste0(collapse = " + ")
  code_rprocess <- gsub("code_dN_Ikin_R_o_TO_BE_REPLACED", code_dN_Ikin_R_o, code_rprocess, fixed=T)

  code_dN_Ikin_D_o <- paste0("dN[I", seq_len(num_of_strains), "i", travel_levels, "_D_origin_", rep(seq_len(U)-1, each=num_of_strains*length(travel_levels)), "]") %>% paste0(collapse = " + ")
  code_rprocess <- gsub("code_dN_Ikin_D_o_TO_BE_REPLACED", code_dN_Ikin_D_o, code_rprocess, fixed=T)

  # New imported cases related transitions
  code_C_k_i_infected_origin_o_new <- paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_infected_origin_"), rep(seq_len(U)-1, each=4)), "_new[u] = dN[I", seq_len(num_of_strains), "i0_R_origin_", rep(seq_len(U)-1, each=num_of_strains), "] + dN[I", seq_len(num_of_strains), "i1_R_origin_", rep(seq_len(U)-1, each=num_of_strains), "] + dN[I", seq_len(num_of_strains), "i0_D_origin_", rep(seq_len(U)-1, each=num_of_strains), "] + dN[I", seq_len(num_of_strains), "i1_D_origin_", rep(seq_len(U)-1, each=num_of_strains), "];") %>% paste0(collapse = " ")
  code_rprocess <- gsub("code_C_k_i_infected_origin_o_new_TO_BE_REPLACED", code_C_k_i_infected_origin_o_new, code_rprocess, fixed=T)

  code_C_k_i_detected_origin_o_new <- paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_detected_origin_"), rep(seq_len(U)-1, each=4)), "_new[u] = ", paste0(paste0("C_", seq_len(num_of_strains), "_i_infected_origin_"), rep(seq_len(U)-1, each=4)), "_new[u]*infection_detection[u];") %>% paste0(collapse = " ")
  code_rprocess <- gsub("code_C_k_i_detected_origin_o_new_TO_BE_REPLACED", code_C_k_i_detected_origin_o_new, code_rprocess, fixed=T)

  code_C_k_i_sequenced_origin_o_new <- paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_sequenced_origin_"), rep(seq_len(U)-1, each=4)), "_new[u] = ", paste0(paste0("C_", seq_len(num_of_strains), "_i_detected_origin_"), rep(seq_len(U)-1, each=4)), "_new[u]*detection_sequencing[u];") %>% paste0(collapse = " ")
  code_rprocess <- gsub("code_C_k_i_sequenced_origin_o_new_TO_BE_REPLACED", code_C_k_i_sequenced_origin_o_new, code_rprocess, fixed=T)

  code_C_1_all_sequenced_new <- paste0("C_1_all_sequenced_new[u] = C_1_c_sequenced_new[u] + ", paste0("C_1_i_detected_origin_", rep(seq_len(U)-1), "_new[u]", collapse = " + "))
  code_rprocess <- gsub("code_C_1_all_sequenced_new_TO_BE_REPLACED", code_C_1_all_sequenced_new, code_rprocess, fixed=T)
  code_C_2_all_sequenced_new <- paste0("C_2_all_sequenced_new[u] = C_2_c_sequenced_new[u] + ", paste0("C_2_i_detected_origin_", rep(seq_len(U)-1), "_new[u]", collapse = " + "))
  code_rprocess <- gsub("code_C_2_all_sequenced_new_TO_BE_REPLACED", code_C_2_all_sequenced_new, code_rprocess, fixed=T)
  code_C_3_all_sequenced_new <- paste0("C_3_all_sequenced_new[u] = C_3_c_sequenced_new[u] + ", paste0("C_3_i_detected_origin_", rep(seq_len(U)-1), "_new[u]", collapse = " + "))
  code_rprocess <- gsub("code_C_3_all_sequenced_new_TO_BE_REPLACED", code_C_3_all_sequenced_new, code_rprocess, fixed=T)
  code_C_4_all_sequenced_new <- paste0("C_4_all_sequenced_new[u] = C_4_c_sequenced_new[u] + ", paste0("C_4_i_detected_origin_", rep(seq_len(U)-1), "_new[u]", collapse = " + "))
  code_rprocess <- gsub("code_C_4_all_sequenced_new_TO_BE_REPLACED", code_C_4_all_sequenced_new, code_rprocess, fixed=T)

  code_C_k_i_detected_all <- paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_detected_origin_"), rep(seq_len(U)-1, each=4)), "_new[u]") %>% paste0(collapse = " + ")
  code_rprocess <- gsub("code_C_k_i_detected_all_TO_BE_REPLACED", code_C_k_i_detected_all, code_rprocess, fixed=T)

  # return code_rprocess
  return(code_rprocess)
}
