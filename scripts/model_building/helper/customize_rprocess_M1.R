customize_rprocess_M1 <- function(code_rprocess, U, units_all, num_of_strains = 4, N_introdution){

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
    lapply(seq_len(num_of_strains), function(k) {
      gsub("k", k, routes)
    }) %>% unlist()
    
  }

  transitions_M1 <- c(
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

  num_of_transition <- length(transitions_M1)
  code_rprocess <- gsub("num_of_transition_TO_BE_REPLACED", num_of_transition, code_rprocess, fixed=T)

  code_initialize_transition_order <- paste0(paste0("int ", (transitions_M1), "=", seq_along(transitions_M1)-1, ";"), collapse = " ")
  code_rprocess <- gsub("initialize_transition_order_TO_BE_REPLACED", code_initialize_transition_order, code_rprocess, fixed=T)

  code_reset_rate_and_dN <- paste0("for(int i=0; i<", num_of_transition, "; i++){rate[i]=0; dN[i]=0;}", collapse = " ")
  code_rprocess <- gsub("reset_rate_and_dN_TO_BE_REPLACED", code_reset_rate_and_dN, code_rprocess, fixed=T)

  # codes for different levels of traveling events
  code_direct_travelers <- paste0(readLines("scripts/model_building/helper/code_rprocess_main_M1_direct_travelers.C"), collapse = "\n")
  code_rprocess <- gsub("code_travel_related_transitions_direct_TO_BE_REPLACED", code_direct_travelers, code_rprocess, fixed=T)

  code_transit_first <- paste0(readLines("scripts/model_building/helper/code_rprocess_main_M1_transit_first.C"), collapse = "\n")
  code_rprocess <- gsub("code_travel_related_transitions_transit_first_TO_BE_REPLACED", code_transit_first, code_rprocess, fixed=T)

  # return code_rprocess
  return(code_rprocess)
}
