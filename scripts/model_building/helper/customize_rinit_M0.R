# Also refer to "scripts/model_building/helper/unit_state_names_M0.R"
## For traveling-related units, we assume that 14% of the traveling population is in the transit units
## For transit units, we assume that the first 80% of the population is in the first transit unit, and the rest is in the second transit unit
## For the initial condition, we assume the traveling population is uniformly distributed among the origin units

get_code_rinit_M0 <- function(U){
  code_rinit_M0 <- paste0(readLines("scripts/model_building/helper/code_rinit_M0.C"), collapse = "\n")

  code_rinit_M0 <- gsub("E_by_direct_in_TO_BE_REPLACED", paste0(paste0(
    "E_1_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(E_1_i_init[u]*0.86/", U, ");\n",
    "E_2_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(E_2_i_init[u]*0.86/", U, ");\n",
    "E_3_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(E_3_i_init[u]*0.86/", U, ");\n",
    "E_4_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(E_4_i_init[u]*0.86/", U, ");\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("I_by_direct_in_TO_BE_REPLACED", paste0(paste0(
    "I_1_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(I_1_i_init[u]*0.86/", U, ");\n",
    "I_2_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(I_2_i_init[u]*0.86/", U, ");\n",
    "I_3_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(I_3_i_init[u]*0.86/", U, ");\n",
    "I_4_i_direct_in_origin_", seq_len(U)-1, "_[u] = nearbyint(I_4_i_init[u]*0.86/", U, ");\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("E_by_transit_first_TO_BE_REPLACED", paste0(paste0(
    "E_1_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(E_1_i_init[u]*0.14/", U, ");\n",
    "E_2_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(E_2_i_init[u]*0.14/", U, ");\n",
    "E_3_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(E_3_i_init[u]*0.14/", U, ");\n",
    "E_4_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(E_4_i_init[u]*0.14/", U, ");\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("I_by_transit_first_TO_BE_REPLACED", paste0(paste0(
    "I_1_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(I_1_i_init[u]*0.14/", U, ");\n",
    "I_2_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(I_2_i_init[u]*0.14/", U, ");\n",
    "I_3_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(I_3_i_init[u]*0.14/", U, ");\n",
    "I_4_i_transit_first_origin_", seq_len(U)-1, "_[u] = nearbyint(I_4_i_init[u]*0.14/", U, ");\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("New_cases_infected_TO_BE_REPLACED", paste0(paste0(
    "C_1_i_infected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_2_i_infected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_3_i_infected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_4_i_infected_origin_", seq_len(U)-1, "_new[u] = 0;\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("New_cases_detected_TO_BE_REPLACED", paste0(paste0(
    "C_1_i_detected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_2_i_detected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_3_i_detected_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_4_i_detected_origin_", seq_len(U)-1, "_new[u] = 0;\n"
  ), collapse = ""), code_rinit_M0)

  code_rinit_M0 <- gsub("New_cases_sequenced_TO_BE_REPLACED", paste0(paste0(
    "C_1_i_sequenced_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_2_i_sequenced_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_3_i_sequenced_origin_", seq_len(U)-1, "_new[u] = 0;\n",
    "C_4_i_sequenced_origin_", seq_len(U)-1, "_new[u] = 0;\n"
  ), collapse = ""), code_rinit_M0)

  return(code_rinit_M0)
}

