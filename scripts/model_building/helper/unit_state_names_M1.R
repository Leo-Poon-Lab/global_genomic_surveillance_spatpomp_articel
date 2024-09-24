get_unit_state_names <- function(num_of_strains, U){
  unit_state_names <- c(
    "S", # susceptible
    "V", # vaccinated or immunized by past infection
    "R", # total recovered
    "D", # total deceased
    
    # Aggregated E and I community
    paste0("E_", seq_len(num_of_strains), "_c"), # community exposed
    paste0("I_", seq_len(num_of_strains), "_c"), # community infected
    
    # Aggregated E and I imported (income and stay)
    paste0("E_", seq_len(num_of_strains), "_i"), # imported exposed
    paste0("I_", seq_len(num_of_strains), "_i"), # imported infected

    # E and I separated by direct and transit
    # We only allow transiting for twice, avoiding infinite loops
    paste0(paste0("E_", seq_len(num_of_strains), "_i_direct_in")), # imported exposed
    paste0(paste0("I_", seq_len(num_of_strains), "_i_direct_in")), # imported infected
    paste0(paste0("E_", seq_len(num_of_strains), "_i_transit_first")), # imported exposed transit
    paste0(paste0("I_", seq_len(num_of_strains), "_i_transit_first")), # imported infected transit  

    # New cases which will be reset at each round
    ## Community cases separated by infected and detected
    "C_1_c_infected_new", "C_2_c_infected_new", "C_3_c_infected_new", "C_4_c_infected_new", # community infected cases of each strain
    "C_1_c_detected_new", "C_2_c_detected_new", "C_3_c_detected_new", "C_4_c_detected_new", # community detected cases of each strain
    "C_1_c_sequenced_new", "C_2_c_sequenced_new", "C_3_c_sequenced_new", "C_4_c_sequenced_new", # community sequenced cases of each strain, will be fitted, N=4

    ## Imported cases separated by infected, detected and sequenced
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_infected_new"))), # infected cases
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_detected_new"))), # detected imported cases
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_i_sequenced_new"))), # sequenced imported cases

    ## Add community and imported cases
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_all_infected_new"))), # infected cases
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_all_detected_new"))), # detected imported cases
    paste0(paste0(paste0("C_", seq_len(num_of_strains), "_all_sequenced_new"))), # sequenced imported cases, these will be fitted, N=4

    "C_reported_new", # new cases reported, this will be fitted, N=1
    "D_new", # new deceased, this will be fitted, N=1
    NULL
  )
}