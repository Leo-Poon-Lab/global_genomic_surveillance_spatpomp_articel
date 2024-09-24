# fixed parameters
fixed_shared_params_values <- c(
  eta_one = 1, # Relative increase/decrease in transmissibility due to strain 1 (e.g. Delta) set to 1
  epsilon_one = 1 / (4.4/365), # Inverse of the latent/incubation period for strain 1 (e.g. Delta)
  epsilon_two = 1 / (3.6/365), # Inverse of the latent/incubation period for strain 2 (e.g. Omicron BA.1)
  epsilon_three = 1 / (3.4/365), # Inverse of the latent/incubation period for strain 3 (e.g. Omicron BA.2)
  # epsilon4 = mean(c(1 / (4.4/365), 1 / (3.6/365))), # take the average of the latent/incubation periods for Others
  NULL
)

if(exists("data_fitting")) {
  fixed_unit_specific_params_values <- list(
    q = data_fitting$transit_rates # fraction of travelers that are transit passenger in unit `u`
  )
}
