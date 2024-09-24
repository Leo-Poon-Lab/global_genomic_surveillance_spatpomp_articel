# This script contains the names of the parameters *to be estimated* in the model fitting process.

shared_param_names_est = c( # shared parameters to be estimated
  "eta_two", # Relative increase/decrease in susceptibility due to strain 2 (e.g. Omicron BA.1)
  "eta_three", # Relative increase/decrease in susceptibility due to strain 3 (e.g. Omicron BA.2)
  "eta_four", # Relative increase/decrease in susceptibility due to strain 4 (e.g. Others)
  "epsilon_four", # Inverse of the latent/incubation period for strain 4 (e.g. Others)
  "lambda_one", # Relative decrease in infection rate due to effective immunity for strain 1 (e.g. Delta)
  "lambda_two", # Relative decrease in infection rate due to effective immunity for strain 2 (e.g. Omicron BA.1)
  "lambda_three", # Relative decrease in infection rate due to effective immunity for strain 3 (e.g. Omicron BA.2)
  "lambda_four", # Relative decrease in infection rate due to effective immunity for strain 4 (e.g. Others)
  "iota_one", # Relative decrease in transmission of imported cases due to International travel control level one - screening arrivals
  "iota_two", # Relative decrease in transmission of imported cases due to International travel control level two - quarantine arrivals
  "iota_three", # Relative decrease in transmission of imported cases due to International travel control level three - ban arrivals for some countries
  "iota_four", # Relative decrease in transmission of imported cases due to International travel control level four - ban arrivals for all countries
  "gamma", # Rate at which infected individuals become recovered or deceased
  "sigma", # intensity for Gamma white noise
  "tau_cases", # dispersion parameter for detected cases normal distribution observation process
  "tau_deaths", # dispersion parameter for deaths normal distribution observation process
  "tau_seq_all", # dispersion parameter for sequencing all cases normal distribution observation process
  "tau_seq_unit", # dispersion parameter for sequencing community/imported cases normal distribution observation process
  "day_BA_one", # the day when 1000 Omicron BA.1 cases was introduced into both E and I compartments of South Africa (N=1000 approximate to <3% of the total infection in SA on 2021-09-25, based on IHME estimates).
  "day_BA_two" # the day when 1000 Omicron BA.2 cases was introduced)
  )

unit_specific_names_est = c( # unit specific parameters to be estimated
  "beta_naught" # natural transmissibility (infection rate) in unit `u`
  )
