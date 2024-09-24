transform_log_params <- function(data, params_to_est, k_values, b_values){
  stopifnot(length(k_values) == length(b_values))
  stopifnot(length(k_values) == length(params_to_est))

  for(i in seq_along(params_to_est)){
    data[[params_to_est[i]]] <- k_values[i] * data[[params_to_est[i]]] + b_values[i]
  }
  return(data)
}

source("scripts/model_fitting/helper/parameter_names_est.R") # load shared_param_names_est, unit_specific_names_est

params_to_est <- c(unit_specific_names_est, shared_param_names_est)

# translate using y=kx+b
k_values <- c(
  1.99*365, # beta_naught, [0.01*365, 2*365], mid=1*365
  3.5, # eta_two, [1.5, 5], mid=3.25 (Ref: https://pubmed.ncbi.nlm.nih.gov/35458551/)
  # 3.5, # eta_three, [1.5, 5], mid=3.25 (Ref: https://pubmed.ncbi.nlm.nih.gov/35458551/)
  8.5, # eta_three, [1.5, 10], mid=5.75 (Ref: https://pubmed.ncbi.nlm.nih.gov/35458551/)
  1.5, # eta_four, [0.5, 2], mid=1.25
  365/3.6-365/4.4, # epsilon_four, [365/4.4, 365/3.6], mid=365/3.96
  0.28, # lambda_one, (Ref: https://www.nature.com/articles/s41467-022-31838-8, relative reduction of 71% to 99%, i.e. 0.01 to 0.29, mid=0.15)
  0.55, # lambda_two, (Ref: https://www.nature.com/articles/s41467-022-31838-8, relative reduction of 13% to 68%, i.e. 0.32 to 0.87, mid=0.595)
  0.23, # lambda_three, (slightly (5%) higher than lambda_two, i.e. 0.37 to 0.92, mid=0.645)
  0.28, # lambda_four, similar to lambda_one
  0.5, # iota_one, screening arrivals, [0.4, 0.9], mid=0.65
  0.5, # iota_two, quarantine arrivals, [0, 0.5], mid=0.25
  0.7, # iota_three, ban arrivals for some countries, [0.1, 0.8], mid=0.45
  0.3, # iota_four, ban arrivals for all countries, [0, 0.3], mid=0.15
  # 365/6-365/10, # gamma, [365/10, 365/6], mid=365/7.5
  365/3-365/10, # gamma, [365/10, 365/3], mid=365/5
  0.05-0.001, # sigma, [0.001, 0.05], mid=0.0255
  1, # tau_cases, [0, 1], mid=0.5
  1, # tau_deaths, [0, 1], mid=0.5
  1, # tau_seq_all, [0, 1], mid=0.5
  1, # tau_seq_unit, [0, 1], mid=0.5
  42, # day_BA_one, [3, 45], mid=24 (2022-10-09)
  65 # day_BA_two, [10, 75], mid=42.5 (2022-10-27)
)
names(k_values) <- c("beta_naught", "eta_two", "eta_three", "eta_four", "epsilon_four", "lambda_one", "lambda_two", "lambda_three", "lambda_four", "iota_one", "iota_two", "iota_three", "iota_four", "gamma", "sigma", "tau_cases", "tau_deaths", "tau_seq_all", "tau_seq_unit", "day_BA_one", "day_BA_two")
k_values <- k_values[match(params_to_est, names(k_values))]

b_values <- c(
  0.01*365, # beta_naught
  1.5, # eta_two
  1.5, # eta_three
  0.5, # eta_four
  365/4.4, # epsilon_four
  0.01, # lambda_one
  0.32, # lambda_two
  0.37, # lambda_three
  0.01, # lambda_four
  0.4, # iota_one
  0, # iota_two
  0.1, # iota_three
  0, # iota_four
  365/10, # gamma
  0.001, # sigma
  0, # tau_cases
  0, # tau_deaths
  0, # tau_seq_all
  0, # tau_seq_unit
  3, # day_BA_one
  10 # day_BA_two
)
names(b_values) <- c("beta_naught", "eta_two", "eta_three", "eta_four", "epsilon_four", "lambda_one", "lambda_two", "lambda_three", "lambda_four", "iota_one", "iota_two", "iota_three", "iota_four", "gamma", "sigma", "tau_cases", "tau_deaths", "tau_seq_all", "tau_seq_unit", "day_BA_one", "day_BA_two")
b_values <- b_values[match(params_to_est, names(b_values))]
