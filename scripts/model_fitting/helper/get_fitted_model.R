get_fitted_model <- function(
  model_name,
  search_profile,
  archived_model = FALSE,
  write_parameter_tables = FALSE
){
  # Load data
  search_profiles <- readRDS(here::here("results/model_data/search_profiles.rds"))
  if(archived_model){
    model_fitted <- readRDS(paste0("results/model_data/archived/model_", model_name, ".rds"))
  } else {
    model_fitted <- readRDS(paste0("results/model_data/model_", model_name, ".rds"))
  }
  country_under_investigation <- model_fitted@unit_names

  cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% dplyr::filter(level==3)

  n_of_units = length(unit_names(model_fitted))
  source("scripts/model_fitting/helper/parameter_names_est.R")

  source(here::here("scripts/model_simulation/helper/load_best_params.R"))
  output_suffix <- get_output_suffix(search_profile=search_profile, archived_fit=FALSE, model_name=model_name)
  if(search_profile=="again"){
    model_fit <- readRDS(here::here(paste0("results/model_data/fit_", model_name, "_mid", output_suffix, ".rds")))
  } else {
    model_fit <- readRDS(here::here(paste0("results/model_data/fit_", model_name, output_suffix, ".rds")))
  }
  if(grepl("again", output_suffix)){names(model_fit) <- "search2"}
  if(search_profile=="low"){
    search2_profile <- search_profiles$search2_low
  } else if(search_profile=="mid"){
    search2_profile <- search_profiles$search2_mid
  } else if(search_profile=="high"){
    search2_profile <- search_profiles$search2_high
  } else if(search_profile=="again"){
    search2_profile <- search_profiles$search2
  } else {
    stop("search_profile must be one of low, mid, high, again")
  }

  mod_n_params <- model_fit[[1]]$mod_n_params
  print(paste0("The parameter size is ", mod_n_params/112, " times than that in the Hatti3 model."))

  # model_fit contains ll estimate for many parameter sets, so we want to find
  # the highest value
  mod_ll <- max(
    model_fit[[length(model_fit)]]$logLiks$logLik,
    na.rm = TRUE
  )
  print(paste0("The log likelihood is ", mod_ll, "."))

  # model 3 AIC
  mod_aic <- 2 * mod_n_params - 2 * mod_ll
  print(paste0("The AIC is ", mod_aic, "."))

  # Find the index for parameters of model that correspond to the MLE.
  mod_params_fitted <- load_best_params(search_profile=search_profile, archived_fit=FALSE, model_name=model_name)
  source("scripts/model_simulation/helper/put_params_into_model.R")
  model_fitted <- put_params_into_model(model_fitted, mod_params_fitted)

  source("scripts/model_fitting/helper/write_parameters_table.R")
  if(write_parameter_tables){
    write_parameter_tables(
      model_name,
      output_suffix,
      mod_params_fitted,
      country_under_investigation,
      n_of_units,
      cross_check_table
    )
  }

  return(model_fitted)

}