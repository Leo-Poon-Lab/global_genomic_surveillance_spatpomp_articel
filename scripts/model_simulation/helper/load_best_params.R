library(naturalsort)

load_best_params <- function(search_profile, archived_fit=FALSE, model_name="Omicron20"){
  if(archived_fit){
    file_dir <- "results/model_data/archived/"
  } else{
    file_dir <- "results/model_data/"
  }
  if(search_profile=="again"){
    file_again_fit <- list.files(file_dir, paste0("fit_", model_name, "_.+again_.+rds$"), full.names = T) %>% naturalsort() %>% tail(1)
    output_suffix <- paste0("_", file_again_fit %>% str_extract("again_\\d+"))
    print(paste0(
      "The pre-fitted model is: ", 
      model_name, "_",
      file_again_fit
      ))
    model_fit <- readRDS(file_again_fit)
    names(model_fit) <- "search2"
  } else {
    model_fit <- readRDS(here::here(paste0(file_dir, "fit_", model_name, "_", search_profile, ".rds")))
    output_suffix <- paste0("_", search_profile)
    print(paste0(
      "The pre-fitted model is: ", 
      model_name, "_",
      search_profile
      ))
  }
  best_m <- which.max(model_fit$search2$logLiks$logLik)
  mod_params_fitted <- model_fit$search2$params[best_m, ]
  mod_params_fitted
}

get_output_suffix <- function(search_profile, archived_fit=FALSE, model_name="Omicron20"){
  if(archived_fit){
    file_dir <- "results/model_data/archived/"
  } else{
    file_dir <- "results/model_data/"
  }
  if(search_profile=="again"){
    file_again_fit <- list.files(file_dir, "fit_", model_name, "_.+again_.+rds$", full.names = T) %>% naturalsort() %>% tail(1)
    output_suffix <- paste0("_", file_again_fit %>% str_extract("again_\\d+"))
  } else {
    output_suffix <- paste0("_", search_profile)
  }
  output_suffix
}