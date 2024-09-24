put_params_into_model <- function(
    input_model,
    input_params
){
  mod_params <- coef(input_model)
  param_names <- names(mod_params)
  
  for(par_name in param_names){
    base_name <- gsub("\\d+$", "", par_name)
    if(base_name %in% names(input_params)){
      mod_params[par_name] <- input_params[[base_name]]
    } else if(par_name %in% names(input_params)){
      mod_params[par_name] <- input_params[[par_name]]
    }
  }
  
  coef(input_model) <- mod_params
  input_model
}
