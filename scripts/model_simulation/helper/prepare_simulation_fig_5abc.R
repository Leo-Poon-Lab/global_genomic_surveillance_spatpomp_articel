# prepare data and scripts for simulation on HPC
prepare_simulation_fig_5abc <- function(
  best_param_set,
  base_model_name,
  values_change_IDR,
  values_change_DSR,
  values_id_unit_intro = seq(0, 28, 1),
  data_fitting,
  check_redo
){
  if(check_redo){
    file.remove(list.files(dir_model_data, pattern = "df_all_values_M4_1\\d{3}.rds$", full.names = TRUE))
  }
  source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
  # simulate data for Fig. 5a
  df_all_values_TH_SP_ID_IDR <- bind_rows(
    lapply(values_change_IDR, function(change_IDR_i){
      lapply(seq_len(nrow(best_param_set)), function(set_n){
        traveler_weight_i <- as.numeric(as.character(best_param_set$parval_1))[set_n]
        strategy_i = as.character(best_param_set$parval_2)[set_n]

        if(grepl("^P", strategy_i)){
          ranking_scheme_i = "per_capita"
          ST_type = "STp"
        } else if(grepl("^T", strategy_i)){
          ranking_scheme_i = "total"
          ST_type = "STt"
        } else{
          stop("Unknown strategy")
        }

        df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = ranking_scheme_i, travel_data_period = "normal")

        lapply(values_id_unit_intro, function(id_unit_intro){
          return(tibble(
            change_IDR = list(df_tmp$change_IDR*change_IDR_i),
            change_DSR = list(df_tmp$change_DSR),
            traveler_weight = list(df_tmp$traveler_weight),
            id_unit_intro = list(id_unit_intro),
            scenario = paste0(ST_type, "_TW_ID_IDR_", strategy_i, "_", traveler_weight_i, "_", id_unit_intro, "_", change_IDR_i)
          ))
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  )

  df_all_values_TH_SP_ID_DSR <- bind_rows(
    lapply(values_change_DSR, function(change_DSR_i){
      lapply(seq_len(nrow(best_param_set)), function(set_n){
        traveler_weight_i <- as.numeric(as.character(best_param_set$parval_1))[set_n]
        strategy_i = as.character(best_param_set$parval_2)[set_n]
        
        if(grepl("^P", strategy_i)){
          ranking_scheme_i = "per_capita"
          ST_type = "STp"
        } else if(grepl("^T", strategy_i)){
          ranking_scheme_i = "total"
          ST_type = "STt"
        } else{
          stop("Unknown strategy")
        }

        df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_i, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = ranking_scheme_i, travel_data_period = "normal")

        lapply(values_id_unit_intro, function(id_unit_intro){
          return(tibble(
            change_IDR = list(df_tmp$change_IDR),
            change_DSR = list(df_tmp$change_DSR*change_DSR_i),
            traveler_weight = list(df_tmp$traveler_weight),
            id_unit_intro = list(id_unit_intro),
            scenario = paste0(ST_type, "_TW_ID_DSR_", strategy_i, "_", traveler_weight_i, "_", id_unit_intro, "_", change_DSR_i)
          ))
        
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  )

  df_all_values_todo <- bind_rows(df_all_values_TH_SP_ID_IDR, df_all_values_TH_SP_ID_DSR)

  n_files <- 35
  ## cut the df_all_values into n_files files
  split_indices <- cut(seq(1, nrow(df_all_values_todo)), breaks = n_files, labels = FALSE)
  for (i in 1:n_files) {
    df_all_values_i <- df_all_values_todo[split_indices == i, ]
    saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", 1000+i, ".rds"))
  }

  ## generate a bash command for submitting i slurm jobs
  bash_command <- c(
    "#!/bin/bash",
    "",
    paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1001 ", 1000+n_files, " amd", " 32 ", model_name)
  )
  writeLines(bash_command, "scripts/model_simulation/Simulate_M4_3_HPC_command.sh")
  system("chmod +x scripts/model_simulation/Simulate_M4_3_HPC_command.sh")

  print("Done")
  return(df_all_values_todo)
}