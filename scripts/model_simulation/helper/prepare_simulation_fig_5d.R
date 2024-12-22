prepare_simulation_fig_5d <- function(
  best_param_set,
  base_model_name,
  values_id_unit_intro,
  values_change_Reff,
  values_change_Immunity_setting,
  change_IDR,
  change_DSR,
  values_Scenarios,
  data_fitting,
  check_redo = FALSE
){
  if(check_redo){
    file.remove(list.files(dir_model_data, pattern = "df_all_values_M4_2\\d{3}.rds$", full.names = TRUE))
  }
  source("scripts/model_simulation/helper/get_params_for_travel_hubs.R")
  df_all_values_Strategy_Reff_Immunity <- bind_rows(
    lapply(seq_len(nrow(best_param_set)), function(set_n){
      # set_n = 1
      traveler_weight_i <- as.numeric(as.character(best_param_set$parval_1))[set_n]
      strategy_chosen = as.character(best_param_set$parval_2)[set_n]

      lapply(values_Scenarios, function(scenario_i){
        # scenario_i = values_Scenarios[1]
        if(scenario_i=="fewer_diagnostics"){
          change_IDR_i <- change_IDR$parval[change_IDR$strategy == strategy_chosen] %>% 
            gsub("%", "", .) %>% as.numeric()/100
          change_DSR_i <- 1
        } else if(scenario_i=="fewer_sequencing"){
          change_IDR_i <- 1
          change_DSR_i <- change_DSR$parval[change_DSR$strategy == strategy_chosen] %>% 
            gsub("%", "", .) %>% as.numeric()/100
        } else if(scenario_i=="half_resources"){
          change_IDR_i <- 0.5
          change_DSR_i <- 0.5
        } else{
          stop("Unknown strategy")
        }

        if(grepl("^P", strategy_chosen)){
          ranking_scheme_i = "per_capita"
          ST_type = "STp"
        } else if(grepl("^T", strategy_chosen)){
          ranking_scheme_i = "total"
          ST_type = "STt"
        } else{
          stop("Unknown strategy")
        }

        df_tmp <- get_params_for_travel_hubs(strategy_i = strategy_chosen, traveler_weight_i = traveler_weight_i, base_model_name = base_model_name, data_fitting = data_fitting, travel_volume_type = ranking_scheme_i, travel_data_period = "normal")

        lapply(values_change_Reff, function(change_Reff_i){
          lapply(values_change_Immunity_setting, function(change_Immunity_i){
            lapply(values_id_unit_intro, function(id_unit_intro){
              return(tibble(
                change_IDR = list(df_tmp$change_IDR*change_IDR_i),
                change_DSR = list(df_tmp$change_DSR*change_DSR_i),
                change_Reff = change_Reff_i,
                change_Immunity_setting = change_Immunity_i,
                traveler_weight = list(df_tmp$traveler_weight),
                id_unit_intro = list(id_unit_intro),
                scenario_name = scenario_i,
                strategy_chosen = strategy_chosen,
                scenario = paste0(ST_type, "_TW_ID_IDR_DSR_Reff_Immu_Scene_", strategy_chosen, "_", traveler_weight_i, "_", id_unit_intro, "_", change_IDR_i, "_", change_DSR_i, "_", change_Reff_i, "_", change_Immunity_i, "_", scenario_i)
              ))
            }) %>% bind_rows()
          }) %>% bind_rows()
        }) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  )

  ## baseline scenario have traveler weight = 1/100 for all regions
  df_all_values_baseline_Reff_Immunity <- bind_rows(
      lapply(values_change_Reff, function(change_Reff_i){
        lapply(values_change_Immunity_setting, function(change_Immunity_i){
          lapply(values_id_unit_intro, function(id_unit_intro){
            return(tibble(
              change_IDR = list(1),
              change_DSR = list(1),
              change_Reff = change_Reff_i,
              change_Immunity_setting = change_Immunity_i,
              traveler_weight = list(1/100),
              id_unit_intro = list(id_unit_intro),
              scenario_name = "full_resources",
              strategy_chosen = "baseline",
              scenario = paste0("STt_TW_ID_IDR_DSR_Reff_Immu_Scene_", "baseline_", 1/100, "_", id_unit_intro, "_", 1, "_", 1, "_", change_Reff_i, "_", change_Immunity_i, "_", "full_resources")
            ))
          }) %>% bind_rows()
        }) %>% bind_rows()
      }) %>% bind_rows()
  )

  df_all_values_todo <- rbind(df_all_values_Strategy_Reff_Immunity, df_all_values_baseline_Reff_Immunity)
  n_files <- 35
  ## cut the df_all_values into n_files files
  split_indices <- cut(seq(1, nrow(df_all_values_todo)), breaks = n_files, labels = FALSE)
  table(split_indices)
  for (i in 1:n_files) {
    df_all_values_i <- df_all_values_todo[split_indices == i, ]
    saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", 2000+i, ".rds"))
  }

  ## generate a bash command for submitting i slurm jobs
  bash_command <- c(
    "#!/bin/bash",
    "",
    paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 2001 ", 2000+n_files, " amd", " 32 ", model_name)
  )
  writeLines(bash_command, "scripts/model_simulation/Simulate_M4_4_HPC_command.sh")
  system("chmod +x scripts/model_simulation/Simulate_M4_3_HPC_command.sh")

  print("Done")
  return(df_all_values_todo)

}