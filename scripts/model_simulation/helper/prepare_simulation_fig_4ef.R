# prepare data and scripts for simulation on HPC
prepare_simulation_fig_4ef <- function(
  best_param_set,
  base_model_name,
  values_change_IDR = c(1/1000, seq(0.2, 0.8, 0.2), NA, seq(1.2, 1.8, 0.2)),
  values_change_DSR = c(1/1000, seq(0.2, 0.8, 0.2), NA, seq(1.2, 1.8, 0.2)),
  values_id_unit_intro = seq(0, 28, 1),
  data_fitting
){
  # simulate data for Fig. 4c
  df_all_values_TH_SP_ID_IDR <- bind_rows(
    lapply(values_change_IDR, function(change_IDR_i){
      df_tmp <- get_params_for_travel_hubs(top_n = as.numeric(as.character(best_param_set$parval_2)), base_model_name = base_model_name, data_fitting = data_fitting)
      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR*change_IDR_i),
          change_DSR = list(df_tmp$change_DSR),
          sequencing_propensity = list(as.numeric(as.character(best_param_set$parval_1))),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("TH_SP_ID_IDR_", best_param_set$parval_2, "_", best_param_set$parval_1, "_", id_unit_intro, 
          "_", change_IDR_i)
        ))
      }) %>% bind_rows()
    })
  )

  df_all_values_TH_SP_ID_DSR <- bind_rows(
    lapply(values_change_DSR, function(change_DSR_i){
      df_tmp <- get_params_for_travel_hubs(top_n = as.numeric(as.character(best_param_set$parval_2)), base_model_name = base_model_name, data_fitting = data_fitting)
      lapply(values_id_unit_intro, function(id_unit_intro){
        return(tibble(
          change_IDR = list(df_tmp$change_IDR),
          change_DSR = list(df_tmp$change_DSR*change_DSR_i),
          sequencing_propensity = list(as.numeric(as.character(best_param_set$parval_1))),
          id_unit_intro = list(id_unit_intro),
          scenario = paste0("TH_SP_ID_DSR_", best_param_set$parval_2, "_", best_param_set$parval_1, "_", id_unit_intro, 
          "_", change_DSR_i)
        ))
      }) %>% bind_rows()
    })
  )

  df_all_values_todo <- bind_rows(df_all_values_TH_SP_ID_IDR, df_all_values_TH_SP_ID_DSR)
  
  ## save the values to N files with 32 rows each for parallel computing
  n_files <- ceiling(nrow(df_all_values_todo)/32)
  files_rds_out <- lapply(seq_len(n_files), function(i){
    df_all_values_i <- df_all_values_todo[(i*32-31):min(i*32, nrow(df_all_values_todo)), ]
    saveRDS(df_all_values_i, paste0(dir_model_data, "/df_all_values_", model_name, "_", i+1000, ".rds"))
    return(paste0(dir_model_data, "/df_all_values_", model_name, "_", i+1000, ".rds"))
  })

  ## generate a bash command for submitting i slurm jobs
  bash_command <- c(
    "#!/bin/bash",
    "",
    paste0("scripts/model_simulation/helper/Simulate_Model_HPC_slurm.sh 1001 ", 1000+length(files_rds_out), " amd", " 64 ", model_name)
  )
  writeLines(bash_command, "scripts/model_simulation/Simulate_M4_3_HPC_command.sh")
  system("chmod +x scripts/model_simulation/Simulate_M4_3_HPC_command.sh")

  print("Done")
}