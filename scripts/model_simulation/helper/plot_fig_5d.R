library(ggplot2)
library(parallel)
source("scripts/data_processing/helper/color_scheme.R")

plot_fig_5d <- function(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values = df_all_values,
  data_fitting = data_fitting,
  bootstrap_each_param=1024*10,
  n_cores = n_cores
){
  # In Figure 5d, we want to plot one single figure, using facet_grid, where rows are change_Immunity_setting and columns are change_Reff. In each of the grid, we plot boxplots, where x-axis is the scenarios (fewer_dignostic, fewer_sequencing, half_resources) with different colors (two selected stategies). Similar baseline horizontal lines used in figure 4 ef should be reused.

  # the values for the boxplots should be pre-processed, they are the global bootstrapped EDT values for each of the scenarios, for each of the change_Reff and change_Immunity_setting. The boxplots should be colored by the strategy used.

  # there are 3 change_Reff values and 3 change_Immunity_setting values, so we should have 9 girds in total. each grid have 2 strategies and 3 scenarios, so 6 boxplots in each grid. 

  values_strategies <- unique(df_all_values$strategy_chosen) %>% .[.!="baseline"] # remove the baseline 
  values_scenarios <- unique(df_all_values$scenario_name) %>% .[.!="full_resources"] # remove the baseline
  values_change_Reff <- unique(df_all_values$change_Reff)
  values_change_Immunity_setting <- unique(df_all_values$change_Immunity_setting)
  codes_id_all <- unlist(unique(df_all_values$id_unit_intro))
  codes_all <- unlist(unique(data_fitting$country_under_investigation))

  df_EDT_boot <- lapply(values_strategies, function(this_strategy){
    # this_strategy <- values_strategies[3]
    lapply(values_scenarios, function(this_scenario){
      # this_scenario <- values_scenarios[1]
      lapply(values_change_Reff, function(this_change_Reff){
        # this_change_Reff <- values_change_Reff[1]
        mclapply(values_change_Immunity_setting, function(this_change_Immunity_setting){
          # this_change_Immunity_setting <- values_change_Immunity_setting[2]
          print(paste(this_strategy, this_scenario, this_change_Reff, this_change_Immunity_setting))

          ## read the EDT simulations
          df_all_values_subset <- df_all_values %>%
            filter(strategy_chosen == this_strategy, scenario_name == this_scenario, change_Reff == this_change_Reff, change_Immunity_setting == this_change_Immunity_setting)
          df_EDT <- lapply(seq_len(nrow(df_all_values_subset)), function(i){
            # i=1
            this_intro_id <- unlist(df_all_values_subset$id_unit_intro[i])
            this_name <- df_all_values_subset$scenario[i]
            df_sims <- readRDS(paste0(dir_model_data, "df_sims_M4_", this_name, ".rds"))
            df_sims$model_name <- this_name
            df_sims$code_intro <- this_intro_id
            return(df_sims)
          }) %>% bind_rows()

          ## bootstrap sampling per the population scale
          ## For each parameter combination, we resample the simulations to match the population size of the region of introduction.

          df_population <- data_fitting$df_covar %>% group_by(code) %>% summarise(population=mean(population)) %>% ungroup() %>% arrange(population) # load population data for bootstrap sampling
          df_population$scale <- df_population$population/min(df_population$population)
          scaling_factor <- bootstrap_each_param/sum(df_population$scale)

          id_all <- as.numeric(unique(df_EDT$id_sim))

          df_sims_out <- lapply(codes_id_all, function(this_code_id){
            # this_code_id=codes_id_all[2]
            this_code <- codes_all[this_code_id+1]
            df_sims_j <- df_EDT %>% filter(code_intro==this_code_id)
            boot_n_j <- round(df_population$scale[df_population$code==this_code]*scaling_factor)

            set.seed(2023)
            id_sims_j <- as.numeric(sample(id_all, boot_n_j, replace=TRUE)) 
            # keep the one with lowest lag between community and imported cases, and resample

            set.seed(2024)
            tmp <- df_sims_j %>% group_by(id_sim) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if there are multiple cases with the same lowest lag, randomly sample one.
            tmp_sampled <- tmp[id_sims_j,]
            tmp_sampled$change_Immunity_setting <- this_change_Immunity_setting
            tmp_sampled$change_Reff <- this_change_Reff
            tmp_sampled$strategy_chosen <- this_strategy
            tmp_sampled$scenario_name <- this_scenario
            tmp_sampled$lag[is.na(tmp_sampled$lag)] <- Inf
            tmp_sampled
          }) %>% bind_rows()

          df_sims_out$id_sim_new <- seq_len(nrow(df_sims_out))
          return(df_sims_out)

        }, mc.cores = 16) %>% bind_rows()
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()

  df_EDT_boot_baseline <- lapply(values_change_Reff, function(this_change_Reff){
    # this_change_Reff <- values_change_Reff[1]
    mclapply(values_change_Immunity_setting, function(this_change_Immunity_setting){
      # this_change_Immunity_setting <- values_change_Immunity_setting[2]
      this_strategy <- "baseline"
      this_scenario <- "full_resources"
      print(paste(this_strategy, this_scenario, this_change_Reff, this_change_Immunity_setting))

      ## read the EDT simulations
      df_all_values_subset <- df_all_values %>%
        filter(strategy_chosen == this_strategy, scenario_name == this_scenario, change_Reff == this_change_Reff, change_Immunity_setting == this_change_Immunity_setting)
      df_EDT <- lapply(seq_len(nrow(df_all_values_subset)), function(i){
        # i=1
        this_intro_id <- unlist(df_all_values_subset$id_unit_intro[i])
        this_name <- df_all_values_subset$scenario[i]
        df_sims <- readRDS(paste0(dir_model_data, "df_sims_M4_", this_name, ".rds"))
        df_sims$model_name <- this_name
        df_sims$code_intro <- this_intro_id
        return(df_sims)
      }) %>% bind_rows()

      ## bootstrap sampling per the population scale
      ## For each parameter combination, we resample the simulations to match the population size of the region of introduction.

      df_population <- data_fitting$df_covar %>% group_by(code) %>% summarise(population=mean(population)) %>% ungroup() %>% arrange(population) # load population data for bootstrap sampling
      df_population$scale <- df_population$population/min(df_population$population)
      scaling_factor <- bootstrap_each_param/sum(df_population$scale)

      id_all <- as.numeric(unique(df_EDT$id_sim))

      df_sims_out <- lapply(codes_id_all, function(this_code_id){
        # this_code_id=codes_id_all[2]
        this_code <- codes_all[this_code_id+1]
        df_sims_j <- df_EDT %>% filter(code_intro==this_code_id)
        boot_n_j <- round(df_population$scale[df_population$code==this_code]*scaling_factor)

        set.seed(2023)
        id_sims_j <- as.numeric(sample(id_all, boot_n_j, replace=TRUE)) 
        # keep the one with lowest lag between community and imported cases, and resample

        set.seed(2024)
        tmp <- df_sims_j %>% group_by(id_sim) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if there are multiple cases with the same lowest lag, randomly sample one.
        tmp_sampled <- tmp[id_sims_j,]
        tmp_sampled$change_Immunity_setting <- this_change_Immunity_setting
        tmp_sampled$change_Reff <- this_change_Reff
        tmp_sampled$strategy_chosen <- this_strategy
        tmp_sampled$scenario_name <- this_scenario
        tmp_sampled$lag[is.na(tmp_sampled$lag)] <- Inf
        tmp_sampled
      }) %>% bind_rows()

      df_sims_out$id_sim_new <- seq_len(nrow(df_sims_out))
      return(df_sims_out)

    }, mc.cores = 16) %>% bind_rows()
  }) %>% bind_rows()
  
  ## add hline data from the median values 
  df_hline_median <- df_EDT_boot_baseline %>% group_by(change_Reff, change_Immunity_setting, case_type) %>% summarise(lag_median = median(lag), y_position = max(lag)*1.5)
  df_hline_median$Immunity_setting <- factor(df_hline_median$change_Immunity_setting, levels=c("More", "Moderate", "Few"), c("Low VE", "Moderate VE", "High VE")) # corresponding to 0.32, 0.595, and 0.87 of the immuned population are susceptible
  df_hline_median$Reff_type <- paste0("Relative Reff = ", df_hline_median$change_Reff)

  df_EDT_boot$Reff_type <- paste0("Relative Reff = ", df_EDT_boot$change_Reff)
  df_EDT_boot$Immunity_setting <- factor(df_EDT_boot$change_Immunity_setting, levels=c("More", "Moderate", "Few"), c("Low VE", "Moderate VE", "High VE")) # corresponding to 0.32, 0.595, and 0.87 of the immuned population are susceptible
  df_EDT_boot$Scenario <- factor(df_EDT_boot$scenario_name, levels=c("fewer_diagnostics", "fewer_sequencing", "half_resources"), c("Fewer diagnostics", "Fewer sequencing", "Half resources"))

  df_EDT_boot$color_name <- paste0(df_EDT_boot$strategy_chosen, " (", df_EDT_boot$Scenario, ")")
  colors_fig5d <- c(colors_fig5abc[1:2], "#4a5a8a", colors_fig5abc[3:4], "#4a8a76")
  names(colors_fig5d) <- sort(unique(df_EDT_boot$color_name))
  df_EDT_boot$color_name <- factor(df_EDT_boot$color_name, levels=names(colors_fig5d))

  # Calculate the median and max values for each group
  df_medians <- df_EDT_boot %>%
    group_by(color_name, Scenario, Reff_type, Immunity_setting) %>%
    summarise(median_lag = median(lag), max_lag = max(lag))

  p_fig_5d <- ggplot(df_EDT_boot, aes(x=Scenario, y=lag)) + 
    geom_violin(aes(color=color_name), width=1.5, alpha=0.5, position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(color=color_name), width=0.1, alpha=0.4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_hline(data=df_hline_median, aes(yintercept=lag_median, linetype=case_type), color="darkred", show.legend = FALSE) +
    geom_text(aes(y=max_lag, label=round(median_lag, 1), color=color_name), data=df_medians, position=position_dodge(width=0.8), vjust=-0.5) +
    geom_text(data=df_hline_median, aes(y=y_position, label=paste0("Baseline: ", round(lag_median, 1))), x=2, color="darkred") +
    scale_color_manual(name="", values=colors_fig5d)+
    scale_fill_manual(name="", values=colors_fig5d)+
    scale_linetype_manual(name="", values=c("dashed", "solid"))+
    facet_grid(Reff_type~Immunity_setting, scales="free")+
    ylab("Identification lag (days)")+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # extend y-axis range
    theme_minimal()+
    theme(
      panel.border = element_blank(),
      axis.text.x = element_text(angle=10),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "left"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))+
    NULL

  ggsave(p_fig_5d, filename = paste0(dir_rst, "Fig_5d.jpg"), width = 12, height = 8, dpi = 300)
  ggsave(p_fig_5d, filename = paste0(dir_rst, "Fig_5d.pdf"), width = 12, height = 8)

  # Save the plot as a list
  list_fig_5d <- list(p_fig_5d=p_fig_5d)
  saveRDS(list_fig_5d, paste0(dir_rst, "list_p_fig_5d.rds"))
  return(list_fig_5d)
}