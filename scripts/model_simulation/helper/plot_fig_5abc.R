pacman::p_load(tidyverse, ggplot2, patchwork, geofacet, ggflags)
source("scripts/data_processing/helper/color_scheme.R")

# Fig 5a, to the left, is the geo_grid showing the boxplot of the identification lag for change_IDR and change_DSR. Change_IDR and change DSR should be dodged, in pink and blue, respectively. 
# Fig 5bc, to the right, showing the global summary of the identification lag for change_IDR and change_DSR. values are bootstrapped and shown as boxplot.
plot_fig_5abc <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  data_fitting,
  bootstrap_each_param=1024*10,
  n_cores
){
  df_all_values_todo <- df_all_values
  df_all_values_todo$id_intro <- mclapply(df_all_values_todo$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[7]
  }, mc.cores = n_cores) %>% unlist()
  df_all_values_todo <- df_all_values_todo %>% mutate(code_intro = factor(as.numeric(id_intro), 0:28, labels=sort(names(colors_spatial_units)))) %>% mutate(code_intro = as.character(code_intro))

  df_EDT_base_hline <- read_csv(paste0(dir_rst, "Fig4a_median_lag_case_type.csv")) %>% filter(parval_1=="1%", parval_2=="Max") # use the data from Fig 4a stratiied by code_intro and case_type.
  df_EDT_base_hline <- df_EDT_base_hline %>% filter(!is.na(code_intro)) 

  ## load simulation data and get EDT
  df_EDT <- lapply(seq_along(df_all_values_todo$scenario), function(i){
    model_name_i <- paste0("M4_", df_all_values_todo$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    if(nrow(df_sims) == 0){
      df_sims <- tibble(date=NA, unitname=NA, id_sim=as.character(rep(1:256, each=2)), variant = "Omicron BA.1", case_type = rep(c("Community", "Imported"), 256), value=NA, date_intro=NA, lag=Inf, scenario=NA, code_intro=NA)
    }
    df_sims$scenario <- df_all_values_todo$scenario[i]
    df_sims$code_intro <- df_all_values_todo$code_intro[i]
    return(df_sims)
  }) %>% bind_rows()

  # Fig 5a
  ## load geo_grid
  source("scripts/model_fitting/helper/geo_grid.R")
  layout_geo_grid_29 <- geo_grid_29()
  layout_geo_grid_29$code_intro <- layout_geo_grid_29$code
  layout_geo_grid_29$name_intro <- layout_geo_grid_29$name
  
  ## reformat data
  df_EDT$params <- mclapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[4]
  }, mc.cores = n_cores) %>% unlist()
  df_EDT$params <- factor(df_EDT$params, levels = c("IDR", "DSR"), labels = c("Changing diagnostics", "Changing sequencing"))
  df_EDT$parval <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    as.character(tmp[8])
  })
  df_EDT$parval <- as.numeric(df_EDT$parval)
  df_EDT$parval <- factor(df_EDT$parval, levels = sort(unique(df_EDT$parval)), labels = gsub(".0%", "%", scales::percent(sort(unique(df_EDT$parval)), accuracy = 0.1)))

  df_EDT$strategy <- mclapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    as.character(tmp[5])
  }, mc.cores = n_cores) %>% unlist()

  df_EDT$group <- factor(paste0(df_EDT$strategy, " (", df_EDT$params, ")"))
  df_EDT <- df_EDT %>% ungroup()

  names(colors_fig5abc) <- levels(df_EDT$group)

  set.seed(2024)
  df_EDT_plot <- df_EDT %>% group_by(code_intro, params, parval, group, id_sim) %>% slice_min(lag) %>% slice_sample(n=1)

  df_EDT_plot_text <- tibble(parval=c("0.1%", "100%"))
  df_EDT_plot_text <- cross_join(df_EDT_plot_text, layout_geo_grid_29 %>% select(code, name))
  df_EDT_plot_text$y <- 110
  df_EDT_plot_text$code_intro <- df_EDT_plot_text$code
  df_EDT_plot_text$code_flags <- tolower(df_EDT_plot_text$code)
  df_EDT_plot_text$parval <- ordered(df_EDT_plot_text$parval, levels=levels(df_EDT$parval))
  layout_geo_grid_29 %>% select(code, name) %>% arrange(code) %>% mutate(id=row_number()) %>% write_csv(paste0(dir_rst, "df_unitnames_29.csv"))

  df_EDT_base_hline_note <- df_EDT_base_hline %>% group_by(code_intro) %>% summarise(note = paste0("Baseline:\n", paste0(round(lag_median, 1), " (", substr(case_type, 1, 3), ")", collapse = ", ")))

  p_fig_5a <- ggplot(df_EDT_plot) + 
    geom_hline(data=df_EDT_base_hline, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
    geom_boxplot(aes(x=parval, y=lag, color=group, fill=group), alpha=0.9, outlier.shape = NA, size = 0.2) +
    geom_text(data=df_EDT_base_hline_note, aes(x="20%", y=60, label=note), hjust=0, vjust=0, size = 2.5, color="darkred") +
    scale_linetype_manual(name="", values=c("dashed", "solid"))+
    scale_color_manual(name="", values=colors_fig5abc)+
    scale_fill_manual(name="", values=colors_fig5abc)+
    geom_label(data = df_EDT_plot_text %>% filter(parval=="0.1%"), aes(x = parval, y = y, label = name), hjust = 0, size = 2.5, alpha=0.9) +
    geom_flag(data = df_EDT_plot_text %>% filter(!grepl("others", code)) %>% filter(parval!="0.1%"), aes(x = parval, y = y, country = code_flags), size = 5) +
    facet_geo(vars(code_intro), grid = layout_geo_grid_29, label="name") +
    ylab("Identification lag (days)")+
    xlab("Changing diagnostics/sequencing")+
    guides(color = guide_legend(nrow = 1)) +
    theme_minimal()+
    theme(
      panel.border = element_blank(),
      axis.text.x = element_text(angle=30),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "bottom"
      )+
    NULL
  ggsave(p_fig_5a, filename = paste0(dir_rst, "fig_5a.jpg"), width = 12, height = 10, dpi = 300)
  ggsave(p_fig_5a, filename = paste0(dir_rst, "fig_5a.pdf"), width = 12, height = 10, dpi = 300)
  
  # Fig 5bc, for both change_IDR and change_DSR:
  ## bootstrap sampling per the population scale
  ## For each parameter combination, we resample the simulations to match the population size of the region of introduction.

  df_population <- data_fitting$df_covar %>% group_by(code) %>% summarise(population=mean(population)) %>% ungroup() %>% arrange(population) # load population data for bootstrap sampling
  df_population$scale <- df_population$population/min(df_population$population)
  scaling_factor <- bootstrap_each_param/sum(df_population$scale)

  list_fig_fg <- lapply(1:2, function(i){
    this_fig_name <- c("fig_5b", "fig_5c")[i]
    df_EDT_i <- df_EDT %>% filter(params == c("Changing diagnostics", "Changing sequencing")[i])

    par_values <- df_EDT_i %>% pull(parval) %>% unique()
    codes_all <- sort(names(colors_spatial_units))

    registerDoParallel(n_cores)
    df_EDT_boot <- foreach(z = seq_along(par_values)) %dopar% {
      print(par_values[z])
      df_EDT_i_z <- df_EDT_i %>% filter(parval==par_values[z]) %>% select(-date_intro)
      id_all <- as.numeric(unique(df_EDT_i_z$id_sim))

      df_sims_out <- lapply(unique(df_EDT_i_z$group), function(group_i){
        # group_i = unique(df_EDT_i_z$group)[1]
        # for every code_intro, resample the simulations
        lapply(seq_along(codes_all), function(j){
          code_intro_j <- codes_all[j]
          df_sims_j <- df_EDT_i_z %>% filter(group==group_i) %>% filter(code_intro==code_intro_j)
          boot_n_j <- round(df_population$scale[df_population$code==code_intro_j]*scaling_factor)
          
          set.seed(2023)
          id_sims_j <- as.numeric(sample(id_all, boot_n_j, replace=TRUE)) 
          # keep the one with lowest lag between community and imported cases, and resample
          
          set.seed(2024)
          tmp <- df_sims_j %>% group_by(id_sim) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if there are multiple cases with the same lowest lag, randomly sample one.
          tmp[id_sims_j,]
        }) %>% bind_rows()
      })
      df_sims_out <- bind_rows(df_sims_out)

      df_sims_out$id_sim_new <- seq_len(nrow(df_sims_out))
      return(df_sims_out)
    }
    stopImplicitCluster()

    df_EDT_boot <- bind_rows(df_EDT_boot)
    df_EDT_boot <- df_EDT_boot %>% filter(!is.na(parval))

    ## add hline for the median, data from Fig 4c, stratified by community and imported cases.
    df_hline_median <- read_csv(paste0(dir_rst, "Fig4c_median_lag_case_type.csv")) %>% filter(parval_1=="1%", parval_2=="Max") %>% filter(!is.na(case_type)) # use the data from Fig 4c stratiied by code_intro and case_type.

    # Calculate the median and max values for each group
    df_medians <- df_EDT_boot %>%
      group_by(parval, group) %>%
      summarise(median_lag = median(lag), max_lag = quantile(lag, 0.75) + 1.5 * IQR(lag))

    p_fig_i <- ggplot(df_EDT_boot) + 
      geom_boxplot(aes(x=parval, y=lag, fill=group, color=group), alpha=0.9, outlier.shape = NA) +
      geom_hline(data=df_hline_median, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
      geom_label(data=tibble(parval=levels(df_EDT_boot$parval)[1], lag=Inf), aes(parval, lag), label="Global", size=5, alpha=0.9, hjust = 0, vjust = 1) +
      geom_text(data=df_medians, aes(x=parval, y=max_lag, label=round(median_lag, 1), color=group), vjust=-0.5, position=position_dodge(width=0.75)) +
      geom_text(data=df_hline_median, aes(x="40%", label=round(lag_median, 1), color=case_type), y=55, hjust=0, vjust=0, size = 2.5, color="darkred") +
      scale_color_manual(name="", values=colors_fig5abc)+
      scale_fill_manual(name="", values=colors_fig5abc)+
      scale_linetype_manual(name="", values=c("dashed", "solid"))+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # extend y-axis range
      ylab("Identification lag (days)")+
      xlab(c("Changing diagnostics", "Changing sequencing")[i])+
      theme_minimal()+
      theme(
        panel.border = element_blank(),
        axis.text.x = element_text(angle=30),
        legend.position = "none"
      )+
      NULL
    ggsave(p_fig_i, filename = paste0(dir_rst, this_fig_name, ".jpg"), width = 6, height = 5, dpi = 300)
    ggsave(p_fig_i, filename = paste0(dir_rst, this_fig_name, ".pdf"), width = 6, height = 5, dpi = 300)
    return(p_fig_i)
  })

  p_fig_5bc <- (list_fig_fg[[1]] + ggtitle("b")) + (list_fig_fg[[2]] + ggtitle("c")) + plot_layout(nrow=1)
  ggsave(p_fig_5bc, filename = paste0(dir_rst, "fig_5bc.jpg"), width = 12, height = 5, dpi = 300)
  ggsave(p_fig_5bc, filename = paste0(dir_rst, "fig_5bc.pdf"), width = 12, height = 5, dpi = 300)

  # for c and d, testing which parval can achieve shorter lag than the baseline
  df_baseline <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="1%", parval_2=="Max")

  df_fig_5b <- p_fig_5bc[[1]]$data
  values_all <- levels(df_fig_5b %>% pull(parval))
  p_values <- lapply(values_all, function(parval_i){
    lapply(unique(df_fig_5b$strategy), function(strategy_i){
      lags_IDR_i <- df_fig_5b %>% filter(parval==parval_i) %>% filter(strategy==strategy_i) %>% pull(lag) 
      t_test_rst <- t.test(lags_IDR_i, df_baseline$lag, "less")
      return(tibble(strategy = strategy_i, parval = parval_i, p_value=t_test_rst$p.value))
    })
  }) %>% bind_rows()
  p_values$parval <- factor(p_values$parval, levels=levels(df_EDT$parval))
  change_IDR_min_better <- p_values %>% group_by(strategy) %>% arrange(parval) %>% filter(p_value<=0.05) %>% slice(1)

  df_fig_5c <- p_fig_5bc[[2]]$data
  values_all <- levels(df_fig_5c %>% pull(parval))
  p_values <- lapply(values_all, function(parval_i){
    # parval_i = values_all[1]
    lapply(unique(df_fig_5c$strategy), function(strategy_i){
      # strategy_i = unique(df_fig_5c$strategy)[1]
      lags_DSR_i <- df_fig_5c %>% filter(parval==parval_i) %>% filter(strategy==strategy_i) %>% pull(lag)
      if(all(is.infinite(lags_DSR_i))){
        return(tibble(strategy = strategy_i, parval = parval_i, p_value=NA))
      } else{
        t_test_rst <- t.test(lags_DSR_i, df_baseline$lag, "less")
        return(tibble(strategy = strategy_i, parval = parval_i, p_value=t_test_rst$p.value))
      }
    })
  }) %>% bind_rows()
  p_values$parval <- factor(p_values$parval, levels=levels(df_EDT$parval))
  change_DSR_min_better <- p_values %>% group_by(strategy) %>% arrange(parval) %>% filter(p_value<=0.05) %>% slice(1)

  list_fig_5abc <- list(p_fig_5a=p_fig_5a, p_fig_5bc=p_fig_5bc, change_IDR_min_better=change_IDR_min_better, change_DSR_min_better=change_DSR_min_better)
  saveRDS(list_fig_5abc, paste0(dir_rst, "list_p_fig_5abc.rds"))
  return(list_fig_5abc)

}