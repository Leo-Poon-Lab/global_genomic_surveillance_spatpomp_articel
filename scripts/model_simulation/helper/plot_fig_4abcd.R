source("scripts/data_processing/helper/color_scheme.R")
pacman::p_load(tidyverse, ggplot2, patchwork, geofacet, ggflags, data.table)
setDTthreads(threads = n_cores)

# Fig. 4a: Global identification lag of the novel variant in different origins (geo_grid version of fig. 3b).
# Fig. 4b: global summary of Fig. 4a. Aggregating the region of introduction.
plot_fig_4abcd <- function(
  dir_rst,
  dir_model_data,
  df_all_values_todo,
  data_fitting,
  bootstrap_each_param_pair=1024*10,
  n_cores
){
  size_label_a=1.5
  size_label_bcd=2.5
  df_all_values_todo <- df_all_values_todo %>% separate(scenario, into=c("params_1", "params_2", "params_3", "parval_1", "parval_2", "parval_3"), sep="_", remove=FALSE) %>% mutate(code = factor(as.numeric(parval_3), 0:28, labels=sort(names(colors_spatial_units)))) 
  df_all_values_todo$code[is.na(df_all_values_todo$code)] <- "ZA"

  df_population <- data_fitting$df_covar %>% group_by(code) %>% summarise(population=mean(population)) %>% ungroup() %>% arrange(population) # load population data for bootstrap sampling
  df_population$scale <- df_population$population/min(df_population$population)
  scaling_factor <- bootstrap_each_param_pair/sum(df_population$scale)

  df_all_values_todo$parval_1[df_all_values_todo$params_1=="STp"] <- gsub("^T", "P", df_all_values_todo$parval_1[df_all_values_todo$params_1=="STp"])

  df_all_values_todo <- df_all_values_todo %>% mutate(
    params_1 = factor(params_1, levels = c("STp", "STt", "TW", "ID"), labels = c("Strategies", "Strategies", "Traveler weight (travel hubs)", "Region of introduction")),
    parval_1 = as.character(parval_1),
    params_2 = factor(params_2, levels = c("STp", "STt", "TW", "ID"), labels = c("Strategies", "Strategies", "Traveler weight (travel hubs)", "Region of introduction")),
    parval_2 = as.character(parval_2),
    params_3 = factor(params_3, levels = c("STp", "STt", "TW", "ID"), labels = c("Strategies", "Strategies", "Traveler weight (travel hubs)", "Region of introduction")),
    parval_3 = as.character(parval_3)
  )
  if(grepl("Traveler weight", df_all_values_todo$params_2[2])){
    # switch the order of params_1 and params_2 if params_2 is Traveler weight
    df_all_values_todo <- df_all_values_todo %>% mutate(
      tmp_params = params_1,
      tmp_parval = parval_1,
      params_1 = params_2,
      parval_1 = parval_2,
      params_2 = tmp_params,
      parval_2 = tmp_parval
    ) %>% select(-tmp_params, -tmp_parval)
  }
  df_all_values_todo$parval_1 <- as.numeric(df_all_values_todo$parval_1)
  df_all_values_todo$parval_1 <- factor(df_all_values_todo$parval_1, levels = sort(unique(df_all_values_todo$parval_1)), labels = gsub(".0%", "%", scales::percent(sort(unique(df_all_values_todo$parval_1)), accuracy = 0.1)))

  stopifnot(all(df_all_values_todo$params_3[!is.na(df_all_values_todo$params_3)]=="Region of introduction"))
  params_names <- c(df_all_values_todo$params_1[2], df_all_values_todo$params_2[2], df_all_values_todo$params_3[2])
  df_all_values_todo <- df_all_values_todo %>% select(-params_1, -params_2, -params_3)
  names(df_all_values_todo)[names(df_all_values_todo)=="parval_3"] <- "id_intro"
  df_all_values_todo <- df_all_values_todo %>% mutate(code_intro = factor(as.numeric(id_intro), 0:28, labels=sort(names(colors_spatial_units)))) %>% mutate(code_intro = as.character(code_intro))

  # fig. 4a
  ## load simulation data and get EDT (by variant, case_type and id_sim)
  source("scripts/model_simulation/helper/get_dates.R")
  df_all_values_todo <- df_all_values_todo %>% filter(scenario!="base")
  df_EDT <- mclapply(seq_along(df_all_values_todo$scenario), function(i){
    # i=8381
    model_name_i <- paste0("M4_", df_all_values_todo$scenario[i])
    df_sims <- tryCatch({
      readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    }, error = function(e) {
      print(i)
      stop(e)
    })
    df_sims$scenario <- df_all_values_todo$scenario[i]
    df_sims$code_intro <- df_all_values_todo$code_intro[i]
    df_sims$parval_1 <- df_all_values_todo$parval_1[i]
    df_sims$parval_2 <- df_all_values_todo$parval_2[i]
    return(df_sims)
  }, mc.cores = n_cores) %>% bind_rows()

  ## use SP=0.01, TH="Max" as the reference for the hline in Fig 4e
  df_EDT %>% filter(parval_1=="1%", parval_2=="Max") %>% group_by(code_intro, parval_1, parval_2, case_type) %>% summarise(lag_median=median(lag), lag_mean=mean(lag), value_median=round(median(value), 4), value_mean=round(mean(value), 4)) %>% ungroup() %>% write_csv(paste0(dir_rst, "Fig4a_median_lag_case_type.csv")) # this is useful for Fig 4e

  source("scripts/model_fitting/helper/geo_grid.R")
  layout_geo_grid_29 <- geo_grid_29()
  layout_geo_grid_29$code_intro <- layout_geo_grid_29$code
  layout_geo_grid_29$name_intro <- layout_geo_grid_29$name

  set.seed(2024)
  table(df_EDT$parval_2)
  df_EDT$parval_2 <- factor(df_EDT$parval_2, levels = c(values_strategies, gsub("^T", "P", values_strategies), "Max"))

  df_EDT_fig_4a <- df_EDT %>% filter(scenario!="base") %>% group_by(code_intro, id_sim, parval_1, parval_2) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if two case types have with the same lowest lag, randomly sample one.
  df_EDT_fig_4a_median <- df_EDT_fig_4a %>% group_by(code_intro, parval_1, parval_2) %>% summarise(lag_median=median(lag), lag_mean=mean(lag)) %>% ungroup() # showing the median lag for each tile (heatmap background)
  set.seed(2024)
  df_EDT_fig_4a_earliest <- df_EDT_fig_4a %>% group_by(code_intro, parval_1, parval_2) %>% mutate(lag_median=median(lag)) %>% ungroup() %>% group_by(code_intro, id_sim) %>% slice_min(lag_median) %>% ungroup() # choose lowest lag_median, if there are ties, keep all.
  write_csv(df_EDT_fig_4a_earliest, paste0(dir_rst, "/Fig4a_label_source.csv"))
  df_EDT_fig_4a_min <- df_EDT_fig_4a_earliest %>% ungroup() %>% group_by(code_intro) %>% mutate(n_sim_region=n()) %>% group_by(code_intro, parval_1, parval_2) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_sim_region=n()) %>% ungroup() %>% mutate(percent_Community=round(N_Community/n_sim_region*100), percent_imported=round(N_imported/n_sim_region*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag for each tile (label/text)

  ## Boot the data to match the population size of the region of introduction
  df_params_pairs <- df_all_values_todo %>% filter(scenario!="base") %>% select(parval_1, parval_2) %>% unique()
  codes_all <- sort(names(colors_spatial_units))
  registerDoParallel(n_cores)
  df_EDT_fig_4bcd <- foreach(i = seq_len(nrow(df_params_pairs))) %dopar% {
    # i=289
    print(i)

    df_sims <- df_EDT %>% filter(parval_1==df_params_pairs$parval_1[i], parval_2==df_params_pairs$parval_2[i]) %>% select(-date_intro)
    id_all <- as.numeric(unique(df_sims$id_sim))

    df_sims_out <- lapply(seq_along(codes_all), function(j){
      code_intro_j <- codes_all[j]
      df_sims_j <- df_sims %>% filter(code_intro==code_intro_j)
      set.seed(2024)
      df_sims_j <- df_sims_j %>% group_by(id_sim) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if there are multiple cases with the same lowest lag, randomly sample one.
      boot_n_j <- round(df_population$scale[df_population$code==code_intro_j]*scaling_factor)

      set.seed(2023)
      id_sims_j <- as.numeric(sample(id_all, boot_n_j, replace=TRUE))
      df_sims_j[id_sims_j,]
    })
    df_sims_out <- bind_rows(df_sims_out)

    df_sims_out$id_sim_new <- seq_len(nrow(df_sims_out))
    return(df_sims_out)
  }
  stopImplicitCluster()
  df_EDT_fig_4bcd <- bind_rows(df_EDT_fig_4bcd)
  df_EDT_fig_4bcd <- df_EDT_fig_4bcd %>% filter(scenario!="base")
  saveRDS(df_EDT_fig_4bcd, paste0(dir_rst, "df_EDT_fig_4bcd.rds"))

  ## Use SP=0.01, TH="Max" as the midpoint reference 
  midpoint_value <- df_EDT_fig_4bcd %>% filter(parval_1=="1%", parval_2=="Max") %>% pull(lag) %>% median()
  range_color <- round(range(df_EDT_fig_4a_median$lag_median))
  writeLines(paste0("midpoint_value: ", midpoint_value, ", range_color: ", range_color), con=paste0(dir_rst, "Fig4a_midpoint_range.txt"))

  df_EDT_fig_4a_median_plot <- df_EDT_fig_4a_median %>% filter(parval_2!="Max")
  df_EDT_fig_4a_min_plot <- df_EDT_fig_4a_min %>% filter(parval_2!="Max")

  p_fig_4a <- ggplot(data=df_EDT_fig_4a_median_plot)+
    geom_tile(aes(x=parval_1, y=parval_2, fill=lag_median))+
    geom_text(data=df_EDT_fig_4a_min_plot %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", nudge_x = -0.2, size=size_label_a)+
    geom_label(data=df_EDT_fig_4a_min_plot %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", nudge_x = 0.2, size=size_label_a, alpha=0.9, label.padding=unit(0.1, "lines"))+
    geom_text(data=df_EDT_fig_4a_min_plot %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", size=size_label_a)+
    geom_label(data=df_EDT_fig_4a_min_plot %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", size=size_label_a, alpha=0.9, label.padding=unit(0.1, "lines"))+
    xlab(params_names[1])+
    ylab(params_names[2])+
    scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color)+
    facet_geo(vars(code_intro), grid = layout_geo_grid_29, label="name") +
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle=30),
      legend.position = "bottom"
      )+
    guides(size="none")+
    NULL
  ggsave(p_fig_4a, filename=paste0(dir_rst, "fig_4a.jpg"), width=10, height=20, dpi=300)
  ggsave(p_fig_4a, filename=paste0(dir_rst, "fig_4a.pdf"), width=10, height=20)

  # Fig. 4bcd: global summary of Fig. 4a. Aggregating the region of introduction.
  ## bootstrap sampling per the population scale
  ## For each parameter combination, we resample the simulations to match the population size of the region of introduction.
  ## Fig. 4b: showing the 25th percentile
  ## Fig. 4c: showing the median
  ## Fig. 4d: showing the 75th percentile

  registerDoParallel(n_cores)
  df_EDT_boot_case_type <- foreach(i = seq_len(nrow(df_params_pairs))) %dopar% {
    print(i)
    df_sims <- df_EDT %>% filter(parval_1==df_params_pairs$parval_1[i], parval_2==df_params_pairs$parval_2[i]) %>% select(-date_intro)
    id_all <- as.numeric(unique(df_sims$id_sim))

    df_sims_out <- lapply(seq_along(codes_all), function(j){
      code_intro_j <- codes_all[j]
      df_sims_j <- df_sims %>% filter(code_intro==code_intro_j)
      set.seed(2024)
      df_sims_j <- df_sims_j %>% group_by(id_sim) %>% slice_min(lag) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if there are multiple cases with the same lowest lag, randomly sample one.
      boot_n_j <- round(df_population$scale[df_population$code==code_intro_j]*scaling_factor)

      set.seed(2023)
      id_sims_j <- as.numeric(sample(id_all, boot_n_j, replace=TRUE))
      bind_rows(df_sims_j %>% filter(case_type=="Community") %>% .[id_sims_j,], df_sims_j %>% filter(case_type=="Imported") %>% .[id_sims_j,])
    })
    df_sims_out <- bind_rows(df_sims_out)

    df_sims_out$id_sim_new <- seq_len(nrow(df_sims_out))
    return(df_sims_out)
  }
  stopImplicitCluster()
  df_EDT_boot_case_type <- bind_rows(df_EDT_boot_case_type)
  ## Use SP=0.01, TH="Max" as the reference for the hline in Fig 4fg
  df_EDT_boot_case_type %>% filter(parval_1=="1%", parval_2=="Max") %>% group_by(parval_1, parval_2, case_type) %>% summarise(lag_median=median(lag), lag_mean=mean(lag), value_median=round(median(value), 4), value_mean=round(mean(value), 4)) %>% ungroup() %>% write_csv(paste0(dir_rst, "Fig4c_median_lag_case_type.csv")) # this is useful for Fig 4fg
  rm(df_EDT_boot_case_type)

  ## reuse the boot data
  df_EDT_fig_4bcd_percentile <- df_EDT_fig_4bcd %>% filter(!parval_2=="Max")%>% group_by(parval_1, parval_2) %>% summarise(lag_median=median(lag), lag_mean=mean(lag), lag_25th=quantile(lag, 0.25), lag_75th=quantile(lag, 0.75)) %>% ungroup() # showing the median lag for each tile (heatmap background)

  ## Fig. 4b: showing the 25th percentile
  df_EDT_fig_4b_earliest <- df_EDT_fig_4bcd %>% filter(!parval_2=="Max")%>% group_by(parval_1, parval_2) %>% mutate(lag_25th=quantile(lag, 0.25), value_25th=quantile(value, 0.25)) %>% ungroup() %>% group_by(id_sim_new) %>% slice_min(lag_25th)
  saveRDS(df_EDT_fig_4b_earliest, paste0(dir_rst, "/Fig4b_label_source.rds"))
  df_EDT_fig_4b_min <- df_EDT_fig_4b_earliest %>% ungroup() %>% group_by(parval_1, parval_2) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_total=n()) %>% ungroup()
  df_EDT_fig_4b_min <- df_EDT_fig_4b_min %>% mutate(percent_Community=round(N_Community/n_total*100), percent_imported=round(N_imported/n_total*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag for each tile (label/text)
  df_EDT_fig_4bcd_percentile$facet_4b <- "Global (25th percentile)"
  df_EDT_fig_4bcd_percentile$facet_4c <- "Global (Median)"
  df_EDT_fig_4bcd_percentile$facet_4d <- "Global (75th percentile)"
  
  p_fig_4b <- ggplot(data=df_EDT_fig_4bcd_percentile)+
    geom_tile(aes(x=parval_1, y=parval_2, fill=lag_25th))+
    geom_text(data=df_EDT_fig_4b_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", nudge_x = -0.2, size=size_label_bcd)+
    geom_label(data=df_EDT_fig_4b_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", nudge_x = 0.2, size=size_label_bcd, alpha=0.9, label.padding=unit(0.1, "lines"))+
    geom_text(data=df_EDT_fig_4b_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", size=size_label_bcd)+
    geom_label(data=df_EDT_fig_4b_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", size=size_label_bcd, alpha=0.9, label.padding=unit(0.1, "lines"))+
    facet_wrap(~facet_4b, nrow=1)+
    xlab(params_names[1])+
    ylab(params_names[2])+
    scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color)+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle=30),
      legend.position = "none"
      )+
    guides(size="none")+
    NULL
  ggsave(p_fig_4b, filename=paste0(dir_rst, "fig_4b.jpg"), width=4, height=4, dpi=300)
  ggsave(p_fig_4b, filename=paste0(dir_rst, "fig_4b.pdf"), width=4, height=4)

  ## Fig. 4c: showing the median
  df_EDT_fig_4c_earliest <- df_EDT_fig_4bcd %>% filter(!parval_2=="Max")%>% group_by(parval_1, parval_2) %>% mutate(lag_50th=quantile(lag, 0.50), value_50th=quantile(value, 0.50)) %>% ungroup() %>% group_by(id_sim_new) %>% slice_min(lag_50th)
  saveRDS(df_EDT_fig_4c_earliest, paste0(dir_rst, "/Fig4c_label_source.rds"))
  df_EDT_fig_4c_min <- df_EDT_fig_4c_earliest %>% ungroup() %>% group_by(parval_1, parval_2) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_total=n()) %>% ungroup()

  df_EDT_fig_4c_min <- df_EDT_fig_4c_min %>% mutate(percent_Community=round(N_Community/n_total*100), percent_imported=round(N_imported/n_total*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag for each tile (label/text)
  p_fig_4c <- ggplot(data=df_EDT_fig_4bcd_percentile)+
    geom_tile(aes(x=parval_1, y=parval_2, fill=lag_median))+
    geom_text(data=df_EDT_fig_4c_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", nudge_x = -0.2, size=size_label_bcd)+ 
    geom_label(data=df_EDT_fig_4c_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", nudge_x = 0.2, alpha=0.9, label.padding=unit(0.1, "lines"), size=size_label_bcd)+
    geom_text(data=df_EDT_fig_4c_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", size=size_label_bcd)+
    geom_label(data=df_EDT_fig_4c_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", alpha=0.9, label.padding=unit(0.1, "lines"), size=size_label_bcd)+
    facet_wrap(~facet_4c, nrow=1)+
    xlab(params_names[1])+
    ylab(params_names[2])+
    scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color)+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle=30),
      legend.position = "none"
      )+
    guides(size="none")+
    NULL
  ggsave(p_fig_4c, filename=paste0(dir_rst, "fig_4c.jpg"), width=4, height=4, dpi=300)
  ggsave(p_fig_4c, filename=paste0(dir_rst, "fig_4c.pdf"), width=4, height=4)

  ## Fig. 4d: showing the 75th percentile
  df_EDT_fig_4d_earliest <- df_EDT_fig_4bcd %>% filter(!parval_2=="Max")%>% group_by(parval_1, parval_2) %>% mutate(lag_75th=quantile(lag, 0.75), value_75th=quantile(value, 0.75)) %>% ungroup() %>% group_by(id_sim_new) %>% slice_min(lag_75th)
  saveRDS(df_EDT_fig_4d_earliest, paste0(dir_rst, "/Fig4d_label_source.rds"))
  df_EDT_fig_4d_min <- df_EDT_fig_4d_earliest %>% ungroup() %>% group_by(parval_1, parval_2) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_total=n()) %>% ungroup()
  df_EDT_fig_4d_min <- df_EDT_fig_4d_min %>% mutate(percent_Community=round(N_Community/n_total*100), percent_imported=round(N_imported/n_total*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag for each tile (label/text)
  p_fig_4d <- ggplot(data=df_EDT_fig_4bcd_percentile)+
    geom_tile(aes(x=parval_1, y=parval_2, fill=lag_75th))+
    geom_text(data=df_EDT_fig_4d_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", nudge_x = -0.2, size=size_label_bcd)+
    geom_label(data=df_EDT_fig_4d_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", nudge_x = 0.2, alpha=0.9, label.padding=unit(0.1, "lines"), size=size_label_bcd)+
    geom_text(data=df_EDT_fig_4d_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", size=size_label_bcd)+
    geom_label(data=df_EDT_fig_4d_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", alpha=0.9, label.padding=unit(0.1, "lines"), size=size_label_bcd)+
    facet_wrap(~facet_4d, nrow=1)+
    xlab(params_names[1])+
    ylab(params_names[2])+
    scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color)+
    theme_minimal()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle=30),
      legend.position = "none"
      )+
    guides(size="none")+
    NULL
  ggsave(p_fig_4d, filename=paste0(dir_rst, "fig_4d.jpg"), width=4, height=4, dpi=300)
  ggsave(p_fig_4d, filename=paste0(dir_rst, "fig_4d.pdf"), width=4, height=4)

  p_fig_4bcd <- wrap_plots(p_fig_4b+ggtitle("b"), p_fig_4c+ggtitle("c"), p_fig_4d+ggtitle("d"), ncol=1)
  p_fig_4abcd <- (p_fig_4a+ggtitle("a")) + p_fig_4bcd + plot_layout(nrow=1, width = c(3.5,1))
  ggsave(p_fig_4abcd, filename=paste0(dir_rst, "fig_4abcd.jpg"), width=16, height=26, dpi=300)
  ggsave(p_fig_4abcd, filename=paste0(dir_rst, "fig_4abcd.pdf"), width=16, height=26)

  list_p_fig_4abcd <- list(
    p_fig_4a=p_fig_4a, 
    p_fig_4b=p_fig_4b,
    p_fig_4c=p_fig_4c,
    p_fig_4d=p_fig_4d,
    p_fig_4bcd=p_fig_4bcd,
    p_fig_4abcd=p_fig_4abcd,
    df_EDT_fig_4bcd_percentile=df_EDT_fig_4bcd_percentile)
  saveRDS(list_p_fig_4abcd, paste0(dir_rst, "list_p_fig_4abcd.rds"))

  return(list_p_fig_4abcd)
}
