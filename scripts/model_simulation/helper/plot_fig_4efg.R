source("scripts/data_processing/helper/color_scheme.R")
pacman::p_load(tidyverse, ggplot2, patchwork, geofacet, ggflags)

# Fig 4e, to the left, is the geo_grid showing the boxplot of the identification lag for change_IDR and change_DSR. Change_IDR and change DSR should be dodged, in pink and blue, respectively. 
# Fig 4fg, to the right, showing the global summary of the identification lag for change_IDR and change_DSR. values are bootstrapped and shown as boxplot.
plot_fig_4efg <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  data_fitting,
  bootstrap_each_param=1024*50,
  n_cores
){
  df_all_values_todo <- df_all_values
  df_all_values_todo$id_intro <- mclapply(df_all_values_todo$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[7]
  }, mc.cores = n_cores) %>% unlist()
  df_all_values_todo <- df_all_values_todo %>% mutate(code_intro = factor(as.numeric(id_intro), 0:28, labels=sort(names(colors_spatial_units)))) %>% mutate(code_intro = as.character(code_intro))

  df_EDT_base_hline <- read_csv(paste0(dir_rst, "Fig4a_median_lag_case_type.csv")) %>% filter(parval_1==0.01, parval_2=="Max") # use the data from Fig 4a stratiied by code_intro and case_type.
  df_EDT_base_hline <- df_EDT_base_hline %>% filter(!is.na(code_intro)) 

  ## load simulation data and get EDT
  df_EDT <- mclapply(seq_along(df_all_values_todo$scenario), function(i){
    model_name_i <- paste0("M4_", df_all_values_todo$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    df_sims$scenario <- df_all_values_todo$scenario[i]
    df_sims$code_intro <- df_all_values_todo$code_intro[i]
    return(df_sims)
  }, mc.cores = n_cores) %>% bind_rows()

  # Fig 4e
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
  df_EDT$params <- factor(df_EDT$params, levels = c("IDR", "DSR"), labels = c("Change IDR", "Change DSR"))
  df_EDT$parval <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    as.character(tmp[8])
  })
  df_EDT$parval[df_EDT$parval == "NA"] <- "1"
  df_EDT$parval <- naturalsort::naturalfactor(df_EDT$parval)

  df_EDT$group <- paste0(df_EDT$params, "_", df_EDT$case_type)
  df_EDT$group <- factor(df_EDT$group, levels = c("Change IDR_Community","Change IDR_Imported","Change DSR_Community","Change DSR_Imported"), labels = c("Change IDR (Community)","Change IDR (Imported)","Change DSR (Community)","Change DSR (Imported)"))
  df_EDT <- df_EDT %>% group_by(code_intro, case_type)

  df_EDT_plot_text <- tibble(parval=c("0.001", "1"))
  df_EDT_plot_text <- cross_join(df_EDT_plot_text, layout_geo_grid_29 %>% select(code, name))
  df_EDT_plot_text$y <- 70
  df_EDT_plot_text$code_intro <- df_EDT_plot_text$code
  df_EDT_plot_text$code_flags <- tolower(df_EDT_plot_text$code)
  df_EDT_plot_text$parval <- ordered(df_EDT_plot_text$parval, levels=levels(df_EDT$parval))
  layout_geo_grid_29 %>% select(code, name) %>% arrange(code) %>% mutate(id=row_number()) %>% write_csv(paste0(dir_rst, "df_unitnames_29.csv"))

  p_fig_4e <- ggplot(df_EDT) + 
    geom_hline(data=df_EDT_base_hline, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
    # geom_violin(aes(x=parval, y=lag, fill=variant), alpha=0.5) +
    geom_boxplot(aes(x=parval, y=lag, color=group, fill=group), alpha=0.5, outlier.shape = NA, size = 0.2) +
    scale_linetype_manual(name="", values=c("dashed", "solid"))+
    scale_color_manual(name="", values=colors_fig4e)+
    scale_fill_manual(name="", values=colors_fig4e)+
    geom_label(data = df_EDT_plot_text %>% filter(parval=="0.001"), aes(x = parval, y = y, label = name), hjust = 0, size = 2.5, alpha=0.8) +
    geom_flag(data = df_EDT_plot_text %>% filter(!grepl("others", code)) %>% filter(parval!="0.001"), aes(x = parval, y = y, country = code_flags), size = 5) +
    facet_geo(vars(code_intro), grid = layout_geo_grid_29, label="name") +
    ylab("Identification lag (days)")+
    xlab("Change IDR/DSR")+
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
  ggsave(p_fig_4e, filename = paste0(dir_rst, "fig_4e.jpg"), width = 12, height = 10, dpi = 300)
  ggsave(p_fig_4e, filename = paste0(dir_rst, "fig_4e.pdf"), width = 12, height = 10, dpi = 300)
  
  # Fig 4fg, for both change_IDR and change_DSR:
  ## bootstrap sampling per the population scale
  ## For each parameter combination, we resample the simulations to match the population size of the region of introduction.

  df_population <- data_fitting$df_covar %>% group_by(code) %>% summarise(population=mean(population)) %>% ungroup() %>% arrange(population) # load population data for bootstrap sampling
  df_population$scale <- df_population$population/min(df_population$population)
  scaling_factor <- bootstrap_each_param/sum(df_population$scale)

  list_fig_fg <- lapply(1:2, function(i){
    this_fig_name <- c("fig_4f", "fig_4g")[i]
    df_EDT_i <- df_EDT %>% filter(params == c("Change IDR", "Change DSR")[i])

    par_values <- df_EDT_i %>% pull(parval) %>% unique()
    codes_all <- sort(names(colors_spatial_units))

    registerDoParallel(n_cores)
    df_EDT_boot <- foreach(z = seq_along(par_values)) %dopar% {
      print(z)
      df_EDT_i_z <- df_EDT_i %>% filter(parval==par_values[z]) %>% select(-date_intro)
      id_all <- as.numeric(unique(df_EDT_i_z$id_sim))

      df_sims_out <- lapply(seq_along(codes_all), function(j){
        code_intro_j <- codes_all[j]
        df_sims_j <- df_EDT_i_z %>% filter(code_intro==code_intro_j)
        # resample both community and imported cases
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

    df_EDT_boot <- bind_rows(df_EDT_boot)
    df_EDT_boot <- df_EDT_boot %>% filter(!is.na(parval))

    ## add hline for the median, data from Fig 4c, stratified by community and imported cases.
    df_hline_median <- read_csv(paste0(dir_rst, "Fig4c_median_lag_case_type.csv")) %>% filter(parval_1==0.01, parval_2=="Max") %>% filter(!is.na(case_type)) # use the data from Fig 4c stratiied by code_intro and case_type.

    p_fig_i <- ggplot(df_EDT_boot) + 
      geom_boxplot(aes(x=parval, y=lag, fill=group, color=group), alpha=0.5, outlier.shape = NA) +
      geom_hline(data=df_hline_median, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
      geom_label(data=tibble(x=ordered(0.001, levels=levels(df_EDT_boot$parval)), y=70), aes(x, y), label="Global", size=5, alpha=0.8, hjust = 0) +
      scale_color_manual(name="", values=colors_fig4e)+
      scale_fill_manual(name="", values=colors_fig4e)+
      scale_linetype_manual(name="", values=c("dashed", "solid"))+
      ylab("Identification lag (days)")+
      xlab(c("Change IDR", "Change DSR")[i])+
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

  p_fig_4fg <- (list_fig_fg[[1]] + ggtitle("f")) + (list_fig_fg[[2]] + ggtitle("g")) + plot_layout(ncol=1)
  ggsave(p_fig_4fg, filename = paste0(dir_rst, "fig_4fg.jpg"), width = 5, height = 12, dpi = 300)
  ggsave(p_fig_4fg, filename = paste0(dir_rst, "fig_4fg.pdf"), width = 5, height = 12, dpi = 300)

  p_fig_4efg <- (p_fig_4e + ggtitle("e")) + p_fig_4fg + plot_layout(nrow=1, widths = c(3.5, 1)) 
  ggsave(p_fig_4efg, filename = paste0(dir_rst, "fig_4efg.jpg"), width = 12, height = 10, dpi = 300)
  ggsave(p_fig_4efg, filename = paste0(dir_rst, "fig_4efg.pdf"), width = 12, height = 10, dpi = 300)

  list_fig_4efg <- list(p_fig_4e=p_fig_4e, p_fig_4fg=p_fig_4fg, p_fig_4efg=p_fig_4efg)
  saveRDS(list_fig_4efg, paste0(dir_rst, "list_p_fig_4efg.rds"))
  return(list_fig_4efg)

}