source("scripts/data_processing/helper/color_scheme.R")

# Fig. 3a: Differences in earliest detection dates in four scenarios comparing to the base scenario. Four scenarios: (1) Increase detection capacity; (2) Increase overall sequencing capacity; (3) Adjust Traveler weight (SP) for imported cases; (4) Same resources but increase sequencing capacity in travel hubs. ([ggstatsplot](https://r-graph-gallery.com/web-violinplot-with-ggstatsplot.html)). For BA.1 and BA.2 separately.## upper panel: BA.1; lower panel: BA.2

plot_fig_3a <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores=64
){
  df_all_values_fig_3a <- df_all_values %>% filter(fig=="fig_3a")

  # load simulation data and get EDT
  df_EDT <- mclapply(seq_along(df_all_values_fig_3a$scenario), function(i){
    model_name_i <- paste0("M3_", df_all_values_fig_3a$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    df_sims$scenario <- df_all_values_fig_3a$scenario[i]
    return(df_sims)
  }, mc.cores = ncores) %>% bind_rows()

  df_EDT_base <- df_EDT %>% filter(scenario=="sequencing_1")
  df_EDT_base_hline <- df_EDT_base %>% group_by(variant, case_type) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(round(mean(lag), 4), 4)) %>% mutate(scenario="base")
  df_EDT_base_hline$lineage_type <- paste0(df_EDT_base_hline$variant, " (", tolower(df_EDT_base_hline$case_type), ")")
  print(df_EDT_base_hline)

  df_EDT$lineage_type <- paste0(df_EDT$variant, " (", tolower(df_EDT$case_type), ")")

  df_EDT$params <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[1]
  })
  df_EDT$params <- factor(df_EDT$params, levels = c("detection", "sequencing", "SP", "TH"), labels = c("Change IDR", "Change DSR", "Traveler weight", "Travel hubs"))
  df_EDT$parval <- sapply(df_EDT$scenario, function(scenario_i){
    tmp <- strsplit(scenario_i, "_")[[1]]
    tmp[2]
  })
  df_EDT$parval[df_EDT$parval=="NA" & df_EDT$params=="Change IDR"] <- "1"
  df_EDT$parval[df_EDT$parval=="NA" & df_EDT$params=="Change DSR"] <- "1"
  df_EDT$parval[df_EDT$parval=="NA" & df_EDT$params=="Travel hubs"] <- "Max"
  order_values <- c(replace_na(values_change_IDR,1), replace_na(values_change_DSR,1), values_sequencing_propensity, replace_na(values_top_n, 29)) %>% unique() %>% sort()
  order_values[order_values=="29"] <- "Max"
  df_EDT$parval <- factor(df_EDT$parval, levels = as.character(order_values))

  ylimit <- df_EDT %>% group_by(variant, case_type, params, parval) %>% summarise(lag_95th=quantile(lag, 0.95)) %>% summarise(ylim=max(lag_95th)) %>% pull(ylim) %>% max()

  # plot
  p_fig_3a <- ggplot(df_EDT) + 
    geom_hline(data=df_EDT_base_hline, aes(yintercept=lag_median, linetype=case_type), color="darkred") +
    # geom_violin(aes(x=parval, y=lag, fill=variant), alpha=0.5) +
    geom_boxplot(aes(x=parval, y=lag, color=lineage_type), alpha=0.5, outlier.shape = NA) +
    facet_grid(rows = vars(variant), cols = vars(params), scales="free_x", space="free_x") +
    scale_color_manual(name="", values=colors_lineage_type) +
    scale_linetype_manual(name="", values=c("dashed", "solid"))+
    scale_y_continuous(limits=c(0, ylimit)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1),
      # strip.text.y = element_text(angle = 0),
      axis.title.x = element_blank(),
      legend.position = "bottom"
      # legend.box="vertical"
    )+
    ylab("Identification lag (days)")+
    guides(color = guide_legend(nrow = 1)) +
    NULL
  ggsave(paste0(dir_rst, "fig_3a.jpg"), p_fig_3a, width=12, height=6)
  ggsave(paste0(dir_rst, "fig_3a.pdf"), p_fig_3a, width=12, height=6)
  writexl::write_xlsx(df_EDT, paste0(dir_rst, "fig_3a_source.xlsx"))
  writexl::write_xlsx(df_EDT_base_hline, paste0(dir_rst, "fig_3a_base_hline.xlsx"))
  return(p_fig_3a)
}

# Fig. 3b - 3g: 2-D density plots of the frequencies of achieving lowest detection lags, for 2 parameters. 3b: Traveler weight with detection capacity (BA.1/BA.2); 3c: Traveler weight with sequencing capacity (BA.1/BA.2); 3d: Traveler weight with travel hubs (BA.1/BA.2)

plot_fig_3bcdefg <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores
){
  name_figs_all <- c("fig_3b", "fig_3c", "fig_3d", "fig_3e", "fig_3f", "fig_3g")
  df_all_values_todo <- df_all_values %>% filter((fig %in% name_figs_all) | (scenario=="sequencing_1"))

  # load simulation data and get EDT
  df_EDT <- mclapply(seq_along(df_all_values_todo$scenario), function(i){
    model_name_i <- paste0("M3_", df_all_values_todo$scenario[i])
    df_sims <- readRDS(paste0(dir_model_data, "df_sims_", model_name_i, ".rds"))
    df_sims$scenario <- df_all_values_todo$scenario[i]
    df_sims$fig <- df_all_values_todo$fig[i]
    return(df_sims)
  }, mc.cores = ncores) %>% bind_rows()

  df_EDT_base <- df_EDT %>% filter(scenario=="sequencing_1")
  df_EDT_base_median <- df_EDT_base %>% group_by(variant, case_type) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(mean(lag), 4)) %>% mutate(scenario="base")
  df_EDT_base_median$lineage_type <- paste0(df_EDT_base_median$variant, " (", tolower(df_EDT_base_median$case_type), ")")

  list_p_fig_3bcdefg <- lapply(seq_along(name_figs_all), function(i){
    title_sub_fig <- c("b", "c", "d", "e", "f", "g")[i]
    name_fig_i <- name_figs_all[i]
    df_EDT_fig_i <- df_EDT %>% filter(fig==name_fig_i)
    df_EDT_fig_i <- df_EDT_fig_i %>% separate(scenario, into=c("params_1", "params_2", "parval_1", "parval_2"), sep="_") %>% mutate(
      params_1 = factor(params_1, levels = c("detection", "sequencing", "SP", "TH"), labels = c("Change IDR", "Change DSR", "Traveler weight", "Travel hubs")),
      parval_1 = as.character(parval_1),
      params_2 = factor(params_2, levels = c("detection", "sequencing", "SP", "TH"), labels = c("Change IDR", "Change DSR", "Traveler weight", "Travel hubs")),
      parval_2 = as.character(parval_2)
    )

    if(df_EDT_fig_i$params_2[1] == "Traveler weight"){
      # switch the order of params_1 and params_2 if params_2 is Traveler weight
      df_EDT_fig_i <- df_EDT_fig_i %>% mutate(
        tmp_params = params_1,
        tmp_parval = parval_1,
        params_1 = params_2,
        parval_1 = parval_2,
        params_2 = tmp_params,
        parval_2 = tmp_parval
      ) %>% select(-tmp_params, -tmp_parval)
    }

    if (df_EDT_fig_i$params_1[1] == "Travel hubs") {
      df_EDT_fig_i$parval_1 <- factor(df_EDT_fig_i$parval_1, levels = replace_na(as.character(values_top_n), "29"), labels = replace_na(as.character(values_top_n), "Max"))
    }
    if (df_EDT_fig_i$params_2[1] == "Travel hubs") {
      df_EDT_fig_i$parval_2 <- factor(df_EDT_fig_i$parval_2, levels = replace_na(as.character(values_top_n), "NA"), labels = replace_na(as.character(values_top_n), "Max"))
    }

    list_p_fig_i <- lapply(c("Omicron BA.1", "Omicron BA.2"), function(this_variant){
      df_EDT_fig_i_j <- df_EDT_fig_i %>% filter(variant==this_variant) # filter by variant
      set.seed(2024)
      df_EDT_fig_i_j <- df_EDT_fig_i_j %>% group_by(id_sim, variant, params_1, params_2, parval_1, parval_2) %>% arrange(lag, desc(value)) %>% slice_min(lag) %>% slice_max(value) %>% slice_sample() %>% ungroup() # for each id_sim should only sample either the imported or community cases, depending on the lowest lag, if two case types have with the same lowest lag, we will then choose the one with highest number of sequenced cases, if still ties, randomly sample one.
      df_EDT_fig_i_j_median <- df_EDT_fig_i_j %>% group_by(variant, params_1, params_2, parval_1, parval_2) %>% summarise(lag_median=round(median(lag), 4), lag_mean=round(mean(lag), 4)) %>% ungroup() # showing the median lag for each tile (heatmap background)
      write_csv(df_EDT_fig_i_j_median, paste0(dir_rst, name_fig_i, "_", gsub(" ", "_", this_variant), "_heatmap_source.csv"))
      df_EDT_fig_i_j_earliest <- df_EDT_fig_i_j %>% group_by(variant, params_1, params_2, parval_1, parval_2) %>% mutate(lag_median=round(median(lag), 4), value_median=round(median(value), 4)) %>% ungroup() %>% group_by(variant, id_sim) %>% slice_min(lag_median) %>% slice_max(value_median) %>% ungroup() # choose lowest lag_median and highest value_meidna, if there are ties, keep all.
      write_csv(df_EDT_fig_i_j_earliest, paste0(dir_rst, name_fig_i, "_", gsub(" ", "_", this_variant), "_label_source.csv"))
      df_EDT_fig_i_j_min <- df_EDT_fig_i_j_earliest %>% ungroup() %>% group_by(parval_1, parval_2) %>% summarise(N_Community=sum(case_type=="Community"), N_imported=sum(case_type=="Imported"), n_sim=nrow(.)) %>% ungroup() %>% mutate(percent_Community=round(N_Community/n_sim*100), percent_imported=round(N_imported/n_sim*100)) # showing the frequency (in N rounds of simulations) of achieving the lowest lag for each tile (label/text)

      range_color <- round(range(df_EDT_fig_i %>% group_by(id_sim, variant, params_1, params_2, parval_1, parval_2) %>% arrange(lag, desc(value)) %>% slice_min(lag) %>% slice_max(value) %>% ungroup() %>% group_by(variant, params_1, params_2, parval_1, parval_2) %>% summarise(lag_median=round(median(lag), 4)) %>% pull(lag_median)))

      midpoint_value <- df_EDT_base_median %>% filter(variant==this_variant) %>% pull(lag_median) %>% min()

      p_fig_i_j <- ggplot()+
        geom_tile(data=df_EDT_fig_i_j_median, aes(x=parval_1, y=parval_2, fill=lag_median))+
        geom_text(data=df_EDT_fig_i_j_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black", nudge_x = -0.2)+
        geom_label(data=df_EDT_fig_i_j_min %>% filter(percent_Community>=1 & percent_imported>=1) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", nudge_x = 0.2, alpha=0.9, label.padding=unit(0.1, "lines"))+
        geom_text(data=df_EDT_fig_i_j_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_Community>=1), aes(x=parval_1, y=parval_2, label=percent_Community), color="black")+
        geom_label(data=df_EDT_fig_i_j_min %>% filter(!(percent_Community>=1 & percent_imported>=1)) %>% filter(percent_imported>=1), aes(x=parval_1, y=parval_2, label=percent_imported), color="black", alpha=0.9, label.padding=unit(0.1, "lines"))+
        xlab(df_EDT_fig_i$params_1[1])+
        ylab(df_EDT_fig_i$params_2[1])+
        scale_fill_gradient2(name="Identification lag (days)", midpoint = midpoint_value, low = ("darkred"), mid = "white", high = ("darkblue"), limits = range_color)+
        scale_size_continuous(range=c(2, 3))+
        facet_grid(rows = vars(variant), scales="free_x", space="free_x")+
        theme_minimal()+
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.key.height = unit(0.3, 'cm'),
          # legend.key.width = unit(1, 'cm'),
          legend.position = "bottom"
          )+
        guides(size="none")+
        NULL

      return(p_fig_i_j)
      })

    p_fig_i <- (list_p_fig_i[[1]] + ggtitle(title_sub_fig)) + list_p_fig_i[[2]] + plot_layout(ncol=1, axes = "keep", guides = "keep") & theme(legend.position = "bottom")
    ggsave(paste0(dir_rst, name_fig_i, ".jpg"), p_fig_i, width=6, height=6)
    ggsave(paste0(dir_rst, name_fig_i, ".pdf"), p_fig_i, width=6, height=6)
    writexl::write_xlsx(df_EDT_fig_i, paste0(dir_rst, name_fig_i, "_source.xlsx"))
    return(p_fig_i)
  })

  p_fig_3bcdefg <- wrap_plots(list_p_fig_3bcdefg[[1]], list_p_fig_3bcdefg[[2]], list_p_fig_3bcdefg[[3]], list_p_fig_3bcdefg[[4]], list_p_fig_3bcdefg[[5]], list_p_fig_3bcdefg[[6]], nrow=2, byrow=TRUE)
  ggsave(paste0(dir_rst, "fig_3bcdefg.jpg"), p_fig_3bcdefg, width=12, height=16)
  ggsave(paste0(dir_rst, "fig_3bcdefg.pdf"), p_fig_3bcdefg, width=12, height=16)
  return(p_fig_3bcdefg)
}

plot_fig_3 <- function(
  dir_rst,
  dir_model_data,
  df_all_values,
  df_vline,
  ncores
){
  p_fig_3a <- plot_fig_3a(dir_rst, dir_model_data, df_all_values, df_vline, ncores)+ggtitle("a")
  p_fig_3bcdefg <- plot_fig_3bcdefg(dir_rst, dir_model_data, df_all_values, df_vline, ncores)

  p_fig_3 <- p_fig_3a + p_fig_3bcdefg + plot_layout(ncol=1, heights = c(0.8, 2))
  ggsave(paste0(dir_rst, "fig_3.jpg"), p_fig_3, width=15, height=20)
  ggsave(paste0(dir_rst, "fig_3.pdf"), p_fig_3, width=15, height=20)

}
