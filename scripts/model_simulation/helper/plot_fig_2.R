source("scripts/data_processing/install_prerequisite.R")
pacman::p_load(dplyr, tidyr, reshape, tidyverse, fuzzyjoin, ggstream, colorspace, ggtext, cowplot, lubridate, ggrepel, patchwork, ggflags, ggalt)
source("scripts/data_processing/helper/color_scheme.R")

scientific_10 <- function(x) {
  scales::scientific_format()(x) %>% 
    str_remove_all("e\\+00") %>% # strip e+00 as we don't need it
    str_replace_all("e", " %*% 10^") %>%
    parse(text = .)
}

plot_fig_2a <- function(
  model_sims,
  dir_rst
){
  # aggregate the imported cases in simulation results
  df_sims <- model_sims %>%
    mutate(
      C_1_i_infected_new = rowSums(select(., starts_with("C_1_i_infected_origin_"))),
      C_1_i_detected_new = rowSums(select(., starts_with("C_1_i_detected_origin_"))),
      C_1_i_sequenced_new = rowSums(select(., starts_with("C_1_i_sequenced_origin_"))),
      C_2_i_infected_new = rowSums(select(., starts_with("C_2_i_infected_origin_"))),
      C_2_i_detected_new = rowSums(select(., starts_with("C_2_i_detected_origin_"))),
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_sequenced_origin_"))),
      C_3_i_infected_new = rowSums(select(., starts_with("C_3_i_infected_origin_"))),
      C_3_i_detected_new = rowSums(select(., starts_with("C_3_i_detected_origin_"))),
      C_3_i_sequenced_new = rowSums(select(., starts_with("C_3_i_sequenced_origin_"))),
      C_4_i_infected_new = rowSums(select(., starts_with("C_4_i_infected_origin_"))),
      C_4_i_detected_new = rowSums(select(., starts_with("C_4_i_detected_origin_"))),
      C_4_i_sequenced_new = rowSums(select(., starts_with("C_4_i_sequenced_origin_")))
    )

  prefixes <- paste0("C_", 1:4, "_", rep(c("i", "c"), each=4), "_")
  suffixes <- c("infected_new", "detected_new", "sequenced_new")
  sum_cols <- c(paste0(rep(prefixes, each=length(suffixes)), rep(suffixes, times=length(prefixes))))
  df_sims <- subset(df_sims, select = c("date_decimal", ".id", "unitname", sum_cols))

  # take the median of the simulation results for all N replicates, then add up different spatial units
  df_sims <- df_sims %>% group_by(date_decimal, unitname) %>% summarise_all(median) %>% ungroup()
  df_sims <- df_sims %>% select(-`.id`) %>% group_by(date_decimal) %>% summarise_at(vars(starts_with("C_")),sum) %>% ungroup()

  # change to long format
  df_sims_long <- df_sims %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value") %>% filter(value >= 1) %>% mutate(value=round(value))
  df_sims_long <- df_sims_long %>%
    separate(variable, into = c("cases", "variants", "cases_type", "cases_level", "new"), sep = "_") %>% 
    mutate(cases_type = ifelse(cases_type == "c", "community", "imported")) %>% 
    select(date_decimal, variants, cases_type, cases_level, value)

  df_sims_long$variants <- factor(df_sims_long$variants, levels=1:4, labels=c("Delta", "Omicron BA.1", "Omicron BA.2", "Others"))

  df_sims_long$date <- date_decimal(df_sims_long$date_decimal)
  df_sims_long$log_value <- log10(df_sims_long$value+1)

  # factorize different cases levels
  df_sims_long$stream_type <- paste(df_sims_long$variants," (",df_sims_long$cases_type,")",sep="")
  df_sims_long$cases_level <- factor(df_sims_long$cases_level,levels=c("infected", "detected", "sequenced"), labels = c("Infected", "Diagnostic", "Sequenced"))
  df_sims_long$stream_type <- factor(df_sims_long$stream_type, 
    levels=c("Delta (community)", "Delta (imported)", 
            "Omicron BA.1 (community)", "Omicron BA.1 (imported)", 
            "Omicron BA.2 (community)", "Omicron BA.2 (imported)", 
            "Others (community)", "Others (imported)")
    )

  # Add introduction dates as vline
  lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
  df_vline <- tibble(lag=lags_intro)
  df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(df_sims_long$date_decimal)
  df_vline$date <- date_decimal(df_vline$date_decimal)
  df_vline$v_line_label <- c("Omicron BA.1 emergence", "Omicron BA.2 emergence")
  df_vline$color <- c("#25567d", "#357933")
  df_vline$cases_level <- factor("Infected", levels=c("Infected", "Diagnostic", "Sequenced"))
  df_vline$v_line_label_y <- max(df_sims_long$value)

  # change alpha for Delta and others
  colors_lineage_type_fig2a <- colors_lineage_type
  colors_lineage_type_fig2a[grepl("Delta", names(colors_lineage_type_fig2a))] <- paste0(colors_lineage_type_fig2a[grepl("Delta", names(colors_lineage_type_fig2a))], "80")
  colors_lineage_type_fig2a[grepl("Others", names(colors_lineage_type_fig2a))] <- paste0(colors_lineage_type_fig2a[grepl("Others", names(colors_lineage_type_fig2a))], "80")

  saveRDS(df_sims_long, "results/figs/model_simulation/Omicron20/df_fig_2a_source.rds")

  p_fig_2a <- df_sims_long %>% filter(date_decimal>min(date_decimal)) %>% ggplot(
      aes(
        x=date, y=value, 
        color = stream_type, 
        fill = stream_type
      )
    ) +
    geom_area(
      stat = "identity",
      size=0
      )+
    facet_wrap( ## needs facet_grid for space argument
      vars(cases_level), 
      ncol=1,
      scales = "free",
      strip.position = "top"
    ) +
    scale_color_manual(
      expand = c(0, 0),
      values = colors_lineage_type_fig2a,
      guide = "none"
    ) +
    scale_fill_manual(
      values = colors_lineage_type_fig2a,
      name = NULL
    ) +
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    # Appearence type label on each panel
    geom_label_repel(
      data = df_vline %>% filter(color=="#25567d"),
      aes(x=date, y=v_line_label_y, label = v_line_label),
      color = "#25567d",
      inherit.aes = FALSE,
      size = 3
    ) +
    geom_label_repel(
      data = df_vline %>% filter(color=="#357933"),
      aes(x=date, y=v_line_label_y, label = v_line_label),
      color = "#357933",
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_y_continuous(
      name="Cases",
      label=scientific_10
    )+
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "", base_size = 12) +
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = .5,
        margin = margin(10, 0, 30, 0)
      ),
      plot.caption = element_text(
        size = 9,
        color = "black",
        hjust = .5,
        margin = margin(20, 0, 5, 0)
      ),
      # axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
      axis.title.x = element_blank(),
      # plot.background = element_rect(fill = "grey88", color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      # strip.text.y = element_blank(),
      legend.position = "none"
      # legend.box.margin = margin(t = 30), 
      # legend.text = element_text(size = 9, color = "grey40"),
      # legend.key.height = unit(.25, "lines"),
      # legend.key.width = unit(2.5, "lines"),
      # plot.margin = margin(rep(20, 4))
    ) +
    NULL

  ggsave(paste0(dir_rst, "/fig_2a.jpg"), plot=p_fig_2a, width=6*2, height=6, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2a.pdf"), plot=p_fig_2a, width=6*2, height=6, dpi=400)

  return(p_fig_2a)
}

# Fig 2b shows the distribution of the earliest date of detection of the first case in each spatial unit
plot_fig_2b <- function(
  model_sims,
  dir_rst
){
  # Add introduction dates as vline
  lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
  df_vline <- tibble(lag=lags_intro)
  df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(model_sims$date_decimal)
  df_vline$date <- date_decimal(df_vline$date_decimal)
  df_vline$v_line_label <- c("Omicron BA.1 emergence", "Omicron BA.2 emergence")
  df_vline$color <- c("#25567d", "#357933")
  df_vline$cases_level <- factor("Infected", levels=c("Infected", "Diagnostic", "Sequenced"))

  # aggregate the imported cases in simulation results
  df_sims <- model_sims %>%
    mutate(
      C_1_i_infected_new = rowSums(select(., starts_with("C_1_i_infected_origin_"))),
      C_1_i_detected_new = rowSums(select(., starts_with("C_1_i_detected_origin_"))),
      C_1_i_sequenced_new = rowSums(select(., starts_with("C_1_i_sequenced_origin_"))),
      C_2_i_infected_new = rowSums(select(., starts_with("C_2_i_infected_origin_"))),
      C_2_i_detected_new = rowSums(select(., starts_with("C_2_i_detected_origin_"))),
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_sequenced_origin_"))),
      C_3_i_infected_new = rowSums(select(., starts_with("C_3_i_infected_origin_"))),
      C_3_i_detected_new = rowSums(select(., starts_with("C_3_i_detected_origin_"))),
      C_3_i_sequenced_new = rowSums(select(., starts_with("C_3_i_sequenced_origin_"))),
      C_4_i_infected_new = rowSums(select(., starts_with("C_4_i_infected_origin_"))),
      C_4_i_detected_new = rowSums(select(., starts_with("C_4_i_detected_origin_"))),
      C_4_i_sequenced_new = rowSums(select(., starts_with("C_4_i_sequenced_origin_")))
    )

  prefixes <- paste0("C_", 2:3, "_", rep(c("i", "c"), each=2), "_")
  suffixes <- c("infected_new", "detected_new", "sequenced_new")
  sum_cols <- c(paste0(rep(prefixes, each=length(suffixes)), rep(suffixes, times=length(prefixes))))
  df_sims <- subset(df_sims, select = c("date_decimal", ".id", "unitname", sum_cols))
  
  # get the earliest date and value 
  detect_cols <- names(df_sims)[grepl("_new$", names(df_sims))]
  result <- lapply(detect_cols, function(col){
    df_sims %>%
      arrange(date_decimal) %>%
      group_by(`.id`, unitname) %>%
      filter(get(col) >= 1) %>%
      slice_head() %>%
      transmute(`.id`, unitname, earliest_date_decimal = date_decimal, num_cases=get(col), measured_metric = col) %>%
      ungroup()
  })
  result <- bind_rows(result)

  # add types and levels
  result <- result %>%
    separate(measured_metric, into = c("cases", "variants", "lineage_type", "cases_level", "new"), sep = "_") %>%
    mutate(lineage_type = ifelse(lineage_type == "c", "community", "imported"),
           variants = ifelse(variants == "2", "Omicron BA.1", "Omicron BA.2")) %>%
    select(.id, unitname, variants, lineage_type, cases_level, earliest_date_decimal, num_cases)

  # add date in date format
  result$earliest_date <- date_decimal(result$earliest_date_decimal)

  # # sum the simulation times
  # df_fig_2b <- result %>%
  #   group_by(earliest_date, unitname, lineage_type, variants, cases_level) %>%
  #   summarise(n_simulation = sum(num_cases>0, na.rm = T)) %>% ungroup()
  
  # if use boxplot, no need to summarize simulation times
  df_fig_2b <- result
  # factorize different cases levels
  df_fig_2b$color_type <- paste0(df_fig_2b$variants, " (", df_fig_2b$lineage_type, ")")
  df_fig_2b$color_type <- factor(df_fig_2b$color_type, levels = rev(c(
    "Delta (imported)", "Delta (community)", 
    "Omicron BA.1 (imported)", "Omicron BA.1 (community)", 
    "Omicron BA.2 (imported)", "Omicron BA.2 (community)",
    "Others (imported)", "Others (community)"
  )))
  df_fig_2b$lineage_type <- factor(df_fig_2b$lineage_type,levels = c("imported","community"), labels = c("Imported", "Community"))
  df_fig_2b$cases_level <- factor(df_fig_2b$cases_level, levels = c("sequenced","detected","infected"), labels = c("Sequenced", "Diagnostic", "Infected"))
  df_fig_2b$variants_level <- factor(df_fig_2b$variants,levels = c("Omicron BA.1","Omicron BA.2"), labels = c("Omicron BA.1", "Omicron BA.2"))

  # add formal unitname 
  df_fig_2b <- left_join(df_fig_2b, cross_check_table %>% transmute(unitname=code, loc_name), by="unitname")
  df_fig_2b$loc_name[is.na(df_fig_2b$loc_name)] <- df_fig_2b$unitname[is.na(df_fig_2b$loc_name)]
  df_fig_2b$loc_name[df_fig_2b$loc_name=="Hong Kong Special Administrative Region of China"] <- "Hong Kong SAR"

  # reorder the spatial units by earliest date
  loc_name_order <- df_fig_2b %>% filter(cases_level=="Infected", variants_level=="Omicron BA.1", lineage_type=="Imported") %>% group_by(loc_name) %>% mutate(min_date=min(earliest_date_decimal), median_date=median(earliest_date_decimal)) %>% arrange(median_date, min_date) %>% pull(loc_name) %>% unique()
  loc_name_order <- c("South Africa", as.character(loc_name_order)) %>% unique()
  df_fig_2b$loc_name <- factor(df_fig_2b$loc_name, levels = loc_name_order)

  # change alpha for Delta and others
  colors_lineage_type_fig2a <- colors_lineage_type
  colors_lineage_type_fig2a[grepl("Delta", names(colors_lineage_type_fig2a))] <- paste0(colors_lineage_type_fig2a[grepl("Delta", names(colors_lineage_type_fig2a))], "80")
  colors_lineage_type_fig2a[grepl("Others", names(colors_lineage_type_fig2a))] <- paste0(colors_lineage_type_fig2a[grepl("Others", names(colors_lineage_type_fig2a))], "80")

  # plotting 29 units vertically will be too long, 
  # we consider to plot 13 to the left and 16 to the right
  saveRDS(df_fig_2b, "results/figs/model_simulation/Omicron20/df_fig_2b_source.rds")

  p_fig_2b_left <- ggplot(df_fig_2b %>% filter(loc_name %in% loc_name_order[1:13])) +
    # geom_point(aes(x = earliest_date, y = cases_level, color = color_type, fill = color_type, size = n_simulation, group = color_type), position = position_dodge(width = 1), alpha=0.7) +
    geom_boxplot(aes(x = earliest_date, y = cases_level, color = color_type, fill = color_type, group = interaction(cases_level, color_type)), position = position_dodge(width = 1, preserve = "single"), outlier.shape = NA, show.legend = FALSE)+
    stat_summary(aes(x = earliest_date, y = cases_level, group = interaction(cases_level, color_type)), position = position_dodge(width = 1, preserve = "single"), geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    theme_minimal(base_family = "", base_size = 12) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_size(range = c(0.7, 2)) +
    scale_color_manual(name="", values = colors_lineage_type_fig2a, drop = FALSE) +
    scale_fill_manual(name="", values = colors_lineage_type_fig2a, drop = FALSE) +
    scale_alpha(range = c(0.5, 0.9)) +
    facet_wrap(
      vars(loc_name),
      strip.position = "top",
      ncol = 1
    )+
    theme(
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1),
      strip.text.y = element_text(angle = 0),
      axis.title = element_blank(),
      legend.position = "none"
    )+
    guides(alpha = "none")
  ggsave(paste0(dir_rst, "/fig_2b_left.jpg"), plot=p_fig_2b_left, width=6*2, height=6*3.5, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2b_left.pdf"), plot=p_fig_2b_left, width=6*2, height=6*3.5, dpi=400)

  p_fig_2b_right <- ggplot(df_fig_2b %>% filter(loc_name %in% loc_name_order[14:29])) +
    # geom_point(aes(x = earliest_date, y = cases_level, color = color_type, fill = color_type, size = n_simulation, group = color_type), position = position_dodge(width = 1), alpha=0.7) +
    geom_boxplot(aes(x = earliest_date, y = cases_level, color = color_type, fill = color_type, group = interaction(cases_level, color_type)), position = position_dodge(width = 1, preserve = "single"), outlier.shape = NA, show.legend = FALSE)+
    stat_summary(aes(x = earliest_date, y = cases_level, group = interaction(cases_level, color_type)), position = position_dodge(width = 1, preserve = "single"), geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% select(-cases_level) %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    theme_minimal(base_family = "", base_size = 12) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_size(range = c(0.7, 2)) +
    scale_color_manual(name="", values = colors_lineage_type_fig2a, drop = FALSE) +
    scale_fill_manual(name="", values = colors_lineage_type_fig2a, drop = FALSE) +
    scale_alpha(range = c(0.5, 0.9)) +
    facet_wrap(
      vars(loc_name),
      strip.position = "top",
      ncol = 1
    )+
    theme(
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1),
      strip.text.y = element_text(angle = 0),
      axis.title = element_blank(),
      # legend.box.margin = margin(t = 30),
      legend.position="none"
    )+
    guides(alpha = "none")
  ggsave(paste0(dir_rst, "/fig_2b_right.jpg"), plot=p_fig_2b_right, width=6*2, height=6*3.5, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2b_right.pdf"), plot=p_fig_2b_right, width=6*2, height=6*3.5, dpi=400)

  return(list(p_fig_2b_left, p_fig_2b_right))
}

# Fig 2c is to show the origins of infected cases
plot_fig_2c <- function(
  model_sims,
  dir_rst
){
  # select the imported infected cases in simulation results
  df_sims <- model_sims %>% select(date_decimal, `.id`, unitname, matches("C_[[:digit:]]_i_infected_origin_[[:digit:]]"))
  df_sims <- df_sims %>% filter(date_decimal <= decimal_date(ymd("2022-02-28"))) # only show the results before 2022-02-28

  # take the median of the simulation results for all N replicates
  df_sims <- df_sims %>% group_by(date_decimal, unitname) %>% summarise_all(median) %>% ungroup()
  df_sims_long <- df_sims %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value")
  df_sims_long <- df_sims_long %>% filter(value >= 1) %>% mutate(value=round(value))
  df_sims_long <- df_sims_long %>% separate(variable, into = c("cases", "variants", "cases_type", "cases_level", "origin", "origin_id"), sep = "_") 
  df_sims_long <- df_sims_long %>% filter(cases_type == "i", cases_level == "infected", value>=1) # only show the infected cases
  df_sims_long <- df_sims_long %>% filter(variants %in% c("2", "3")) %>% mutate(variants = factor(variants, levels=2:3, labels=c("Omicron BA.1", "Omicron BA.2")))
  df_sims_long$origin_id <- factor(as.numeric(df_sims_long$origin_id), levels=0:28, labels = sort(unique(model_sims$unitname))) %>% as.character()

  df_sims_long <- df_sims_long %>% group_by(date_decimal, variants, origin_id) %>% summarise(value = sum(value)) %>% ungroup()
  df_sims_long$date <- date_decimal(df_sims_long$date_decimal)

  saveRDS(df_sims_long, "results/figs/model_simulation/Omicron20/df_plot_2c_source.rds")

  # only show the top 10 origins for each unit, the rest origins are grouped as "Others"
  df_plot_2c_upper <- df_sims_long %>% filter(variants=="Omicron BA.1")
  df_plot_2c_upper_new_origin_id <- df_plot_2c_upper %>% group_by(variants, origin_id) %>% summarise(value = sum(value)) %>% ungroup() %>% arrange(variants, desc(value)) %>% group_by(variants) %>% mutate(rank = row_number(), new_origin_id = ifelse(rank>10, "Others", origin_id)) %>% ungroup() %>% select(variants, origin_id, new_origin_id) %>% arrange(variants, origin_id)
  df_plot_2c_upper <- df_plot_2c_upper %>% left_join(df_plot_2c_upper_new_origin_id, by=c("variants", "origin_id"))
  df_plot_2c_upper <- df_plot_2c_upper %>% group_by(date, variants, new_origin_id) %>% summarise(value = sum(value)) %>% ungroup()

  # add flag code
  df_plot_2c_upper$code_flag <- tolower(df_plot_2c_upper$new_origin_id)

  df_plot_2c_upper$new_origin_id <- factor(df_plot_2c_upper$new_origin_id, levels=df_plot_2c_upper %>% group_by(new_origin_id) %>% summarise(value = sum(value)) %>% ungroup() %>% arrange(desc(value)) %>% pull(new_origin_id))
  dates_flag_labels <- df_plot_2c_upper %>% group_by(new_origin_id) %>% filter(value==max(value)) %>% pull(date) %>% unique()
  df_plot_2c_flag_upper <- df_plot_2c_upper %>% filter(date %in% dates_flag_labels) %>% group_by(date) %>% arrange(desc(new_origin_id)) %>% mutate(y_flag=cumsum(value)-value/2) %>% ungroup() %>% group_by(new_origin_id) %>% slice_max(value, n=1, with_ties=FALSE) %>% ungroup() %>% select(date, new_origin_id, y_flag, code_flag)

  # Add introduction dates as vline
  lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
  df_vline <- tibble(lag=lags_intro)
  df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(df_sims$date_decimal)
  df_vline$date <- date_decimal(df_vline$date_decimal)
  df_vline$v_line_label <- c("Omicron BA.1 emergence", "Omicron BA.2 emergence")
  df_vline$color <- c("#25567d", "#357933")
  df_vline$v_line_label_y <- max(df_plot_2c_flag_upper$y_flag)
  
  p_fig_2c_upper <- df_plot_2c_upper %>% ggplot(
      aes(
        x=date, y=value, 
        color = new_origin_id, 
        fill = new_origin_id
      )
    ) +
    geom_area(
      stat = "identity",
      size=0,
      show.legend = FALSE
    )+
    geom_flag(
      data = df_plot_2c_flag_upper %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, y=y_flag),
      size=3
    )+
    geom_label(
      data = df_plot_2c_flag_upper %>% filter(grepl("others", code_flag)),
      aes(x=date-days(1), y=y_flag, label = new_origin_id),
      alpha=0.8,
      fill="white",
      size=2
    )+
    facet_wrap( ## needs facet_grid for space argument
      vars(variants), 
      ncol=1,
      scales = "fixed",
      strip.position = "top"
    ) +
    scale_color_manual(
      expand = c(0, 0),
      values = c("Others"="grey", colors_spatial_units),
      guide = "none"
    ) +
    scale_fill_manual(
      values = c("Others"="grey", colors_spatial_units),
      name = NULL,
      guide = "none"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_label_repel(
      data = df_vline %>% filter(color=="#25567d"),
      aes(x=date, y=v_line_label_y, label = v_line_label),
      color = "#25567d",
      size = 3,
      hjust = 0,
      inherit.aes = FALSE
    ) +
    geom_label_repel(
      data = df_vline %>% filter(color=="#357933"),
      aes(x=date, y=v_line_label_y, label = v_line_label),
      color = "#357933",
      size = 3,
      hjust = 1,
      inherit.aes = FALSE
    ) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_y_continuous(name="Cases", label=scientific_10) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "", base_size = 12) +
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = .5,
        margin = margin(10, 0, 30, 0)
      ),
      plot.caption = element_text(
        size = 9,
        color = "black",
        hjust = .5,
        margin = margin(20, 0, 5, 0)
      ),
      # axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
      axis.title.x = element_blank(),
      # plot.background = element_rect(fill = "grey88", color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      legend.position = 'none'
      # strip.text.y = element
      )
  ggsave(paste0(dir_rst, "/fig_2c_upper.jpg"), plot=p_fig_2c_upper, width=6*2, height=6, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2c_upper.pdf"), plot=p_fig_2c_upper, width=6*2, height=6, dpi=400)

  # only show the top 10 origins for each unit, the rest origins are grouped as "Others"
  df_plot_2c_lower <- df_sims_long %>% filter(variants=="Omicron BA.2")
  df_plot_2c_lower_new_origin_id <- df_plot_2c_lower %>% group_by(variants, origin_id) %>% summarise(value = sum(value)) %>% ungroup() %>% arrange(variants, desc(value)) %>% group_by(variants) %>% mutate(rank = row_number(), new_origin_id = ifelse(rank>10, "Others", origin_id)) %>% ungroup() %>% select(variants, origin_id, new_origin_id) %>% arrange(variants, origin_id)
  df_plot_2c_lower <- df_plot_2c_lower %>% left_join(df_plot_2c_lower_new_origin_id, by=c("variants", "origin_id"))
  df_plot_2c_lower <- df_plot_2c_lower %>% group_by(date, variants, new_origin_id) %>% summarise(value = sum(value)) %>% ungroup()

  # add flag code
  df_plot_2c_lower$code_flag <- tolower(df_plot_2c_lower$new_origin_id)

  df_plot_2c_lower$new_origin_id <- factor(df_plot_2c_lower$new_origin_id, levels=df_plot_2c_lower %>% group_by(new_origin_id) %>% summarise(value = sum(value)) %>% ungroup() %>% arrange(desc(value)) %>% pull(new_origin_id))
  dates_flag_labels <- df_plot_2c_lower %>% group_by(new_origin_id) %>% filter(value==max(value)) %>% pull(date) %>% unique()
  df_plot_2c_flag_lower <- df_plot_2c_lower %>% filter(date %in% dates_flag_labels) %>% group_by(date) %>% arrange(desc(new_origin_id)) %>% mutate(y_flag=cumsum(value)-value/2) %>% ungroup() %>% group_by(new_origin_id) %>% slice_max(value, n=1, with_ties=FALSE) %>% ungroup() %>% select(date, new_origin_id, y_flag, code_flag)

  p_fig_2c_lower <- df_plot_2c_lower %>% ggplot(
      aes(
        x=date, y=value, 
        color = new_origin_id, 
        fill = new_origin_id
      )
    ) +
    geom_area(
      stat = "identity",
      size=0,
      show.legend = FALSE
    )+
    geom_flag(
      data = df_plot_2c_flag_lower %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, y=y_flag),
      size=3
    )+
    geom_label(
      data = df_plot_2c_flag_lower %>% filter(grepl("others", code_flag)),
      aes(x=date-days(1), y=y_flag, label = new_origin_id),
      alpha=0.8,
      fill="white",
      size=2
    )+
    facet_wrap( ## needs facet_grid for space argument
      vars(variants), 
      ncol=1,
      scales = "fixed",
      strip.position = "top"
    ) +
    scale_color_manual(
      expand = c(0, 0),
      values = c("Others"="grey", colors_spatial_units),
      guide = "none"
    ) +
    scale_fill_manual(
      values = c("Others"="grey", colors_spatial_units),
      name = NULL,
      guide = "none"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_y_continuous(name="Cases", label=scientific_10) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "", base_size = 12) +
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = .5,
        margin = margin(10, 0, 30, 0)
      ),
      plot.caption = element_text(
        size = 9,
        color = "black",
        hjust = .5,
        margin = margin(20, 0, 5, 0)
      ),
      # axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
      axis.title.x = element_blank(),
      # plot.background = element_rect(fill = "grey88", color = NA),
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid = element_blank(),
      panel.spacing.y = unit(0, "lines"),
      legend.position = 'none'
      # strip.text.y = element
      )
  ggsave(paste0(dir_rst, "/fig_2c_lower.jpg"), plot=p_fig_2c_lower, width=6*2, height=6, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2c_lower.pdf"), plot=p_fig_2c_lower, width=6*2, height=6, dpi=400)

  return(list(p_fig_2c_upper, p_fig_2c_lower))
}

# Fig 2d is to show the origins of infected cases for every units, use linerange
## x-axis is date, base plot y-axis is the proportion of imported cases in all infections
## then overlay another linerange plot showing the dominant origins of imported cases
plot_fig_2d <- function(
  model_sims,
  dir_rst
){
  # aggregate the imported cases in simulation results
  df_sims <- model_sims %>%
    mutate(
      C_1_i_infected_new = rowSums(select(., starts_with("C_1_i_infected_origin_"))),
      C_1_i_detected_new = rowSums(select(., starts_with("C_1_i_detected_origin_"))),
      C_1_i_sequenced_new = rowSums(select(., starts_with("C_1_i_sequenced_origin_"))),
      C_2_i_infected_new = rowSums(select(., starts_with("C_2_i_infected_origin_"))),
      C_2_i_detected_new = rowSums(select(., starts_with("C_2_i_detected_origin_"))),
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_sequenced_origin_"))),
      C_3_i_infected_new = rowSums(select(., starts_with("C_3_i_infected_origin_"))),
      C_3_i_detected_new = rowSums(select(., starts_with("C_3_i_detected_origin_"))),
      C_3_i_sequenced_new = rowSums(select(., starts_with("C_3_i_sequenced_origin_"))),
      C_4_i_infected_new = rowSums(select(., starts_with("C_4_i_infected_origin_"))),
      C_4_i_detected_new = rowSums(select(., starts_with("C_4_i_detected_origin_"))),
      C_4_i_sequenced_new = rowSums(select(., starts_with("C_4_i_sequenced_origin_")))
    )

  df_sims <- df_sims %>% select(date_decimal, `.id`, unitname, matches("C_[[:digit:]]_i_infected_origin_[[:digit:]]"), matches("C_[[:digit:]]_[ic]_infected_new$")) # only select the infected cases

  df_sims <- df_sims %>% filter(date_decimal <= decimal_date(ymd("2022-02-28")))

  # take the median the simulation results for all N replicates
  df_sims <- df_sims %>% group_by(date_decimal, unitname) %>% summarise_all(median) %>% ungroup()

  # Add introduction dates as vline
  lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
  df_vline <- tibble(lag=lags_intro)
  df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(df_sims$date_decimal)
  df_vline$color <- c("#25567d", "#357933")
  df_vline$date <- date_decimal(df_vline$date_decimal)

  # 1. Base plot: the proportion of imported/exportation cases in all infections
  ## 1.1. Importation proportion
  df_fig_2d_base_all_inf <- df_sims %>% select(date_decimal, unitname, matches("C_[[:digit:]]_[ic]_infected_new$")) %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value") %>% filter(grepl("^C_2_", variable) | grepl("^C_3_", variable)) %>% filter(value>=1) %>% mutate(value=round(value)) %>% pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>% filter((C_2_i_infected_new>=1) | (C_3_i_infected_new>=1)) %>% transmute(date_decimal, unitname, C_2_all_infected_new = (C_2_i_infected_new+C_2_c_infected_new), C_3_all_infected_new = (C_3_i_infected_new+C_3_c_infected_new))
  
  df_fig_2d_base_importation <- df_sims %>% select(date_decimal, unitname, matches("C_[[:digit:]]_[ic]_infected_new$")) %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value") %>% filter(value>=1) %>% mutate(value=round(value)) %>% filter(grepl("^C_2_", variable) | grepl("^C_3_", variable)) %>% pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>% filter((C_2_i_infected_new>=1) | (C_3_i_infected_new>=1)) %>% transmute(date_decimal, unitname, BA_1_imported_prop = C_2_i_infected_new/(C_2_i_infected_new+C_2_c_infected_new), BA_2_imported_prop = C_3_i_infected_new/(C_3_i_infected_new+C_3_c_infected_new), BA_1_imported = C_2_i_infected_new, BA_2_imported = C_3_i_infected_new)
  df_fig_2d_base_importation_long <- df_fig_2d_base_importation %>% select(-BA_1_imported, -BA_2_imported) %>% pivot_longer(cols = starts_with("BA_"), names_to = "group", values_to = "value") %>% filter(value>0)

  df_fig_2d_base_exportation <- df_sims %>% select(date_decimal, unitname, matches("C_[[:digit:]]_i_infected_origin_[[:digit:]]")) %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value") %>% separate(variable, into = c("cases", "variants", "cases_type", "cases_level", "origin", "origin_id"), sep = "_") %>% filter(value>=1) %>% mutate(value=round(value)) %>% filter(cases_type == "i", cases_level == "infected") %>% filter(variants %in% c("2", "3")) %>% mutate(variants = factor(variants, levels=2:3, labels=c("Omicron BA.1", "Omicron BA.2"))) %>% group_by(date_decimal, variants, origin_id) %>% summarise(value = sum(value)) %>% ungroup()
  df_fig_2d_base_exportation$unitname <- factor(df_fig_2d_base_exportation$origin_id, levels=0:28, labels = sort(unique(model_sims$unitname))) %>% as.character()
  df_fig_2d_base_exportation <- df_fig_2d_base_exportation %>% select(-origin_id) %>% pivot_wider(names_from = variants, values_from = value, values_fill = 0) %>% left_join(df_fig_2d_base_all_inf, by=c("date_decimal", "unitname")) %>% mutate(BA_1_exported_prop = `Omicron BA.1`/C_2_all_infected_new, BA_2_exported_prop = `Omicron BA.2`/C_3_all_infected_new) %>% select(date_decimal, unitname, BA_1_exported_prop, BA_2_exported_prop, BA_1_exported=`Omicron BA.1`, BA_2_exported=`Omicron BA.2`)
  df_fig_2d_base_exportation_long <- df_fig_2d_base_exportation %>% select(-BA_1_exported, -BA_2_exported) %>% pivot_longer(cols = starts_with("BA_"), names_to = "group", values_to = "value") %>% filter(value>0)

  df_fig_2d_base <- bind_rows(df_fig_2d_base_importation_long, df_fig_2d_base_exportation_long)
  df_fig_2d_base$date <- date_decimal(df_fig_2d_base$date_decimal)
  df_fig_2d_base <- df_fig_2d_base %>% mutate(variant = factor(str_extract(group, "BA_[[:digit:]]"), levels=c("BA_1", "BA_2"), labels=c("Omicron BA.1", "Omicron BA.2")), group = factor(group, levels=c("BA_1_imported_prop", "BA_1_exported_prop", "BA_2_imported_prop", "BA_2_exported_prop"), labels=c("BA.1 Ori.", "BA.1 Dest.", "BA.2 Ori.", "BA.2 Dest.")))

  df_fig_2d_base$adjusted_y <- df_fig_2d_base$value
  df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.2 Ori."] <- df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.2 Ori."] + 1
  df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.1 Dest."] <- df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.1 Dest."] + 2
  df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.1 Ori."] <- df_fig_2d_base$adjusted_y[df_fig_2d_base$group=="BA.1 Ori."] + 3

  # split the figure into two parts, the left part should contain 13 units, the right part should contain 16 units
  ## use the same order as used in figure 2b
  df_fig_2b <- readRDS("results/figs/model_simulation/Omicron20/df_fig_2b_source.rds")
  order_units <- unique(df_fig_2b %>% arrange(loc_name) %>% pull(unitname))
  ## add formal unitname 
  df_loc_name <- left_join(tibble(unitname=order_units), cross_check_table %>% transmute(unitname=code, loc_name), by="unitname")
  df_loc_name$loc_name[is.na(df_loc_name$loc_name)] <- df_loc_name$unitname[is.na(df_loc_name$loc_name)]
  df_loc_name$loc_name[df_loc_name$loc_name=="Hong Kong Special Administrative Region of China"] <- "Hong Kong SAR"

  df_fig_2d_base <- df_fig_2d_base %>% group_by(unitname, group) %>% filter(n()>=1) %>% ungroup() # for the geom_xsplines to work, the data should have at least 2 points
  df_fig_2d_base_left <- df_fig_2d_base %>% filter(unitname %in% order_units[1:13])
  df_fig_2d_base_left$unitname <- factor(df_fig_2d_base_left$unitname, levels=df_loc_name$unitname[1:13], labels=df_loc_name$loc_name[1:13])
  df_fig_2d_base_right <- df_fig_2d_base %>% filter(unitname %in% order_units[14:29])
  df_fig_2d_base_right$unitname <- factor(df_fig_2d_base_right$unitname, levels=df_loc_name$unitname[14:29], labels=df_loc_name$loc_name[14:29])

  # 2. Overlay plot: the dominant origins of imported/exported cases
  df_fig_2d_overlay <- df_sims %>% select(date_decimal, unitname, matches("C_[[:digit:]]_i_infected_origin_[[:digit:]]")) %>% pivot_longer(cols = starts_with("C_"), names_to = "variable", values_to = "value") %>% separate(variable, into = c("cases", "variants", "cases_type", "cases_level", "origin", "origin_id"), sep = "_") %>% filter(cases_type == "i", cases_level == "infected", value>=1) %>% filter(variants %in% c("2", "3")) %>% mutate(variants = factor(variants, levels=2:3, labels=c("Omicron BA.1", "Omicron BA.2")))
  df_fig_2d_overlay$origin_id <- factor(as.numeric(df_fig_2d_overlay$origin_id), levels=0:28, labels = sort(unique(model_sims$unitname))) %>% as.character()
  df_fig_2d_overlay$date <- date_decimal(df_fig_2d_overlay$date_decimal) %>% round_date()
  df_fig_2d_overlay <- df_fig_2d_overlay %>% select(date, variants, unitname, origin_id, value)
  
  ## summarize the non-overlap time periods for each unit
  concatenate_intervals <- function(dates, date_type="start"){
    if(length(dates)==1){
      if(date_type=="start"){
      return(dates)
    } else if(date_type=="end"){
      return(dates+days(1))
    }
    }
    check_overlap <- sapply(seq_len(length(dates)-1), function(i){
      int_overlaps(interval(dates[i], dates[i]+days(1)), interval(dates[i+1], dates[i+1]+days(1)))
    })

    rle_rst <- rle(check_overlap)
    df_intervals <- lapply(seq_along(rle_rst$values), function(j){
      if(rle_rst$values[j]){ # with overlap
        if(j==1){
          dates_starts <- dates[1]
          dates_ends <- dates_starts + days(rle_rst$lengths[j])
        } else {
          dates_starts <- dates[cumsum(rle_rst$lengths)[j-1]+1]
          dates_ends <- dates_starts + days(rle_rst$lengths[j])
        }
        return(data.frame(start=dates_starts, end=dates_ends))
      } else { # without overlap
        if(j==1){
          dates_starts <- dates[1:rle_rst$lengths[j]]
        } else {
          dates_starts <- dates[(cumsum(rle_rst$lengths)[j-1]+1):cumsum(rle_rst$lengths)[j]]
        }
        dates_ends <- dates_starts + days(1)
        return(data.frame(start=dates_starts, end=dates_ends))
      }
    }) %>% bind_rows()
    if(date_type=="start"){
      return(df_intervals$start)
    } else if(date_type=="end"){
      return(df_intervals$end)
    }
  }

  ## Dominant origins of imported cases
  df_fig_2d_overlay_importation <- df_fig_2d_overlay %>% group_by(date, variants, unitname) %>% slice_max(value, with_ties=F) %>% ungroup() %>% transmute(date, variants, unitname, origin_id, value_importation=value)
  df_fig_2d_overlay_importation <- df_fig_2d_overlay_importation %>% group_by(variants, unitname, origin_id) %>% reframe(date_start_importation=concatenate_intervals(date, "start"), date_end_importation=concatenate_intervals(date, "end")) %>% ungroup()
  df_fig_2d_overlay_importation$adjusted_y <- 0
  df_fig_2d_overlay_importation$adjusted_y[df_fig_2d_overlay_importation$variants=="Omicron BA.1"] <- 3.5
  df_fig_2d_overlay_importation$adjusted_y[df_fig_2d_overlay_importation$variants=="Omicron BA.2"] <- 1.5
  df_fig_2d_overlay_importation$code_flag <- tolower(df_fig_2d_overlay_importation$origin_id)
  df_fig_2d_overlay_importation <- df_fig_2d_overlay_importation %>% group_by(variants, unitname) %>% arrange(variants, unitname, date_start_importation) %>% ungroup() 
  rle_overlay_imporatation <- df_fig_2d_overlay_importation %>% group_by(variants, unitname) %>% reframe(rle_value=rle(code_flag)[[1]], rle_code=rle(code_flag)[[2]]) # only show the first occurance if the next introduction is also the same unit
  df_fig_2d_overlay_importation_flags <- df_fig_2d_overlay_importation[c(1, cumsum(rle_overlay_imporatation$rle_value)[-nrow(rle_overlay_imporatation)]+1),]

  df_fig_2d_overlay_importation_left <- df_fig_2d_overlay_importation %>% filter(unitname %in% order_units[1:13])
  df_fig_2d_overlay_importation_left$unitname <- factor(df_fig_2d_overlay_importation_left$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_importation_flags_left <- df_fig_2d_overlay_importation_flags %>% filter(unitname %in% order_units[1:13])
  df_fig_2d_overlay_importation_flags_left$unitname <- factor(df_fig_2d_overlay_importation_flags_left$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_importation_right <- df_fig_2d_overlay_importation %>% filter(unitname %in% order_units[14:29])
  df_fig_2d_overlay_importation_right$unitname <- factor(df_fig_2d_overlay_importation_right$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_importation_flags_right <- df_fig_2d_overlay_importation_flags %>% filter(unitname %in% order_units[14:29])
  df_fig_2d_overlay_importation_flags_right$unitname <- factor(df_fig_2d_overlay_importation_flags_right$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)

  ## Dominant origins of exported cases
  df_fig_2d_overlay_exportation <- df_fig_2d_overlay %>% group_by(date, variants, origin_id) %>% slice_max(value, n=1, with_ties=F) %>% ungroup() %>% transmute(date, variants, code=origin_id, dest_id=unitname, value_exportation = value)
  df_fig_2d_overlay_exportation <- df_fig_2d_overlay_exportation %>% group_by(variants, code, dest_id) %>% reframe(date_start_exportation=concatenate_intervals(date, "start"), date_end_exportation=concatenate_intervals(date, "end")) %>% ungroup()
  df_fig_2d_overlay_exportation$adjusted_y <- 0
  df_fig_2d_overlay_exportation$adjusted_y[df_fig_2d_overlay_exportation$variants=="Omicron BA.1"] <- 2.5
  df_fig_2d_overlay_exportation$adjusted_y[df_fig_2d_overlay_exportation$variants=="Omicron BA.2"] <- 0.5
  df_fig_2d_overlay_exportation <- df_fig_2d_overlay_exportation %>% mutate(unitname=code) %>% select(-code)
  df_fig_2d_overlay_exportation$code_flag <- tolower(df_fig_2d_overlay_exportation$dest_id)
  df_fig_2d_overlay_exportation <- df_fig_2d_overlay_exportation %>% group_by(variants, unitname) %>% arrange(variants, unitname, date_start_exportation) %>% ungroup()
  rle_overlay_exportation <- df_fig_2d_overlay_exportation %>% group_by(variants, unitname) %>% reframe(rle_value=rle(code_flag)[[1]], rle_code=rle(code_flag)[[2]]) # only show the first occurance if the next introduction is also the same unit
  df_fig_2d_overlay_exportation_flags <- df_fig_2d_overlay_exportation[c(1, cumsum(rle_overlay_exportation$rle_value)[-nrow(rle_overlay_exportation)]+1),]

  df_fig_2d_overlay_exportation_left <- df_fig_2d_overlay_exportation %>% filter(unitname %in% order_units[1:13])
  df_fig_2d_overlay_exportation_left$unitname <- factor(df_fig_2d_overlay_exportation_left$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_exportation_flags_left <- df_fig_2d_overlay_exportation_flags %>% filter(unitname %in% order_units[1:13])
  df_fig_2d_overlay_exportation_flags_left$unitname <- factor(df_fig_2d_overlay_exportation_flags_left$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_exportation_right <- df_fig_2d_overlay_exportation %>% filter(unitname %in% order_units[14:29])
  df_fig_2d_overlay_exportation_right$unitname <- factor(df_fig_2d_overlay_exportation_right$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)
  df_fig_2d_overlay_exportation_flags_right <- df_fig_2d_overlay_exportation_flags %>% filter(unitname %in% order_units[14:29])
  df_fig_2d_overlay_exportation_flags_right$unitname <- factor(df_fig_2d_overlay_exportation_flags_right$unitname, levels=df_loc_name$unitname, labels=df_loc_name$loc_name)

  # plot left part
  p_fig_2d_overlay_left <- ggplot()+
    geom_linerange(
      data = df_fig_2d_overlay_importation_left,
      aes(xmin=date_start_importation, xmax=date_end_importation, y=adjusted_y, color=origin_id),
      size=4,
      alpha=1,
      show.legend = FALSE
    )+
    geom_linerange(
      data = df_fig_2d_overlay_exportation_left,
      aes(xmin=date_start_exportation, xmax=date_end_exportation, y=adjusted_y, color=dest_id),
      size=4,
      alpha=1,
      show.legend = FALSE
    )+
    # geom_xspline(
    #   data=df_fig_2d_base_left %>% filter(variant=="Omicron BA.1"), aes(x=date, y=adjusted_y, group=group), color=colors_lineage[["Omicron BA.1"]]
    # )+
    # geom_xspline(
    #   data=df_fig_2d_base_left %>% filter(variant=="Omicron BA.2"), aes(x=date, y=adjusted_y, group=group), color=colors_lineage[["Omicron BA.2"]]
    # )+
    facet_wrap(vars(unitname), ncol=1)+
    scale_y_continuous(
      name="Proportion in total infection cases", breaks=c(0, 1, 2, 3)+0.5, labels=rev(levels(df_fig_2d_base$group))
      )+
    geom_vline(
      data = df_vline %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_flag(
      data = df_fig_2d_overlay_importation_flags_left %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, x=date_start_importation, y=adjusted_y),
      size=3
    )+
    geom_flag(
      data = df_fig_2d_overlay_exportation_flags_left %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, x=date_start_exportation, y=adjusted_y),
      size=3
    )+
    geom_label(
      data = df_fig_2d_overlay_importation_flags_left %>% filter(grepl("others", code_flag)) %>% mutate(origin_id = gsub("_others", "", origin_id)),
      aes(x=date_start_importation, y=adjusted_y, label = origin_id),
      alpha=0.8,
      fill="white",
      label.padding = unit(0.1, "lines"),
      size=1.5
    )+
    geom_label(
      data = df_fig_2d_overlay_exportation_flags_left %>% filter(grepl("others", code_flag)) %>% mutate(dest_id = gsub("_others", "", dest_id)),
      aes(x=date_start_exportation, y=adjusted_y, label = dest_id),
      alpha=0.8,
      fill="white",
      label.padding = unit(0.1, "lines"),
      size=1.5
    )+
    theme_minimal(base_family = "", base_size = 12) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_size(range = c(0.7, 2)) +
    scale_color_manual(name="", values = colors_spatial_units, drop = FALSE) +
    scale_alpha(range = c(0.5, 0.9)) +
    theme(
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1),
      axis.title = element_blank(),
      strip.text.y = element_text(angle = 0),
      # legend.box.margin = margin(t = 30),
      panel.background = element_rect(fill = NA, color = NA),
      panel.spacing.y = unit(0, "lines"),
      legend.spacing.y = unit(0, "lines"),
      legend.position="none"
    )+
    guides(alpha = "none", color = guide_legend(nrow = 2))
  ggsave(paste0(dir_rst, "/fig_2d_overlay_left.jpg"), plot=p_fig_2d_overlay_left, width=6*2, height=6*3, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2d_overlay_left.pdf"), plot=p_fig_2d_overlay_left, width=6*2, height=6*3)

  # plot right part
  df_fig_2d_overlay_importation_right$origin_id <- factor(df_fig_2d_overlay_importation_right$origin_id, levels=sort(df_loc_name$unitname))
  p_fig_2d_overlay_right <- ggplot()+
    geom_linerange(
      data = df_fig_2d_overlay_importation_right,
      aes(xmin=date_start_importation, xmax=date_end_importation, y=adjusted_y, color=origin_id),
      size=4,
      alpha=1,
      show.legend=TRUE
    )+
    geom_linerange(
      data = df_fig_2d_overlay_exportation_right,
      aes(xmin=date_start_exportation, xmax=date_end_exportation, y=adjusted_y, color=dest_id),
      size=4,
      alpha=1,
      show.legend=TRUE
    )+
    # geom_xspline(
    #   data=df_fig_2d_base_right %>% filter(variant=="Omicron BA.1"), aes(x=date, y=adjusted_y, group=group), color=colors_lineage[["Omicron BA.1"]]
    #   )+
    # geom_xspline(
    #   data=df_fig_2d_base_right %>% filter(variant=="Omicron BA.2"), aes(x=date, y=adjusted_y, group=group), color=colors_lineage[["Omicron BA.2"]]
    #   )+
    facet_wrap(vars(unitname), ncol=1)+
    scale_y_continuous(
      name="Proportion in total infection cases", breaks=c(0, 1, 2, 3)+0.5, labels=rev(levels(df_fig_2d_base$group))
      )+
    geom_vline(
      data = df_vline %>% filter(color=="#25567d"),
      aes(xintercept = date),
      color = "#25567d", 
      size = 0.5,
      linetype = "dotted"
    ) +
    geom_vline(
      data = df_vline %>% filter(color=="#357933"),
      aes(xintercept = date),
      color = "#357933", 
      size = 0.5,
      linetype = "dotted"
    ) +

    geom_flag(
      data = df_fig_2d_overlay_importation_flags_right %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, x=date_start_importation, y=adjusted_y),
      size=3
    )+
    geom_flag(
      data = df_fig_2d_overlay_exportation_flags_right %>% filter(!grepl("others", code_flag)),
      aes(country=code_flag, x=date_start_exportation, y=adjusted_y),
      size=3
    )+
    geom_label(
      data = df_fig_2d_overlay_importation_flags_right %>% filter(grepl("others", code_flag)) %>% mutate(origin_id = gsub("_others", "", origin_id)),
      aes(x=date_start_importation, y=adjusted_y, label = origin_id),
      alpha=0.8,
      fill="white",
      label.padding = unit(0.1, "lines"),
      size=1.5
    )+
    geom_label(
      data = df_fig_2d_overlay_exportation_flags_right %>% filter(grepl("others", code_flag)) %>% mutate(dest_id = gsub("_others", "", dest_id)),
      aes(x=date_start_exportation, y=adjusted_y, label = dest_id),
      alpha=0.8,
      fill="white",
      label.padding = unit(0.1, "lines"),
      size=1.5
    )+
    theme_minimal(base_family = "", base_size = 12) +
    scale_x_datetime(
      date_breaks = "1 week",
      date_labels = "%y %b %d",
      expand = c(0, 0),
      limits = c(min(df_vline$date)-days(1), ymd("2022-02-28"))
    ) +
    scale_size(range = c(0.7, 2)) +
    scale_color_manual(name="", values = colors_spatial_units, drop = FALSE) +
    scale_alpha(range = c(0.5, 0.9)) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1),
      strip.text.y = element_text(angle = 0),
      # legend.box.margin = margin(t = 30),
      # panel.background = element_rect(fill = NA, color = NA),
      panel.spacing.y = unit(0, "lines"),
      legend.position="bottom"
    )+
    guides(
      alpha = "none", 
      color = guide_legend(
        nrow = 2, 
        # keywidth=0.1,
        keyheight=0.1,
        default.unit="inch")) +
    NULL
  ggsave(paste0(dir_rst, "/fig_2d_overlay_right.jpg"), plot=p_fig_2d_overlay_right, width=6*2, height=6*3, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2d_overlay_right.pdf"), plot=p_fig_2d_overlay_right, width=6*2, height=6*3)

  return(list(p_fig_2d_overlay_left, p_fig_2d_overlay_right))
}

plot_fig_2 <- function(
  model_sims,
  dir_rst
){
  p_fig_2a <- plot_fig_2a(model_sims, dir_rst)
  list_fig_2b <- plot_fig_2b(model_sims, dir_rst)

  p_fig_2_left <- ((p_fig_2a / list_fig_2b[[1]] + plot_layout(heights = c(3,14.5))) | list_fig_2b[[2]]) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = 'bottom')
  ggsave(paste0(dir_rst, "/fig_2_left.jpg"), plot=p_fig_2_left, width=5*3, height=5*4, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2_left.pdf"), plot=p_fig_2_left, width=5*3, height=5*4, dpi=400)

  list_fig_2c <- plot_fig_2c(model_sims, dir_rst)
  p_fig_2c <- list_fig_2c[[1]] + list_fig_2c[[2]] + plot_layout(ncol=1, guides = "collect", axes = "collect") & theme(legend.position = 'none')
  ggsave(paste0(dir_rst, "/fig_2c.jpg"), plot=p_fig_2c, width=3*4, height=3*3, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2c.pdf"), plot=p_fig_2c, width=5*3, height=5*4, dpi=400)

  list_fig_2d <- plot_fig_2d(model_sims, dir_rst)
  p_fig_2_right <- ((p_fig_2c / list_fig_2d[[1]] + plot_layout(heights = c(1.3, 1.3, 14.5))) | list_fig_2d[[2]]) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = 'bottom')
  ggsave(paste0(dir_rst, "/fig_2_right.jpg"), plot=p_fig_2_right, width=5*3, height=5*4, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2_right.pdf"), plot=p_fig_2_right, width=5*3, height=5*4, dpi=400)

  p_fig_2 <- (p_fig_2_left + p_fig_2_right) + plot_layout(nrow = 1, widths = c(1, 1, 2.1))
  ggsave(paste0(dir_rst, "/fig_2.jpg"), plot=p_fig_2, width=5*6, height=5*4, dpi=400)
  ggsave(paste0(dir_rst, "/fig_2.pdf"), plot=p_fig_2, width=5*6, height=5*4, dpi=400)

}
