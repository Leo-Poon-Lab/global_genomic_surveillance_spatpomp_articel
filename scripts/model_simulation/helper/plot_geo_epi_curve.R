library(geofacet)
library(ggflags)
source("scripts/model_fitting/helper/sim_cases.R")

plot_geo_epi_curve <- function(
  data_quants,
  df_meas,
  events,
  lineages_for_events,
  colors_events,
  colors_labels,
  geo_grid=geo_grid,
  add_flags=FALSE
){
  lineages_for_events <- gsub(" ", "_", lineages_for_events)
  names(colors_events) <- gsub(" ", "_", names(colors_events))

  df_meas <- aggregating_measurements(df_meas)
  df_meas$date <- as.Date(lubridate::date_decimal(df_meas$date_decimal))
  if(any(names(df_meas)=="daily_cases")){
    names(df_meas)[names(df_meas)=="daily_cases"] <- "Cases"
  }

  data_quants$code <- data_quants$unitname

  df_meas_plot <- df_meas %>% select(date, code, loc_name, all_of(lineages_for_events)) %>% pivot_longer(cols = c(-date, -code, -loc_name), names_to = "lineage", values_to = "Cases")

  # single line plot
  p <- ggplot(data = df_meas_plot %>% filter(Cases>1)) +
    geom_line(aes(x = date, y = Cases, group = lineage, color = lineage), linewidth = 0.7, linetype="solid", alpha=0.9) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10(  # Plot on log-scale
      limits = c(1, 10^7.5),
      breaks = c(1, 10, 10^3, 10^5, 10^7),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = c(0, 0)
    ) +
    scale_color_manual(values = colors_events, labels = colors_labels) +
    NULL

  # ribbon bands
  for (i in seq_along(events)){
    this_event <- events[i]
    data_quants_tmp <- data_quants %>% 
      dplyr::select(code, date, quantile_level, all_of(this_event)) %>% 
      filter(date!=min(data_quants$date)) %>% 
      mutate_at(vars(all_of(this_event)), ~ifelse(.<1, 1, .)) %>%
      pivot_wider(names_from = quantile_level, values_from = this_event, values_fn = mean)
    p <- p +
      geom_ribbon(data = data_quants_tmp, aes(x = date, ymin = get("0.05"), ymax = get("0.95")), fill = colors_events[i], alpha = 0.3) +
      NULL
  }

  # text of region
  df_meas_plot_text <- df_meas_plot %>% filter(date==min(df_meas_plot$date), lineage=="Cases")
  df_meas_plot_text <- left_join(df_meas_plot_text, geo_grid %>% select(code, name), by="code")
  p <- p +
    geom_label(data = df_meas_plot_text, aes(x = date, y = 10^6.5, label = name), hjust = 0, vjust = 0, size = 2.5, nudge_x = 5, color=colors_spatial_units[order(names(colors_spatial_units))]) +
    NULL

  df_meas_plot_text$code_flags <- tolower(df_meas_plot_text$code)
  if(add_flags){
    p <- p +
      geom_flag(data = df_meas_plot_text %>% filter(!grepl("others", code)), aes(x = date+170, y = 10^7, country = code_flags), size = 5) +
      NULL
  }

  p <- p + 
    facet_geo(~code, grid = geo_grid, label="name") +
    theme_minimal() +
    theme(
      # Top-right position
      legend.position = "bottom",
      # Elements within a guide are placed one next to the other in the same row
      legend.direction = "horizontal",
      # Different guides are stacked vertically
      legend.box = "vertical",
      # No legend title
      legend.title = element_blank(),
      # Light background color
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 30, 20, 30),
      # # Customize the title. Note the new font family and its larger size.
      # plot.title = element_text(
      #   margin = margin(0, 0, -100, 0), 
      #   size = 26, 
      #   family = "KyivType Sans", 
      #   face = "bold", 
      #   vjust = 0, 
      #   color = "grey25"
      # ),
      plot.caption = element_text(size = 11),
      # Remove titles for x and y axes.
      axis.title = element_blank(),
      # Specify color for the tick labels along both axes 
      axis.text = element_text(color = "grey40"),
      # Specify face and color for the text on top of each panel/facet
      # strip.text = element_text(color = "grey20")
      strip.background = element_blank(),
      strip.text = element_blank()
    ) + 
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

}

get_geo_grid <- function(
  num_country=86,
  regions_code_3,
  other_codes
  ){
  # get geo grid
  if(num_country == 86){
    data_geo_grid <- read_csv("https://raw.githubusercontent.com/hafen/grid-designer/master/grids/world_86countries_grid.csv")
  } else if(num_country == 192){
    data_geo_grid <- read_csv("https://raw.githubusercontent.com/hafen/grid-designer/master/grids/world_countries_grid1.csv")
    names(data_geo_grid)[2] <- "code"
  } else{
    stop("num_country must be either 86 or 192")
  }
  data_geo_grid_sub <- data_geo_grid %>% filter(code %in% regions_code_3)
  if(any(!regions_code_3 %in% data_geo_grid$code)){
    code_to_add <- regions_code_3[!regions_code_3 %in% data_geo_grid$code]
    data_geo_grid_sub <- bind_rows(data_geo_grid_sub, tibble(code=code_to_add, name=code_to_add, row=max(data_geo_grid$row)+1, col=max(data_geo_grid$col)+seq_along(code_to_add)))
  }
  if(all(!is.na(other_codes))){
    data_geo_grid_sub <- bind_rows(data_geo_grid_sub, tibble(code=other_codes, name=other_codes, row=max(data_geo_grid$row)+2, col=max(data_geo_grid$col)+seq_along(other_codes)))
  }
  ## change row to be sequential
  old_row <- sort(unique(data_geo_grid_sub$row))
  new_row <- 1:length(old_row)
  data_geo_grid_sub <- data_geo_grid_sub %>% mutate(row = new_row[match(row, old_row)])
  ## change col to be sequential
  old_col <- sort(unique(data_geo_grid_sub$col))
  new_col <- 1:length(old_col)
  data_geo_grid_sub <- data_geo_grid_sub %>% mutate(col = new_col[match(col, old_col)])

  return(data_geo_grid_sub)
}