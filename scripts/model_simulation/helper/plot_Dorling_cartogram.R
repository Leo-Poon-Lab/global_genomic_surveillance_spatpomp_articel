library(rnaturalearth)
library(sf)
library(tidyverse)
library(cartogram)
library(ggforce)
library(packcircles)
library(ggrepel)

cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

# change names for C cases:
change_compartment_names <- function(x){
  switch(x,
    "C_1_sequenced"="Total Delta sequenced",
    "C_1c_detected"="Community Delta detected",
    "C_1c_new"="Community Delta new",
    "C_1c_sequenced"="Community Delta sequenced",
    "C_1i_detected"="Imported Delta detected",
    "C_1i_new"="Imported Delta new",
    "C_1i_sequenced"="Imported Delta sequenced",
    "C_2_sequenced"="Total BA.1 sequenced",
    "C_2c_detected"="Community BA.1 detected",
    "C_2c_new"="Community BA.1 new",
    "C_2c_sequenced"="Community BA.1 sequenced",
    "C_2i_detected"="Imported BA.1 detected",
    "C_2i_new"="Imported BA.1 new",
    "C_2i_sequenced"="Imported BA.1 sequenced",
    "C_3_sequenced"="Total BA.2 sequenced",
    "C_3c_detected"="Community BA.2 detected",
    "C_3c_new"="Community BA.2 new",
    "C_3c_sequenced"="Community BA.2 sequenced",
    "C_3i_detected"="Imported BA.2 detected",
    "C_3i_new"="Imported BA.2 new",
    "C_3i_sequenced"="Imported BA.2 sequenced",
    "C_4_sequenced"="Total Others sequenced",
    "C_4c_detected"="Community Others detected",
    "C_4c_new"="Community Others new",
    "C_4c_sequenced"="Community Others sequenced",
    "C_4i_detected"="Imported Others detected",
    "C_4i_new"="Imported Others new",
    "C_4i_sequenced"="Imported Others sequenced",
    "C_new"="Total cases",
    "C_reported_new"="Total cases reported",
    "C_new_sequenced"="Total cases sequenced"
  )
}

plot_dorling <- function(
  data,
  events,
  colors_events,
  quants_range,
  rel_circle_size = 0.8
){
  # data = df_EDT_quants_lag
  # events = c("C_2c_sequenced", "C_2i_sequenced", "C_3c_sequenced", "C_3i_sequenced")

  data_median <- data %>% filter(quantile=="median") %>% select(unitname, any_of(events))
  data_upper <- data %>% filter(quantile=="upper") %>% select(unitname, any_of(events))
  data_lower <- data %>% filter(quantile=="lower") %>% select(unitname, any_of(events))

  max_lag <- data %>% select(starts_with("C_")) %>% unlist() %>% max(na.rm=T)
  data_median <- data_median %>% mutate_at(vars(starts_with("C_")), function(x){max_lag-x}) # inverse lag by max_lag
  data_inner <- data_upper %>% mutate_at(vars(starts_with("C_")), function(x){max_lag-x}) # inverse lag by max_lag
  data_outer <- data_lower %>% mutate_at(vars(starts_with("C_")), function(x){max_lag-x}) # inverse lag by max_lag

  world <- ne_countries(scale = 110, type = "countries", returnclass = "sf")%>%
  # Convert WGS84 to projected crs (here Robinson)
  sf::st_transform(world_ne, crs="ESRI:54030")

  list_new_region <- c("AF_others", "AS_others", "EU_others", "NA_others", "OC_others", "SA_others", "HK")
  all_codes_2 <- sort(unique(data$unitname))
  all_position <- lapply(seq_along(all_codes_2), function(i){
      # Get geospatial data
      if(all_codes_2[i]=="AF_others"){
        return(c(28, 11))
      } else if(all_codes_2[i]=="AS_others"){
        return(c(99, 15))
      } else if(all_codes_2[i]=="EU_others"){
        return(c(15, 62))
      } else if(all_codes_2[i]=="NA_others"){
        return(c(-87, 24))
      } else if(all_codes_2[i]=="OC_others"){
        return(c(161, -22))
      } else if(all_codes_2[i]=="SA_others"){
        return(c(-66, -31))
      } else if(all_codes_2[i]=="HK"){
        return(c(114.17, 22.29))
      }
      data_country <- world[which(world$iso_a2_eh == all_codes_2[i]), ]
      point <- c(data_country$label_x[1], data_country$label_y[1])
  })
  data_pos <- tibble(code_2=all_codes_2, lon=sapply(all_position, function(x) {x[1]}), lat=sapply(all_position, function(x) {x[2]}))
  sf_others <- data_pos %>% filter(code_2 %in% list_new_region) %>% .[,c("lon", "lat")] %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% sf::st_transform(crs="ESRI:54030")
  sf_others$iso_a2_eh <- data_pos$code_2[data_pos$code_2 %in% list_new_region]

  map_data <- world %>% left_join(data_pos %>% mutate(iso_a2_eh=code_2) %>% select(-code_2))
  map_data <- map_data %>% filter(!is.na(lon))
  map_data <- map_data %>% select(iso_a2_eh, geometry)
  map_data <- rbind(map_data, sf_others)

  centr <- map_data%>%
    st_centroid()%>%
    st_coordinates()

  map_data <- bind_cols(map_data, centr)

  # reorder data_outer and data_inner per map_data
  data_outer <- data_outer[match(map_data$iso_a2_eh, data_outer$unitname),]
  data_inner <- data_inner[match(map_data$iso_a2_eh, data_inner$unitname),]
  data_median <- data_median[match(map_data$iso_a2_eh, data_median$unitname),]

  ## circle parking for the background circle
  # Prepare data for circle packing
  map_data <- left_join(map_data, data_outer %>% mutate(iso_a2_eh=unitname) %>% select(-unitname))
  packing_data <- data.frame(x = map_data$X, y = map_data$Y, radius = apply(map_data %>% as_tibble() %>% select(starts_with("C_")), 1, function(x){max(x, na.rm=T)}))
  packing_data$radius <- packing_data$radius*12000*rel_circle_size

  # Use the packcircles algorithm to find non-overlapping positions for the circles
  layout <- circleRepelLayout(packing_data, sizetype = "radius")$layout

  # Update map_data with new coordinates
  map_data$X_new <- map_data$X
  map_data$Y_new <- map_data$Y
  map_data$radius <- packing_data$radius
  map_data$X_new[!is.na(layout$x)] <- layout$x[!is.na(layout$x)]
  map_data$Y_new[!is.na(layout$y)] <- layout$y[!is.na(layout$y)]

  ## Set colors and theme
  col_world <- "#F5F4EF"
  col_back <- "#ffffff"
  theme_custom <- theme_void()+
    theme(plot.background = element_rect(fill=col_back,color=NA),
    panel.grid.major = element_line(color = alpha("grey", 0.5)))

  data_example_cricle = tibble(
    x = -12000000, y = -6000000, r = (max_lag-c(10, 50, 100))*12000*rel_circle_size, lags=paste0(c(10, 50, 100), " days"))
  data_example_cricle$x <- data_example_cricle$x + data_example_cricle$r
  data_example_cricle$x_label <- data_example_cricle$x + data_example_cricle$r

  p_base <- ggplot()+
    # World basemap
    geom_sf(
      world,mapping=aes(geometry=geometry),
      fill=col_world,color=alpha("dimgrey",0.25)
    )+
    # Draw Dorling cartogram with geom_circle()
    ggforce::geom_circle(
      data = map_data, aes(x0 = X_new, y0 = Y_new, r = radius),
      fill=alpha("dimgrey",0.15),color=NA, linewidth = 0.3
    )+
    ggforce::geom_circle(
      data = data_example_cricle, aes(x0 = x, y0 = y, r = r),
      color=alpha("black",0.5)
    )+
    geom_text_repel(data = data_example_cricle, aes(x = x_label, y = y, label=lags), size=2, color="black", alpha=0.8)+
    scale_x_continuous(breaks = seq(-180, 180, 25))+
    theme_custom

  ## all events for all countries
  n_points=100
  circleFun <- function(
    center=c(0,0),   # center of the circle 
    diameter=1,      # diameter 
    npoints,     # number of points to draw the circle
    start=0, end=2   # start point/end point
  ){
      tt <- seq(start*pi, end*pi, length.out=npoints)
      tb <- tibble(
        x = center[1] + diameter / 2 * cos(tt), 
        y = center[2] + diameter / 2 * sin(tt)
      )
    return(tb)
  }

  break_step <- 2/length(events)
  
  data_outer$unitname

  data_crop_outer <- lapply(seq_along(events), function(i){ # outer
    this_event=events[i]
    data_event_crop <- lapply(seq_len(nrow(map_data)), function(j){
      # j=1
      df_points_median <- circleFun( # for median
        c(map_data$X_new[j],map_data$Y_new[j]), data_median[[this_event]][j]*12000*rel_circle_size*2, start=break_step*i, end=break_step*(i+1), npoints=n_points
      )
      df_points_outer <- circleFun( # for outer
        c(map_data$X_new[j],map_data$Y_new[j]), data_outer[[this_event]][j]*12000*rel_circle_size*2, start=break_step*i, end=break_step*(i+1), npoints=n_points
      )
      if(is.na(data_outer[[this_event]][j])){
        return(NA)
      } else {
        return(bind_cols(
          iso_a2_eh = rep(map_data$iso_a2_eh[j],n_points*2),
          bind_rows(df_points_median, df_points_outer[rev(seq_len(nrow(df_points_outer))),])
        ))
      }
    })
    bind_cols(bind_rows(data_event_crop[!is.na(data_event_crop)]), event=this_event)
  }) %>% bind_rows()
  data_crop_outer$event <- factor(data_crop_outer$event, levels=events)

  data_crop_inner <- lapply(seq_along(events), function(i){ # inner
    this_event=events[i]
    data_event_crop <- lapply(seq_len(nrow(map_data)), function(j){
      # j=1
      df_points_median <- circleFun( # for median
        c(map_data$X_new[j],map_data$Y_new[j]), data_median[[this_event]][j]*12000*rel_circle_size*2, start=break_step*i, end=break_step*(i+1), npoints=n_points
      )
      df_points_inner <- circleFun( # for inner
        c(map_data$X_new[j],map_data$Y_new[j]), data_inner[[this_event]][j]*12000*rel_circle_size*2, start=break_step*i, end=break_step*(i+1), npoints=n_points
      )
      if(is.na(data_median[[this_event]][j])){return(NA)}
      if(is.na(data_inner[[this_event]][j])){ # suggesting infinite detection lag
        return(bind_cols(
          iso_a2_eh = rep(map_data$iso_a2_eh[j],n_points+1),
          bind_rows(df_points_median, tibble(x=map_data$X_new[j], y=map_data$Y_new[j]))
        ))
      } else {
        return(bind_cols(
          iso_a2_eh = rep(map_data$iso_a2_eh[j],n_points*2),
          bind_rows(df_points_inner[rev(seq_len(nrow(df_points_inner))),], df_points_median)
        ))
      }
    })
    bind_cols(bind_rows(data_event_crop[!is.na(data_event_crop)]), event=this_event)
  }) %>% bind_rows()
  data_crop_inner$event <- factor(data_crop_inner$event, levels=events)

  data_crop_median <- lapply(seq_along(events), function(i){ # inner
    this_event=events[i]
    data_event_crop <- lapply(seq_len(nrow(map_data)), function(j){
      # j=1
      df_points_median <- circleFun( # for median
        c(map_data$X_new[j],map_data$Y_new[j]), data_median[[this_event]][j]*12000*rel_circle_size*2, start=break_step*i, end=break_step*(i+1), npoints=n_points
      )
      if(is.na(data_median[[this_event]][j])){return(NA)}
      return(bind_cols(
          iso_a2_eh = rep(map_data$iso_a2_eh[j],n_points),
          df_points_median
        ))
    })
    bind_cols(bind_rows(data_event_crop[!is.na(data_event_crop)]), event=this_event)
  }) %>% bind_rows()
  data_crop_median$event <- factor(data_crop_median$event, levels=events)

  lengend_labels <- sapply(names(colors_events), change_compartment_names)
  map_data <- left_join(map_data %>% mutate(code=iso_a2_eh), cross_check_table %>% select(code, loc_name))

  map_data$loc_name[is.na(map_data$loc_name)] <- map_data$iso_a2_eh[is.na(map_data$loc_name)]
  map_data$loc_name[map_data$iso_a2_eh=="HK"] <- "Hong Kong SAR"

  p_add <- p_base + 
    geom_polygon(data=data_crop_outer, aes(x=x, y=y, group=interaction(iso_a2_eh, event), fill=event), color=NA, alpha=0.7)+
    geom_polygon(data=data_crop_inner, aes(x=x, y=y, group=interaction(iso_a2_eh, event), fill=event), color=NA, alpha=0.7)+
    geom_line(data=data_crop_median, aes(x=x, y=y, group=interaction(iso_a2_eh, event), color=event), alpha=0.7)+
    scale_fill_manual(values = colors_events, labels = lengend_labels, name = "Inverse day lag")+
    scale_color_manual(values = colors_events, labels = lengend_labels, name = "Inverse day lag")+
    shadowtext::geom_shadowtext(data = map_data, aes(x = X_new, y = Y_new, label=loc_name), size=2, color="black", alpha=0.8, bg.color = "white",
    bg.alpha = 0.5)+
    theme(
            # legend.justification defines the edge of the legend that the legend.position coordinates refer to
            legend.justification = c(0, 1),
            # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
            legend.position = c(0.1, 0.48)
            )+
    guides(
      colour = guide_legend(override.aes = list(alpha = 1)),
      fill = guide_legend(override.aes = list(alpha = 1)))+
    NULL
  p_add
}