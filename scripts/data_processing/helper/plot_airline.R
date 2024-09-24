library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library(ggrepel)
library(geosphere)

plot_global_country <- function(mov_mat=mov_mat){
    # mov_mat <- readRDS(here::here(paste0("data/estimated_movement_matrix/move_mat_", 2021, "_", "country", ".rds")))
    # mov_mat <- mov_mat %>% filter(date>="2021-09-01")

    mov_mat_summary <- mov_mat %>% group_by(code_o, code_d, data_source) %>% summarise(total_flight=sum(flow_all)) %>% ungroup()
    # mov_mat_summary <- mov_mat_aggregate %>% group_by(code_o, code_d, data_source) %>% summarise(total_flight=sum(flow_all)) %>% ungroup()
    mov_mat_summary$pair <- mapply(function(x,y){paste0(sort(c(x, y)), collapse="-")}, mov_mat_summary$code_o, mov_mat_summary$code_d)
    mov_mat_summary <- mov_mat_summary %>% group_by(pair, data_source) %>% summarise(total_flight_pair=sum(total_flight)) 
    mov_mat_summary$code_o <- sapply(mov_mat_summary$pair, function(x) {gsub("-\\D+$", "", x)})
    mov_mat_summary$code_d <- sapply(mov_mat_summary$pair, function(x) {gsub("^\\D+-", "", x)})

    spdf_world <- ne_countries()
    # spdf_world <- spdf_world@data
    all_codes_2 <- sort(unique(c(mov_mat_summary$code_o, mov_mat_summary$code_d)))
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
        data_country <- spdf_world[which(spdf_world$iso_a2_eh == all_codes_2[i]), ]
        point <- c(data_country$label_x[1], data_country$label_y[1])
    })
    df_pos <- tibble(code_2=all_codes_2, lon=sapply(all_position, function(x) {x[1]}), lat=sapply(all_position, function(x) {x[2]}))
    mov_mat_summary <- left_join(mov_mat_summary, df_pos %>% transmute(code_o=code_2, lon_o=lon, lat_o=lat), "code_o")
    mov_mat_summary <- left_join(mov_mat_summary, df_pos %>% transmute(code_d=code_2, lon_d=lon, lat_d=lat), "code_d")

    # quantile(mov_mat_summary$total_flight_pair)

    ## Plot flight routes
    world <- ne_countries(scale = "medium", returnclass = "sf")
    world$total_flow <- sapply(world$iso_a2, function(x){
        tmp <- mov_mat_summary %>% filter((code_o==x) | (code_d==x)) %>% .$total_flight_pair
        sum(tmp)
    })
    mov_mat_summary$line_weight <- log10(mov_mat_summary$total_flight_pair)
    mov_mat_summary$line_weight <- mov_mat_summary$line_weight/max(mov_mat_summary$line_weight)*1

    # curve_points <- gcIntermediate(p1 = mov_mat_summary %>% select(lon_o, lat_o) %>% as.matrix(), p2 = mov_mat_summary %>% select(lon_d, lat_d) %>% as.matrix(), n = 20, addStartEnd = TRUE, sp = TRUE)
    # # Convert SpatialPoints to sf object
    # curve_points_sf <- st_as_sf(curve_points)

    p <- ggplot() + 
        geom_sf(data = world, aes(fill = log10(total_flow)), lwd = 0, colour = "white") + 
        scale_fill_gradientn(name="Total flow in each country\n(log 10 scale)", colors = c("#9DBF9E", "#FCB97D", "#A84268"),na.value = "grey80")+
        geom_curve(data = mov_mat_summary, aes(x = lon_o, y = lat_o, xend = lon_d, yend = lat_d, linewidth=I(line_weight), color=data_source), alpha=0.1)+
        geom_text_repel(data=df_pos, aes(x = lon, y = lat, label = code_2), col = "black", size = 2, segment.color = NA) + 
        # geom_sf(data = curve_points_sf, color="black", size = 0.1, alpha=0.5) +
        scale_colour_manual(name="Data source", values=c("dark red", "dark blue", "black"))+
        scale_size_identity()+
        theme_void() +
        theme(
            # legend.justification defines the edge of the legend that the legend.position coordinates refer to
            legend.justification = c(0, 1),
            # Set the legend flush with the left side of the plot, and just slightly below the top of the plot
            legend.position = c(0.05, .65)
            )+
        guides(colour = guide_legend(override.aes = list(alpha = 1)))
    p
}
