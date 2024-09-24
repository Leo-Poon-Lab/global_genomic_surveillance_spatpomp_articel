get_params_for_travel_hubs <- function(
  top_n = NA, 
  base_model_name,
  data_fitting,
  travel_volume_type = "per_capita"
){
  if(is.na(top_n)){
    return(
      tibble(
        code = data_fitting$country_under_investigation,
        change_DSR = NA,
        change_IDR = NA
      )
    )
  }
  if(top_n > data_fitting$num_of_country_under_investigation){
    stop("The number of top travel hubs should be less than the number of countries.")
  }
  stopifnot(travel_volume_type %in% c("per_capita", "total"))

  # read data
  data_measurements <- readRDS(paste0("results/model_data/data_measurements_", base_model_name, ".rds"))

  # list all the travel hubs in order
  df_travel_hubs <- data_fitting$travel$mov_mat %>% group_by(code_d) %>% summarise(total=sum(flow_all)) %>% arrange(desc(total))
  df_travel_hubs <- data_fitting$df_covar %>% select(code, population) %>% unique() %>% mutate(code_d=code) %>% left_join(df_travel_hubs, by="code_d")
  df_travel_hubs$travel_per_capita <- df_travel_hubs$total / df_travel_hubs$population
  # write_csv(df_travel_hubs, "results/model_data/data_travel_hubs_order.csv")
  if(travel_volume_type == "total"){
    travel_hubs_selected <- df_travel_hubs %>% arrange(desc(total)) %>% head(top_n) %>% pull(code_d)
  } else{
    travel_hubs_selected <- df_travel_hubs %>% arrange(desc(travel_per_capita)) %>% head(top_n) %>% pull(code_d)
  }

  # calculate the total sequencing amount globally
  sequencing_amount_total <- sum(data_fitting$gisaid$data_GISAID$N)

  # calculate the sequencing amount of the top travel hubs
  sequencing_amount_travel_hubs <- sum(data_fitting$gisaid$data_GISAID %>% filter(code %in% travel_hubs_selected) %>% pull(N))

  # calculate the change_DSR for the top travel hubs
  change_DSR_travel_hubs <- sequencing_amount_total / sequencing_amount_travel_hubs

  # calculate the total detection amount globally
  detection_amount_total <- sum(data_measurements$daily_cases + data_measurements$daily_deaths)

  # calculate the detection amount of the top travel hubs
  detection_amount_travel_hubs <- sum(data_measurements$daily_cases[data_measurements$code %in% travel_hubs_selected] + data_measurements$daily_deaths[data_measurements$code %in% travel_hubs_selected])

  # calculate the change_IDR for the top travel hubs
  change_IDR_travel_hubs <- detection_amount_total / detection_amount_travel_hubs

  # return results
  ## make sure the order are the same as the data_fitting$country_under_investigation
  df_out <- left_join(tibble(code=data_fitting$country_under_investigation), tibble(code=travel_hubs_selected), by="code") %>% mutate(change_DSR=ifelse(code %in% travel_hubs_selected, change_DSR_travel_hubs, 0)) %>% mutate(change_IDR=ifelse(code %in% travel_hubs_selected, change_IDR_travel_hubs, 0))

  return(df_out)
}