get_params_for_travel_hubs <- function(
  strategy_i, 
  traveler_weight_i,
  base_model_name,
  data_fitting,
  travel_volume_type = "per_capita",
  travel_data_period = "pandemic"
){
  if(strategy_i == "Max"){ # no resources reallocation
    return(
      tibble(
        code = data_fitting$country_under_investigation,
        change_DSR = 1,
        change_IDR = 1,
        traveler_weight = traveler_weight_i,
        travel_hubs_selected = list(data_fitting$country_under_investigation)
      )
    )
  }
  stopifnot(traveler_weight_i>=0)
  stopifnot(traveler_weight_i<=1)

  strategy_split_i <- strsplit(strategy_i, "-") %>% unlist()
  if(grepl("+", strategy_i, fixed = TRUE)){
    check_adding_others <- TRUE
    rank_to_add <- as.numeric(gsub("^T\\d+\\+", "", strategy_split_i[1]))
    top_n <- as.numeric(gsub("\\D", "", gsub("\\+\\d+", "", strategy_split_i[1])))
    check_top_n_add <- grepl("a", strategy_split_i[1])
    reallocation_type <- strategy_split_i[2]
  } else {
    check_adding_others <- FALSE
    top_n <- as.numeric(gsub("\\D", "", strategy_split_i[1]))
    check_top_n_add <- grepl("a", strategy_split_i[1])
    reallocation_type <- strategy_split_i[2]
  }

  stopifnot(reallocation_type %in% c("R", "H", "M")) # R: regular, H: half, M: more
  if(top_n > data_fitting$num_of_country_under_investigation){
    stop("The number of top travel hubs should be less than the number of countries.")
  }
  stopifnot(travel_volume_type %in% c("per_capita", "total"))
  stopifnot(travel_data_period %in% c("pandemic", "normal"))

  # read data
  data_measurements <- readRDS(paste0("results/model_data/data_measurements_", base_model_name, ".rds"))

  # list all the travel hubs in order
  source("scripts/data_processing/travel_data_2019.R")
  if(travel_data_period == "pandemic"){
    df_travel_hubs <- data_fitting$travel$mov_mat %>% group_by(code_d) %>% summarise(total=sum(flow_all)) %>% arrange(desc(total)) %>% transmute(code=code_d, total=total)
  } else if(travel_data_period == "normal"){
    df_travel_hubs <- get_mov_mat_aggregate_2019() %>% group_by(code_d) %>% summarise(total=sum(flow_all)) %>% arrange(desc(total)) %>% transmute(code=code_d, total=total*365)
  }

  df_travel_hubs <- data_fitting$df_covar %>% select(code, population) %>% unique() %>% left_join(df_travel_hubs, by="code")
  df_travel_hubs$travel_per_capita <- df_travel_hubs$total / df_travel_hubs$population
  df_travel_hubs <- df_travel_hubs %>% left_join(data_fitting$df_new_code %>% transmute(code=new_code, continent) %>% unique(), by = "code")

  if(travel_volume_type == "total"){
    df_travel_hubs <- df_travel_hubs %>% mutate(relative_weight_global=total/sum(total)) %>% arrange(grepl("others", code), desc(total)) %>% mutate(travel_hub_rank_global = row_number()) %>% group_by(continent) %>% arrange(grepl("others", code), desc(total)) %>% mutate(travel_hub_rank_in_continent=row_number()) %>% ungroup() # make sure code "others" is always at the bottom, then sort by total
    write_csv(df_travel_hubs, paste0("results/model_data/data_travel_hubs_order_total_", travel_data_period, ".csv"))
    travel_hubs_selected <- df_travel_hubs %>% filter(travel_hub_rank_global <= top_n) %>% pull(code)
    if(check_top_n_add){
      travel_hubs_selected <- unique(c(travel_hubs_selected, df_travel_hubs %>% filter(travel_hub_rank_in_continent==1) %>% pull(code)))
    }
    if(check_adding_others){
      travel_hubs_selected <- unique(c(travel_hubs_selected, df_travel_hubs %>% filter(travel_hub_rank_global==(rank_to_add+top_n)) %>% pull(code)))
    }
  } else if(travel_volume_type == "per_capita"){
    df_travel_hubs <- df_travel_hubs %>% mutate(relative_weight_global=travel_per_capita/sum(travel_per_capita)) %>% arrange(grepl("others", code), desc(travel_per_capita)) %>% mutate(travel_hub_rank_global = row_number()) %>% group_by(continent) %>% arrange(grepl("others", code), desc(travel_per_capita)) %>% mutate(travel_hub_rank_in_continent=row_number()) %>% ungroup() # make sure code "others" is always at the bottom, then sort by travel_per_capita
    write_csv(df_travel_hubs, paste0("results/model_data/data_travel_hubs_order_percapita_", travel_data_period, ".csv"))
    travel_hubs_selected <- df_travel_hubs %>% filter(travel_hub_rank_global <= top_n) %>% pull(code)
    if(check_top_n_add){
      travel_hubs_selected <- unique(c(travel_hubs_selected, df_travel_hubs %>% filter(travel_hub_rank_in_continent==1) %>% pull(code)))
    }
    if(check_adding_others){
      travel_hubs_selected <- unique(c(travel_hubs_selected, df_travel_hubs %>% filter(travel_hub_rank_global==(rank_to_add+top_n)) %>% pull(code)))
    }
  }

  # calculate the total sequencing amount for each region
  df_sequencing_amount_total <- data_fitting$gisaid$data_GISAID %>% group_by(code) %>% summarise(total_N=sum(N)) %>% mutate(check_in_selected_hubs=code %in% travel_hubs_selected) %>% left_join(data_fitting$df_new_code %>% transmute(code=new_code, continent) %>% unique(), by="code") %>% arrange(code)

  df_diagnostic_amount_total <- data_measurements %>% transmute(code, total_diag = daily_cases + daily_deaths) %>% group_by(code) %>% summarise(total_diag=sum(total_diag))

  if(reallocation_type == "H"){
    prop_reallocation <- 0.5
  } else if(reallocation_type == "M"){
    prop_reallocation <- 0.9
  }

  # calculate the change_DSR for the top travel hubs
  if(reallocation_type == "R"){ # regular, no reallocation of resources
    change_DSR_travel_hubs <- 1
    change_IDR_travel_hubs <- 1
    df_out <- tibble(code=data_fitting$country_under_investigation) %>% 
      mutate(change_DSR=ifelse(code %in% travel_hubs_selected, change_DSR_travel_hubs, 1)) %>% 
      mutate(change_IDR=ifelse(code %in% travel_hubs_selected, change_IDR_travel_hubs, 1)) %>% 
      mutate(traveler_weight=ifelse(code %in% travel_hubs_selected, traveler_weight_i, 999))
  } else { # half of the resources (not in selected hubs) are reallocated
    ## Sequencing
    df_sequencing_amount_add <- df_sequencing_amount_total %>% 
      mutate(total_N_not_selected_half=sum(total_N[!check_in_selected_hubs])*prop_reallocation) %>% 
      ungroup() %>% 
      left_join(df_travel_hubs %>% select(code, relative_weight_global, travel_hub_rank_in_continent), by="code") %>% 
      filter(code %in% travel_hubs_selected) %>% 
      mutate(new_weight=relative_weight_global/sum(relative_weight_global)) %>% 
      ungroup() %>% 
      mutate(new_N = total_N_not_selected_half * new_weight) %>% 
      select(code, new_N) # reallocate half of the resources in the non-selected hubs to the selected hubs by the relative weight (per capita or total depending on the travel_volume_type)
    df_sequencing_amount_total_new <- left_join(df_sequencing_amount_total, df_sequencing_amount_add, by = "code") %>% 
      mutate(total_N_new=ifelse(check_in_selected_hubs, total_N+new_N, total_N*(1-prop_reallocation))) %>% 
      select(code, total_N, total_N_new) 
    stopifnot(abs(sum(df_sequencing_amount_total_new$total_N) - sum(df_sequencing_amount_total_new$total_N_new)) < 1) # make sure the total amount is the same
    change_DSR_codes <- df_sequencing_amount_total_new$total_N_new / df_sequencing_amount_total_new$total_N
    change_DSR_codes[df_sequencing_amount_total_new$total_N==0] <- 1

    ## Diagnostics
    df_diag_amount_total <- df_diagnostic_amount_total %>% 
      mutate(check_in_selected_hubs=code %in% travel_hubs_selected) %>% 
      left_join(data_fitting$df_new_code %>% 
      transmute(code=new_code, continent) %>% unique(), by="code") %>% 
      arrange(code)
    df_diag_amount_add <- df_diag_amount_total %>% 
      mutate(total_diag_not_selected_half=sum(total_diag[!check_in_selected_hubs])*prop_reallocation) %>% 
      ungroup() %>% 
      left_join(df_travel_hubs %>% select(code, relative_weight_global, travel_hub_rank_in_continent), by="code") %>% 
      filter(code %in% travel_hubs_selected) %>%
      mutate(new_weight=relative_weight_global/sum(relative_weight_global)) %>% 
      ungroup() %>% 
      mutate(new_diag = total_diag_not_selected_half * new_weight) %>% 
      select(code, new_diag) # reallocate half of the resources in the non-selected hubs to the selected hubs by the relative weight (per capita or total depending on the travel_volume_type)
    df_diag_amount_total_new <- left_join(df_diag_amount_total, df_diag_amount_add, by = "code") %>% 
      mutate(total_diag_new=ifelse(check_in_selected_hubs, total_diag+new_diag, total_diag*(1-prop_reallocation))) %>% 
      select(code, total_diag, total_diag_new)
    stopifnot(abs(sum(df_diag_amount_total_new$total_diag) - sum(df_diag_amount_total_new$total_diag_new)) < 1) # make sure the total amount is the same
    change_IDR_codes <- df_diag_amount_total_new$total_diag_new / df_diag_amount_total_new$total_diag
    change_IDR_codes[df_diag_amount_total_new$total_diag==0] <- 1

    df_out <- tibble(code=data_fitting$country_under_investigation) %>%
      mutate(change_DSR=change_DSR_codes[match(code, df_sequencing_amount_total_new$code)]) %>% 
      mutate(change_IDR=change_IDR_codes[match(code, df_diag_amount_total_new$code)]) %>% 
      mutate(traveler_weight=ifelse(code %in% travel_hubs_selected, traveler_weight_i, 999))
  } 

  df_out$travel_hubs_selected <- list(travel_hubs_selected)
  return(df_out)
}