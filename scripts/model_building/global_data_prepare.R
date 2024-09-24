library(spatPomp)
library(pomp)
library(lubridate)
library(tidyverse)
library(parallel)
library(slider)

global_data_prepare <- function(
  date_start,
  date_end,
  U = NA, # include top U countries 
  U_each_continent = NA, # include top U countries in each continent
  top_by = NA, # "cases", "flow", sort the top U countries by cases or flow
  threshold = 0, # minimum number of cases/flow to be included
  threshold_country = "ZZZZZZ", # the included countries must have cases/flow higher than that in this country
  global_scale = "country", # "continent", "country", "country_state"
  lineage_top_n = 3,
  add_HK = FALSE,
  check_plot_airline = FALSE
  ) {
  # 0. Sanity check #
  date_start <- ymd(date_start)
  date_start_10d_before <- date_start-10
  stopifnot(date_start_10d_before > ymd("2019-12-31"))
  date_end <- ymd(date_end)
  date_start_d <- decimal_date(date_start)
  date_end_d <- decimal_date(date_end)
  stopifnot(
    (date_start_d > 2019) &&
    (date_start_d < 2023) &&
    (date_end_d > 2019) &&
    (date_end_d < 2023) &&
    ((date_end_d - date_start_d) >= 1/12) # at least one month
  )

  stopifnot(global_scale %in% c("continent", "country", "country_state"))
  if (global_scale == "country") {
    this_level <- 3
  } else if (global_scale == "continent") {
    this_level <- 2
  }

  stopifnot(lineage_top_n<=10)

  stopifnot(!is.na(U) | !is.na(U_each_continent))
  stopifnot(top_by %in% c("cases", "flow", "cases and flow"))

  # 1. Loading appropriate data #
  cross_check_table <- readxl::read_excel(
    here::here(
      "data//our_airports_data/cross_check_table_ihme_input_completed.xlsx"
    )
  )
  cross_check_table <- cross_check_table %>% filter(level==this_level)

  # IHME data # 
  df_ihme <- readRDS(here::here("data/IHME/df_ihme.rds"))
  df_ihme <- df_ihme %>% filter(level==this_level) %>% filter(date>=ymd(date_start_10d_before) & date<=ymd(date_end))
  df_ihme <- df_ihme %>% mutate(date_decimal=decimal_date(date))
  df_ihme <- left_join(df_ihme, cross_check_table %>% select(-level), by="code") 

  # travel policy data
  df_travel_policy <- read_csv("./data/covid-policy-tracker/international-travel-covid.csv")
  names(df_travel_policy)[c(2:4)] <- c("country_iso3c", "date", "travel_control_level")
  df_ihme <- left_join(df_ihme, df_travel_policy, by=c("date", "country_iso3c"))

  # load movement matrix #
  year_start <- year(ymd(date_start))
  year_end <- year(ymd(date_end))

  mov_mat <- lapply(unique(c(year_start, year_end)), function(this_year) {
    tmp <- readRDS(here::here(paste0(
      "data/estimated_movement_matrix/move_mat_",
      this_year,
      "_",
      global_scale,
      ".rds")))
    tmp
  })
  mov_mat <- bind_rows(mov_mat)
  mov_mat <- mov_mat %>% filter(date>=ymd(date_start) & date<=ymd(date_end))
  units_travel <- sort(unique(c(mov_mat$code_o, mov_mat$code_d)))
  stopifnot(all(units_travel %in% df_ihme$code))
  df_ihme <- df_ihme %>% filter(code %in% units_travel)

  # re-determine the countries under investigation according to input U #
  get_new_code <- function(top_by){
    if(top_by == "cases"){ # top_by == "cases"
      if(is.na(U)){
        if(is.na(U_each_continent)){ # include all units
          units_cases <- df_ihme %>% group_by(code) %>% summarise(total_case=sum(daily_cases)) %>% arrange(desc(total_case)) %>% filter(total_case>=threshold) %>% .$code
          if(threshold_country %in% units_cases){
            units_cases <- units_cases[1:which(units_cases==threshold_country)]}
          if(add_HK){units_cases <- unique(c(units_cases, "HK"))}
          units_cases <- units_cases[units_cases %in% units_travel]
          df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code=code) %>% filter(code %in% units_flow)
        } else { # sort by U_each_continent
          units_cases <- df_ihme %>% group_by(continent, code) %>% summarise(total_case=sum(daily_cases)) %>% ungroup() %>% group_by(continent) %>% arrange(desc(total_case)) %>% filter(row_number()<=U_each_continent)%>% filter(total_case>=threshold) %>% .$code
          if(threshold_country %in% units_cases){
            units_cases <- units_cases[1:which(units_cases==threshold_country)]}
          if(add_HK){units_cases <- unique(c(units_cases, "HK"))}
          units_cases <- units_cases[units_cases %in% units_travel]
          df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code=paste0(continent, "_others"))
          df_new_code$new_code[df_new_code$code %in% units_cases] <- df_new_code$code[df_new_code$code %in% units_cases]
        }
      } else { # sort by U
        units_cases <- df_ihme %>% group_by(code) %>% summarise(total_case=sum(daily_cases)) %>% arrange(desc(total_case)) %>% filter(row_number()<=U) %>% filter(total_case>=threshold) %>% .$code
        if(threshold_country %in% units_cases){
            units_cases <- units_cases[1:which(units_cases==threshold_country)]}
        if(add_HK){units_cases <- unique(c(units_cases, "HK"))}
        units_cases <- units_cases[units_cases %in% units_travel]
        df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code="Others")
        df_new_code$new_code[df_new_code$code %in% units_cases] <- df_new_code$code[df_new_code$code %in% units_cases]
      }

    } else if (top_by == "flow") { # top_by == "flow"
      mov_mat_total_flow <- left_join(mov_mat %>% mutate(code=code_o) %>% group_by(code) %>% summarise(flow_all_out=sum(flow_all)) %>% ungroup(), mov_mat %>% mutate(code=code_d) %>% group_by(code) %>% summarise(flow_all_in=sum(flow_all)) %>% ungroup()) %>% mutate(flow_all=flow_all_out+flow_all_in) %>% select(code, flow_all) %>% arrange(desc(flow_all))
      mov_mat_total_flow <- left_join(mov_mat_total_flow, cross_check_table)
      mov_mat_total_flow %>% filter(continent=="AS")
      if(is.na(U)){
        if(is.na(U_each_continent)){ # include all units
          units_flow <- mov_mat_total_flow %>% arrange(desc(flow_all)) %>% ungroup() %>% filter(flow_all>=threshold) %>% .$code
          if(threshold_country %in% units_flow){
            units_flow <- units_flow[1:which(units_flow==threshold_country)]}
          if(add_HK){units_flow <- unique(c(units_flow, "HK"))}
          df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code=code) %>% filter(code %in% units_flow)
        } else { # sort by U_each_continent
          # mov_mat_total_flow %>% filter(country_iso2c=="ZA")
          units_flow <- mov_mat_total_flow %>% group_by(continent) %>% arrange(desc(flow_all)) %>% filter(row_number()<=U_each_continent)%>% filter(flow_all>=threshold) %>% ungroup() %>% .$code
          if(threshold_country %in% units_flow){
            units_flow <- units_flow[1:which(units_flow==threshold_country)]}
          if(add_HK){units_flow <- unique(c(units_flow, "HK"))}
          df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code=paste0(continent, "_others"))
          df_new_code$new_code[df_new_code$code %in% units_flow] <- df_new_code$code[df_new_code$code %in% units_flow]
        }
      } else { # sort by U
        units_flow <- mov_mat_total_flow %>% arrange(desc(flow_all)) %>% filter(row_number()<=U)%>% filter(flow_all>=threshold) %>% ungroup() %>% .$code
        if(threshold_country %in% units_flow){
            units_flow <- units_flow[1:which(units_flow==threshold_country)]}
        if(add_HK){units_flow <- unique(c(units_flow, "HK"))}
        df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code="Others")
        df_new_code$new_code[df_new_code$code %in% units_flow] <- df_new_code$code[df_new_code$code %in% units_flow]
      }
    }
    return(df_new_code)
  }

  if(top_by == "cases") { # cases
    df_new_code <- get_new_code("cases")
  } else if (top_by == "flow") { # flow
    df_new_code <- get_new_code("flow")
  } else if (top_by == "cases and flow") { # cases and flow
    df_new_code_cases <- get_new_code("cases")
    df_new_code_flow <- get_new_code("flow")
    units_selected <- union(unique(df_new_code_flow$new_code),unique(df_new_code_cases$new_code))
    df_new_code <- df_ihme %>% select(continent, code) %>% unique() %>% mutate(new_code=paste0(continent, "_others"))
    df_new_code$new_code[df_new_code$code %in% units_selected] <- df_new_code$code[df_new_code$code %in% units_selected]
  }

  country_under_investigation <- sort(unique(df_new_code$new_code))
  num_of_country_under_investigation <- length(country_under_investigation)
  # cross_check_table %>% filter(code %in%country_under_investigation) %>% arrange(continent,loc_name)

  # assign new code to df_ihme #
  df_ihme <- left_join(df_ihme, df_new_code, by="code")
  ## moving average of daily_cases and daily_deaths
  window_size <- 7
  df_ihme <- df_ihme %>% mutate(code=new_code) %>% group_by(code) %>% mutate(daily_cases = map_dbl(seq_along(daily_cases), ~mean(daily_cases[max(1, .x - window_size + 1):.x])), daily_deaths = map_dbl(seq_along(daily_deaths), ~mean(daily_deaths[max(1, .x - window_size + 1):.x]))) %>% ungroup()
  ## aggregate by date and code
  df_ihme_aggregate <- df_ihme %>% ungroup() %>% 
    select(date, date_decimal, code, population, mandates_mean, infection_detection, infection_fatality, cumulative_all_effectively_vaccinated, cumulative_all_fully_vaccinated, inf_cuml_mean, inf_cuml_mean_unvax, inf_cuml_mean_vax, daily_cases, daily_deaths, inf_mean, travel_control_level) %>%
    group_by(date, date_decimal, code) %>%
    summarise(
      travel_control_level=round(mean(travel_control_level, na.rm = T), 0),
      mandates_mean=sum(mandates_mean*population, na.rm = T)/sum(population, na.rm = T),
      infection_detection=sum(daily_cases+daily_deaths, na.rm = T)/sum(inf_mean, na.rm = T),
      infection_fatality=sum(daily_deaths, na.rm = T)/sum(inf_mean, na.rm = T),
      cumulative_all_effectively_vaccinated=sum(cumulative_all_effectively_vaccinated, na.rm = T),
      cumulative_all_fully_vaccinated=sum(cumulative_all_fully_vaccinated, na.rm = T),
      inf_cuml_mean=sum(inf_cuml_mean, na.rm = T),
      inf_cuml_mean_unvax=sum(inf_cuml_mean_unvax, na.rm = T),
      inf_cuml_mean_vax=sum(inf_cuml_mean_vax, na.rm = T),
      daily_cases=sum(daily_cases, na.rm = T),
      daily_deaths=sum(daily_deaths, na.rm = T),
      population=sum(population, na.rm = T)) %>% 
    ungroup()
  df_ihme_aggregate <- df_ihme_aggregate %>% mutate(code=factor(code, levels = country_under_investigation)) %>% arrange(code, date_decimal)

  df_ihme_aggregate$infection_fatality[df_ihme_aggregate$infection_fatality>1] <- 1
  df_ihme_aggregate$infection_detection[df_ihme_aggregate$infection_detection>1] <- 1

  # assign new code to mov_mat #
  summary_source <- function(x){
    if(length(unique(x))==1){return(x[1])}else{return("Mixed")}
  }
  mov_mat_aggregate <- mov_mat %>% left_join(df_new_code %>% transmute(code_o=code, new_code_o=new_code), by= "code_o") %>% left_join(df_new_code %>% transmute(code_d=code, new_code_d=new_code), by= "code_d") %>% mutate(code_o=new_code_o, code_d=new_code_d) %>% select(date, code_o, code_d, flow_all, data_source) %>% group_by(date, code_o, code_d) %>% summarise(flow_all=sum(flow_all), data_source=summary_source(data_source)) %>% ungroup() %>% arrange(code_o, code_d, date) %>% filter(code_o!=code_d)

  ## moving average of mov_mat
  window_size=7
  mov_mat_aggregate <- cross_join(tibble(date=seq.Date(date_start, date_end, by="day")), tibble(code_o=country_under_investigation)) %>% cross_join(tibble(code_d=country_under_investigation)) %>% filter(code_o!=code_d) %>% left_join(mov_mat_aggregate, by=c("date", "code_o", "code_d")) %>% mutate(flow_all=ifelse(is.na(flow_all), 0, flow_all)) %>% arrange(code_o, code_d, date) 

  mov_mat_aggregate <- mov_mat_aggregate %>% group_by(code_o, code_d) %>% mutate(flow_all = map_dbl(seq_along(flow_all), ~mean(flow_all[max(1, .x - window_size + 1):.x]))) %>% ungroup()

  if(check_plot_airline){
    source(here::here("scripts/data_processing/helper/plot_airline.R"))
    plot_airline <- plot_global_country(mov_mat_aggregate %>% filter(!is.na(data_source)))
    ggsave(here::here("results/figs/flights_map_Omicron20.jpg"), plot=plot_airline, width=12, height=12*0.618, dpi = 450)
  }  
  
  # load GISAID data #
  source(here::here("scripts/data_processing/helper/load_GISAID.R"))
  list_gisaid <- load_GISAID(date_start=date_start_10d_before, date_end=date_end, country_under_investigation, df_new_code)
  if(length(list_gisaid$lineages)>lineage_top_n+1){
    warning("There are more than specified number of lineages in the study period. Please adjust manually before fitting model.")
  }

  # 2. Build essential covariate matrix #
  # Build movement matrix Csnippet #
  U <- num_of_country_under_investigation
  n_days <- as.numeric(date_end-date_start)
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  v_by_g_C_arrays <- mclapply(1:n_days, function(this_day){ # Dimension 1: day
    this_date <- date_start+this_day-1
    # print(this_date)
    mov_mat_day <- mov_mat_aggregate %>% filter(date==this_date)
    v_by_g <- matrix(0,U,U)
    for(u1 in 1:U){ # Dimension 2: from
      for(u2 in 1:U){ # Dimension 3: to
        if(u1==u2){next()}
        flow <- mov_mat_day$flow_all[(mov_mat_day$code_o==country_under_investigation[u1]) & (mov_mat_day$code_d==country_under_investigation[u2])]
        flow <- flow*365 # convert to annual flow
        if(length(flow)==0){
          # print(paste0("no flow for ", country_under_investigation[u1], " to ", country_under_investigation[u2], " on ", this_date))
          v_by_g[u1,u2] <- 0
        }else{
          v_by_g[u1,u2] <- flow
        }
      }
    }
    v_by_g_C_rows <- apply(v_by_g,1,to_C_array) # Dimension 2: from
    v_by_g_C_array <- to_C_array(v_by_g_C_rows) # Dimension 3: to
  }, mc.cores=8)
  v_by_g_C_arrays <- to_C_array(unlist(v_by_g_C_arrays))
  v_by_g_C <- Csnippet(paste0("const double v_by_g[", n_days, "][", U, "][", U, "] = ", v_by_g_C_arrays,"; "))

  # 3. Specify the POMP model #
  # 3.1. Getting a measurements data.frame with columns for times, spatial units and measurements.
  df_meas <- df_ihme_aggregate %>% select(date, date_decimal, code, daily_cases, daily_deaths)
  if(any(df_meas$daily_cases<0, na.rm = TRUE)){
    df_meas$daily_cases[df_meas$daily_cases<0] <- 0
  }
  if(any(df_meas$daily_deaths<0, na.rm = TRUE)){
    df_meas$daily_deaths[df_meas$daily_deaths<0] <- 0
  }

  # 3.2. Getting a covariates data.frame with columns for times, spatial units and covariate data.
  df_covar <- df_ihme_aggregate %>% select(date, date_decimal, code,
    travel_control_level, population, mandates_mean, infection_detection, infection_fatality, cumulative_all_effectively_vaccinated, cumulative_all_fully_vaccinated, inf_cuml_mean, inf_cuml_mean_unvax)

  # 3.3 get the pre-estimated $q_u$
  df_transit <- readxl::read_excel(here::here("data/transit_data/cross_check_table_transit.xlsx"))
  df_transit <- df_transit %>% filter(level==this_level) %>% select(code, transfer_rates)
  df_transit <- left_join(df_transit, df_new_code, by="code")
  df_transit <- df_transit %>% filter(!is.na(new_code))
  df_transit <- df_transit %>% group_by(new_code) %>% summarise(transfer_rates=mean(transfer_rates, na.rm = T))
  transit_rates <- df_transit$transfer_rates[match(df_transit$new_code, country_under_investigation)]
  stopifnot(length(transit_rates)==length(country_under_investigation))

  # return list of global data #
  return(list(
    dates = list(date_start=date_start, date_end=date_end, date_start_d=date_start_d, date_end_d=date_end_d),
    df_new_code=df_new_code,
    df_meas = df_meas,
    df_covar = df_covar,
    country_under_investigation = country_under_investigation,
    num_of_country_under_investigation = num_of_country_under_investigation,
    travel = list(n_days=n_days, order=country_under_investigation, v_by_g_C=v_by_g_C, mov_mat=mov_mat_aggregate),
    transit_rates = transit_rates,
    # plot_airline = plot_airline,
    gisaid = list_gisaid,
    level = this_level))
}
