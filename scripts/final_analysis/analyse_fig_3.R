library(tidyverse)

# Fig 3

## Fig 3a
## Where are the Omicron BA.1 and BA.2 variants first identified?
### For IDR/DSR=1, mainly (>99%) in South Africa

## Baseline conditions
(df_fig_3a_base_hline <- readxl::read_excel("results/figs/model_simulation/M3/fig_3a_base_hline.xlsx"))

## Fig 3c
## Why sequencing more community cases works better in BA.1 than BA.2?
### Check the IDR and DSR
data_covar <- readRDS("results/model_data/data_covariates_Omicron20.rds")

lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
df_vline <- tibble(lag=lags_intro)
df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(data_covar$date_decimal)
df_vline$date <- date_decimal(df_vline$date_decimal)

### Check the number of BA.1 / BA.2 cases in different periods
model_sims <- readRDS("results/figs/model_simulation/Omicron20/model_sims.rds")
model_sims_median <- model_sims %>% filter(unitname=="ZA") %>% transmute(
  date = as_date(date_decimal(date_decimal)),
  tot_cases = daily_cases,
  C_k_c_i = C_1_c_infected_new + C_2_c_infected_new + C_3_c_infected_new + C_4_c_infected_new,
  C_2_c_i = C_2_c_infected_new,
  C_3_c_i = C_3_c_infected_new,
  C_k_c_d = C_1_c_detected_new + C_2_c_detected_new + C_3_c_detected_new + C_4_c_detected_new,
  C_2_c_d = C_2_c_detected_new,
  C_3_c_d = C_3_c_detected_new,
  C_K_c_s = C_1_c_sequenced_new + C_2_c_sequenced_new + C_3_c_sequenced_new + C_4_c_sequenced_new,
  C_2_c_s = C_2_c_sequenced_new,
  C_3_c_s = C_3_c_sequenced_new
  ) %>% group_by(date) %>% summarise_all(median)

(date_BA1_detected <- model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+20) %>% filter(C_2_c_s>=1) %>% slice_head() %>% .$date)
model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+20) %>% transmute(date, C_2_c_i_prop=C_2_c_i/C_k_c_i, C_3_c_i_prop=C_3_c_i/C_k_c_i)%>% filter(date == date_BA1_detected) # propotion of community infections in total infections, community prevalence
data_covar %>% filter(code=="ZA") %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date, infection_detection, detection_sequencing) %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+14) %>% filter(date == date_BA1_detected) # IDR and DSR when the first BA.1 case is detected

(date_BA2_detected <- model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+20) %>% filter(C_3_c_s>=1) %>% slice_head() %>% .$date)
model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+20) %>% transmute(date, C_2_c_i_prop=C_2_c_i/C_k_c_i, C_3_c_i_prop=C_3_c_i/C_k_c_i) %>% filter(date == date_BA2_detected) # propotion of community infections in total infections, community prevalence
data_covar %>% filter(code=="ZA") %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date, infection_detection, detection_sequencing) %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+20) %>% filter(date == date_BA2_detected) # IDR and DSR when the first BA.2 case is detected

## analyse the diff in identification lags of the preferred strategy
df_baseline <- readxl::read_excel("results/figs/model_simulation/M3/fig_3a_base_hline.xlsx")

df_3c_BA2_hm <- read_csv("results/figs/model_simulation/M3/fig_3c_Omicron_BA.2_heatmap_source.csv")
df_3c_BA2_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 9 day max

df_3d_BA2_hm <- read_csv("results/figs/model_simulation/M3/fig_3d_Omicron_BA.2_heatmap_source.csv")
df_3d_BA2_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 9.5 day max

df_3e_BA1_hm <- read_csv("results/figs/model_simulation/M3/fig_3e_Omicron_BA.1_heatmap_source.csv")
df_3e_BA1_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 1 day

df_3e_BA2_hm <- read_csv("results/figs/model_simulation/M3/fig_3e_Omicron_BA.2_heatmap_source.csv")
df_3e_BA2_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 9 day

df_3f_BA1_hm <- read_csv("results/figs/model_simulation/M3/fig_3f_Omicron_BA.1_heatmap_source.csv")
df_3f_BA1_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 1 day

df_3f_BA2_hm <- read_csv("results/figs/model_simulation/M3/fig_3f_Omicron_BA.2_heatmap_source.csv")
df_3f_BA2_hm %>% left_join(df_baseline %>% filter(case_type=="Community") %>% transmute(variant, lag_median_base=lag_median)) %>% mutate(lag_diff=lag_median-lag_median_base) %>% arrange(lag_diff) %>% select(variant, parval_1, parval_2, lag_median, lag_median_base, lag_diff) %>% slice_min(lag_diff) # shortened for 9 day

## Table S2
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)
data_fitting <- readRDS(paste0("results/model_data/data_fitting_", "Omicron20", ".rds"))
df_travel_order_t_nor <- read_csv("results/model_data/data_travel_hubs_order_total_normal.csv")
df_travel_order_t_pan <- read_csv("results/model_data/data_travel_hubs_order_total_pandemic.csv")
df_travel_order_p_nor <- read_csv("results/model_data/data_travel_hubs_order_percapita_normal.csv")
df_travel_order_p_pan <- read_csv("results/model_data/data_travel_hubs_order_percapita_pandemic.csv")

df_s2 <- df_travel_order_t_nor %>% transmute(code, population_m = population/1000000, total_normal_m=total/1000000, travel_hub_rank_total_normal = paste0(travel_hub_rank_global, ifelse(travel_hub_rank_in_continent==1, "*", "")))
df_s2 <- left_join(df_s2, df_travel_order_t_pan %>% transmute(code, total_pandemic_m=total/1000000, travel_hub_rank_total_pandemic = paste0(travel_hub_rank_global, ifelse(travel_hub_rank_in_continent==1, "*", ""))), by="code")
df_s2 <- left_join(df_s2, df_travel_order_p_nor %>% transmute(code, travel_hub_rank_percapita_normal = paste0(travel_hub_rank_global, ifelse(travel_hub_rank_in_continent==1, "*", ""))), by="code")
df_s2 <- left_join(df_s2, df_travel_order_p_pan %>% transmute(code, travel_hub_rank_percapita_pandemic = paste0(travel_hub_rank_global, ifelse(travel_hub_rank_in_continent==1, "*", ""))), by="code")

df_s2 <- left_join(df_s2, cross_check_table %>% transmute(loc_name, code, continent), by="code") 
df_s2 <- left_join(df_s2, tibble(code=data_fitting$country_under_investigation, transit_rates= data_fitting$transit_rates), by="code")

data_covar_m4 <- readRDS("results/model_data/data_covariates_M4.rds")
df_s2 <- data_covar_m4 %>% transmute(code, IDR=infection_detection, DSR=detection_sequencing, diagnostic_capacity_daily=diagnostic_capacity, sequencing_capacity_daily=sequencing_capacity) %>% unique() %>% left_join(df_s2)

df_s2 <- df_s2 %>% select(loc_name, code, continent, population_m, total_pandemic_m, total_normal_m, travel_hub_rank_total_pandemic, travel_hub_rank_total_normal, travel_hub_rank_percapita_pandemic, travel_hub_rank_percapita_normal, transit_rates, IDR, DSR, diagnostic_capacity_daily, sequencing_capacity_daily) %>% arrange(total_pandemic_m)

df_s2$loc_name[df_s2$code=="OC_others"] <- "Oceania (others)"
df_s2$loc_name[df_s2$code=="AF_others"] <- "Africa (others)"
df_s2$loc_name[df_s2$code=="EU_others"] <- "Europe (others)"
df_s2$loc_name[df_s2$code=="AS_others"] <- "Asia (others)"
df_s2$loc_name[df_s2$code=="NA_others"] <- "North America (others)"
df_s2$loc_name[df_s2$code=="SA_others"] <- "South America (others)"
df_s2$continent[df_s2$code=="OC_others"] <- "Oceania"
df_s2$continent[df_s2$code=="AF_others"] <- "Africa"
df_s2$continent[df_s2$code=="EU_others"] <- "Europe"
df_s2$continent[df_s2$code=="AS_others"] <- "Asia"
df_s2$continent[df_s2$code=="NA_others"] <- "North America"
df_s2$continent[df_s2$code=="SA_others"] <- "South America"

df_s2 <- df_s2 %>% arrange(as.numeric(gsub("*", "", travel_hub_rank_total_pandemic, fixed=T)))
df_s2 <- df_s2 %>% mutate(population_m = round(population_m, 2), total_pandemic_m = round(total_pandemic_m, 2), total_normal_m = round(total_normal_m, 2), transit_rates = paste0(round(transit_rates*100,2), "%"), IDR = paste0(round(IDR*100,2), "%"), DSR = paste0(round(DSR*100,2), "%")) 

write_csv(df_s2, "results/figs/model_simulation/M3/Table_S2.csv")

## Table S4
### Table S4 should list two settings: the Omicron period and the 2019 period

df_all_values_fig3 <- list.files("results/model_data/model_simulation/M3/", "df_all_values_M3_\\d{1,2}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows() %>% filter(fig %in% c("fig_3e", "fig_3f"))
df_all_values_fig4 <- list.files("results/model_data/model_simulation/M4/", "df_all_values_M4_\\d{1,2}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows() %>% filter(id_unit_intro == 14)
df_all_values_fig4 <- df_all_values_fig4 %>% filter(!grepl("+", scenario, fixed=TRUE))

df_all_values_fig3$strategy_abbr <- sapply(strsplit(df_all_values_fig3$scenario, "_"),function(x) x[3])
df_all_values_fig3$strategy_abbr[grepl("^STp_", df_all_values_fig3$scenario)] <- gsub("^T", "P", df_all_values_fig3$strategy_abbr[grepl("^STp_", df_all_values_fig3$scenario)])
df_all_values_fig3$weight <- sapply(strsplit(df_all_values_fig3$scenario, "_"),function(x) x[4])
df_all_values_fig3 <- df_all_values_fig3 %>% filter(weight==0.5)

df_all_values_fig4$strategy_abbr <- sapply(strsplit(df_all_values_fig4$scenario, "_"),function(x) x[4])
df_all_values_fig4$strategy_abbr[grepl("^STp_", df_all_values_fig4$scenario)] <- gsub("^T", "P", df_all_values_fig4$strategy_abbr[grepl("^STp_", df_all_values_fig4$scenario)])
df_all_values_fig4$weight <- sapply(strsplit(df_all_values_fig4$scenario, "_"),function(x) x[5])
df_all_values_fig4 <- df_all_values_fig4 %>% filter(weight==0.5)

strategy_abbr_all <- unique(c(df_all_values_fig4$strategy_abbr, df_all_values_fig3$strategy_abbr))

df_s4 <- tibble(st_abbr = strategy_abbr_all)
df_s4$base_travel_hubs_num <- gsub("\\D", "", df_s4$st_abbr)

tmp <- bind_rows(df_all_values_fig3, ) %>% transmute(st_abbr=strategy_abbr, traveler_weight)
tmp$total_travel_hubs_num <- sapply(tmp$traveler_weight, function(x){sum(x!=999)})
tmp <- tmp %>% group_by(st_abbr) %>% select(st_abbr, total_travel_hubs_num) %>% unique() %>% summarise(total_travel_hubs_num_char = paste0(total_travel_hubs_num, collapse = "/")) %>% ungroup() %>% select(st_abbr, total_travel_hubs_num_char) %>% unique()
df_s4 <- left_join(df_s4 , tmp, by="st_abbr") 
df_s4$addition_travel_hubs_num <- mapply(function(x,y){paste0(x-y, collapse = "/")},lapply(strsplit(df_s4$total_travel_hubs_num_char, "/"), as.numeric), as.numeric(df_s4$base_travel_hubs_num))
df_s4$travel_hubs_ranking <- ifelse(grepl("^T", df_s4$st_abbr), "total", "per capita")
df_s4$with_continent <- ifelse(grepl("a-", df_s4$st_abbr), "Y", "N")
df_s4$surveillance_resources_type <- sapply(strsplit(df_s4$st_abbr, "-"),function(x) x[2])

get_mode <- function(v) {
  as.numeric(names(sort(table(unlist(v)), decreasing = TRUE))[1])
}

codes_all <- sort(unique(model_sims$unitname))

df_s4 <- left_join(
  df_s4, 
  df_all_values_fig3 %>% 
    transmute(
      st_abbr=strategy_abbr,
      info_Omicron = paste0(
        "Diagnostic efforts:",
        ifelse(
          sapply(change_IDR, function(x) all(x==1)),
          "All unchanged.\n",
          sapply(change_IDR, function(x) paste0(paste0(codes_all[unlist(x)-get_mode(x)>0.001], ":", round(x[unlist(x)-get_mode(x)>0.001]*100), "%", collapse = ", "), ", and others scaled down to XXX of the original level.\n"))
        ),
        "Sequencing efforts:",
        ifelse(
          sapply(change_DSR, function(x) all(x==1)),
          "All unchanged.",
          sapply(change_DSR, function(x) paste0(paste0(codes_all[unlist(x)-get_mode(x)>0.001], ":", round(x[unlist(x)-get_mode(x)>0.001]*100), "%", collapse = ", "), ", and others scaled down to XXX of the original level."))
        )
      )
    ),
  by="st_abbr"
)
df_s4$info_Omicron <- ifelse(df_s4$surveillance_resources_type=="H", gsub("XXX", "50%", df_s4$info_Omicron), gsub("XXX", "10%", df_s4$info_Omicron))

df_s4 <- left_join(
  df_s4, 
  df_all_values_fig4 %>% 
    transmute(
      st_abbr=strategy_abbr,
      info_2019 = paste0(
        "Diagnostic efforts: ",
        ifelse(
          sapply(change_IDR, function(x) all(x==1)),
          "All unchanged.\n",
          sapply(change_IDR, function(x) paste0(paste0(codes_all[unlist(x)-get_mode(x)>0.001], ":", round(x[unlist(x)-get_mode(x)>0.001]*100), "%", collapse = ", "), ", and others scaled down to XXX of the original level.\n"))
        ),
        "Sequencing efforts: ",
        ifelse(
          sapply(change_DSR, function(x) all(x==1)),
          "All unchanged.",
          sapply(change_DSR, function(x) paste0(paste0(codes_all[unlist(x)-get_mode(x)>0.001], ":", round(x[unlist(x)-get_mode(x)>0.001]*100), "%", collapse = ", "), ", and others scaled down to XXX of the original level."))
        )
      )
    ),
  by="st_abbr"
)
df_s4$info_2019 <- ifelse(df_s4$surveillance_resources_type=="H", gsub("XXX", "50%", df_s4$info_2019), gsub("XXX", "10%", df_s4$info_2019))

df_s4$base_travel_hubs_num <- as.numeric(df_s4$base_travel_hubs_num)
df_s4$info_travelhubs <- paste0(
  df_s4$total_travel_hubs_num,
  ": Top ",
  df_s4$base_travel_hubs_num,
  " (",
  ifelse(df_s4$travel_hubs_ranking=="total", "T-ranked", "P-ranked"),
  ") travel hubs ",
  ifelse(df_s4$with_continent=="Y", paste0("with ", df_s4$addition_travel_hubs_num, " additional continental hubs"), "")
  )

df_s4 %>% select(st_abbr, info_travelhubs, info_Omicron, info_2019) %>% write_csv("results/figs/model_simulation/M3/Table_S4.csv")

# check the IDR for OC_others
df_covar_m4 <- readRDS("results/model_data/data_covariates_M4.rds")
df_covar_m4 %>% filter(code=="OC_others") %>% select(date_decimal, infection_detection, detection_sequencing) %>% arrange(date_decimal) %>% unique() # DSR is 0, so that the change_IDR will not work

