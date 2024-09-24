library(tidyverse)

# Fig 3

## Fig 3a
## Where are the Omicron BA.1 and BA.2 variants first identified?
### For IDR/DSR=1, mainly (>99%) in South Africa

## Fig 3c
## Why sequencing more community cases works better in BA.1 than BA.2?
### Check the IDR and DSR
data_covar <- readRDS("results/model_data/data_covariates_Omicron20.rds")

lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
df_vline <- tibble(lag=lags_intro)
df_vline$date_decimal <- floor(df_vline$lag)/365.25 + min(data_covar$date_decimal)
df_vline$date <- date_decimal(df_vline$date_decimal)

data_covar %>% filter(code=="ZA") %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date, infection_detection, detection_sequencing) %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+14) %>% print(n=200)

data_covar %>% filter(code=="ZA") %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date, infection_detection, detection_sequencing) %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+14) %>% print(n=200)

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

model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+20) %>% print(n=30)
model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_one"]]), date<=as_date(df_vline$date[["day_BA_one"]])+20) %>% transmute(date, C_2_c_i_prop=C_2_c_i/C_k_c_i, C_3_c_i_prop=C_3_c_i/C_k_c_i) %>% print(n=30)

model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+20) %>% print(n=30)
model_sims_median %>% filter(date>=as_date(df_vline$date[["day_BA_two"]]), date<=as_date(df_vline$date[["day_BA_two"]])+20) %>% transmute(date, C_2_c_i_prop=C_2_c_i/C_k_c_i, C_3_c_i_prop=C_3_c_i/C_k_c_i) %>% print(n=30)


### This is because the diagnostic-sequencing ratio is too low during the Omicron BA.2 emergence period ()
