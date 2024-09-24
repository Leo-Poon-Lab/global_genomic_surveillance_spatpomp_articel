library(tidyverse)
library(countrycode)

# IHME estimates
files_ihme_estimates <- list.files(here::here("data/IHME/"), "file_reference", full.names = T)
df_ihme <- lapply(files_ihme_estimates, read_csv)
df_ihme <- bind_rows(df_ihme)
df_ihme <- df_ihme %>% select(
  date, # Date of prediction
  location_id, # Location ID code
  location_name, # Location name
  population, # Population size
  cumulative_all_effectively_vaccinated, # Effectively vaccinated (one and two dose with efficacy)
  cumulative_all_fully_vaccinated, # Fully vaccinated (one of one and†two†of two doses)
  cumulative_all_vaccinated, # Initially vaccinated (one dose of two doses)
  cumulative_cases, # Cumulative cases (raw data)
  infection_detection, # Infection detection ratio
  infection_fatality, # Infection fatality ratio
  daily_cases, # Daily cases (raw data)
  daily_deaths, # Daily deaths (raw data with excess mortality scalar applied)
  daily_deaths_unscaled, # Daily deaths (raw data without excess mortality scalar applied)
  cases_mean, # Daily cases (mean estimate)
  inf_mean, # Daily infections (mean estimate)
  mandates_mean, # Mandate intensity across 17 non-pharmaceutical intervention variables on a scale of 0 to 1 (mean estimate)
  seir_daily_unscaled_mean, # Daily reported deaths (mean estimate)
  inf_cuml_mean, inf_cuml_mean_unvax, inf_cuml_mean_vax, # Cumulative infections (mean estimate)
  NULL
)
df_ihme <- df_ihme %>% filter(location_name!="Global")

df_ihme_supp <- read_csv(here::here("data/IHME/IHME_COVID_19_IES_2019_2021_COVID_19_MORT_Y2022M04D08.CSV"))
df_ihme_supp <- df_ihme_supp %>% select(location_id, location_name, level) %>% unique()

# make the country and location code consistent with the airport data #
# iso country and region name data #
df_all_loc_name <- df_ihme %>% select(location_id, location_name) %>% unique()
df_all_loc_name$location_name[duplicated(df_all_loc_name$location_name)] # "Georgia" "Punjab" are duplicated

df_countries <- read_csv(here::here("data/our_airports_data/countries.csv"), col_types=cols(.default = "c"), na="") %>% mutate(iso_country=code, country=name) %>% select(iso_country, country, everything()) %>% select(-name, -code)
df_regions <- read_csv(here::here("data/our_airports_data/regions.csv"), col_types=cols(.default = "c"), na="") %>% mutate(iso_region=code, region=name) %>% select(iso_country, local_code, iso_region, everything()) %>% select(-name, -code)

# English to ISO #
cross_check_table <- tibble(loc_id = df_all_loc_name$location_id, loc_name = df_all_loc_name$location_name, country_iso2c = countrycode(loc_name, origin = 'country.name', destination = 'iso2c'), country_iso3c = countrycode(loc_name, origin = 'country.name', destination = 'iso3c'), )
cross_check_table <- left_join(cross_check_table, df_countries %>% mutate(country_iso2c=iso_country) %>% select(country_iso2c, continent), "country_iso2c")
cross_check_table <- left_join(cross_check_table, df_regions %>% mutate(loc_name=region) %>% select(loc_name, iso_region), "loc_name")
cross_check_table %>% filter(loc_id %in% names(table(cross_check_table$loc_id)[table(cross_check_table$loc_id)>1]))
sum(!is.na(cross_check_table$country_iso2c))
sum(is.na(cross_check_table$country_iso2c))
cross_check_table <- left_join(cross_check_table, df_ihme_supp %>% mutate(loc_id=location_id, loc_name=location_name) %>% select(loc_id, loc_name, level), c("loc_id", "loc_name"))

writexl::write_xlsx(cross_check_table, here::here("data/our_airports_data/cross_check_table_ihme_input.xlsx"))

# after manual completion of the cross-check table, we load the data back #
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx"))

# remove "China" and only keeps "China (without Hong Kong and Macao)" #
if(any(unique(df_ihme$location_name)=="China (without Hong Kong and Macao)")){
  df_ihme <- df_ihme %>% filter(location_name != "China")
}

# remove "Provincia autonoma di Trento" which has no airport
df_ihme <- df_ihme %>% filter(location_name != "Provincia autonoma di Trento")

# paste the unified iso location code back #
stopifnot(all(unique(df_ihme$location_id) %in% cross_check_table$loc_id))
df_ihme <- left_join(df_ihme, cross_check_table %>% mutate(location_id=loc_id) %>% select(location_id, code, level))
saveRDS(df_ihme, here::here("data/IHME/df_ihme.rds"))

