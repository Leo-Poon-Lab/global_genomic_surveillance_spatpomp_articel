library(tidyverse)

# Load some common data
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

# Fig 1 
dir_rst_fig1 <- "results/figs/model_simulation/Omicron20/"
## Fig 1A
## May need to show the table for total travel volume and travel per capita, as supplementary material.
df_travel_volume <- read_csv("results/model_data/data_travel_hubs_order.csv")
df_travel_volume %>% left_join(cross_check_table %>% select(code, loc_name), "code") %>% select(-code_d) %>% mutate(loc_name = ifelse(is.na(loc_name), code, loc_name)) %>% select(code, loc_name, population, total_travel_volume=total, travel_per_capita) %>% arrange(desc(travel_per_capita)) %>% mutate(travel_hub_rank = row_number()) %>% write_csv(paste0(dir_rst_fig1, "travel_hub_rank.csv"))

## extract transit passenger data
df_fitting <- readRDS("results/model_data/data_fitting_Omicron20.rds")
df_transit <- tibble(code=df_fitting$country_under_investigation, transit_rates=df_fitting$transit_rates)
write_csv(df_transit, paste0(dir_rst_fig1, "transit_rates.csv"))

## Fig 1B, nothing to do here

## Fig 1C, nothing to do here