#!/usr/bin/Rscript

#### Note: This script is used to fit the Omicron model ####
source(here::here("scripts/data_processing/install_prerequisite.R"))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(slider))
suppressPackageStartupMessages(library(doParallel))
setwd(here::here())
check_redo_data_fitting <- FALSE
model_name <- "Omicron20"
variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others")
# colors_lineage <- MetBrewer::met.brewer("Degas", 4, "continuous")
colors_lineage <- RColorBrewer::brewer.pal(4, "Set1")
names(colors_lineage) <- variants

# modeling the Omicron variant emergence and circulation in the first 4 months #
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

# This step outputs the data required for model fitting.
preprocessing <- function(){
  ## Prepare data for model fitting ##
  source(here::here("scripts/model_building/global_data_prepare.R"))
  data_fitting <- global_data_prepare(
    date_start = "2021-09-15",
    date_end = "2022-03-15",
    U = NA,
    U_each_continent = 3, # the number of top U countries to be included in each continent
    top_by = "cases and flow", # "cases" or "flow" or "cases and flow", sort the top U countries by cases or traffic flow
    threshold = 0,
    threshold_country = "HK", # Hong Kong to be used as the threshold
    global_scale = "country", # "continent", "country", "country_state"
    lineage_top_n = 3,
    add_HK = TRUE # force to add Hong Kong to the data
  )

  # print("country_under_investigation: ")
  # print(data_fitting$country_under_investigation)

  ## Flight maps ##
  # if (!exists("support_ggsave")) {
  #   support_ggsave <- TRUE
  # }
  # if (support_ggsave) {
  #   ggsave(here::here(paste0("results/figs/flights_map_", model_name, ".jpg")), plot=data_fitting$plot_airline, width=12, height=12*0.618, dpi = 450)
  # }

  ## Merge lineages ##
  ### you can merge similar lineages e.g., BA.1 and BA.1.1 into one lineage ###
  lineages_all <- names(data_fitting$gisaid$lineages)
  lineages_new <- lineages_all
  lineages_new[grepl("^BA.1", lineages_new)] <- "Omicron BA_one"
  lineages_new[grepl("^BA.2", lineages_new)] <- "Omicron BA_two"
  lineages_new[grepl("^AY", lineages_new)] <- "Delta"
  lineages_new[grepl("^B.1.617", lineages_new)] <- "Delta"
  lineages_new[!lineages_new %in% c("Omicron BA_one", "Omicron BA_two", "Delta")] <- "Others"

  data_fitting$gisaid$data_GISAID <- data_fitting$gisaid$data_GISAID %>% mutate(lineage_new = factor(lineage, levels=lineages_all, labels=lineages_new))
  data_fitting$gisaid$data_GISAID_origin <- data_fitting$gisaid$data_GISAID_origin %>% mutate(lineage_new = factor(lineage, levels=lineages_all, labels=lineages_new))

  data_fitting$gisaid$data_GISAID <- data_fitting$gisaid$data_GISAID %>% group_by(`Collection date`, date_decimal, code, lineage_new) %>% summarise(N=sum(N)) %>% ungroup()
  data_fitting$gisaid$data_GISAID_origin <- data_fitting$gisaid$data_GISAID_origin %>% group_by(`Collection date`, date_decimal, sequencing_country, parent_country, lineage_new) %>% summarise(N=sum(N)) %>% ungroup()

  data_fitting$gisaid$data_GISAID$lineage_new <- as.character(data_fitting$gisaid$data_GISAID$lineage_new)
  data_fitting$gisaid$data_GISAID_origin$lineage_new <- as.character(data_fitting$gisaid$data_GISAID_origin$lineage_new)

  source(here::here("scripts/model_fitting/helper/plots_preprocessing.R"))
  # plot_gisaid(data_fitting, outprefix="GISAID_raw_")

  ## calculate moving average of the lineages
  window_size=7

  data_fitting$gisaid$data_GISAID <- cross_join(tibble(`Collection date`=seq(min(data_fitting$gisaid$data_GISAID$`Collection date`), max(data_fitting$gisaid$data_GISAID$`Collection date`), by="day")), tibble(code=data_fitting$country_under_investigation)) %>% cross_join(tibble(lineage_new=unique(lineages_new))) %>% left_join(data_fitting$gisaid$data_GISAID, by=c("code", "lineage_new", "Collection date")) %>% mutate(N=ifelse(is.na(N), 0, N))
  data_fitting$gisaid$data_GISAID <- data_fitting$gisaid$data_GISAID %>% group_by(code, lineage_new) %>% mutate(N = map_dbl(seq_along(N), ~mean(N[max(1, .x - window_size + 1):.x])))

  data_fitting$gisaid$data_GISAID_origin <- cross_join(tibble(`Collection date`=seq(min(data_fitting$gisaid$data_GISAID_origin$`Collection date`), max(data_fitting$gisaid$data_GISAID_origin$`Collection date`), by="day")), tibble(sequencing_country=data_fitting$country_under_investigation)) %>% cross_join(tibble(parent_country=data_fitting$country_under_investigation)) %>% cross_join(tibble(lineage_new=unique(lineages_new))) %>% left_join(data_fitting$gisaid$data_GISAID_origin, by=c("sequencing_country", "parent_country", "lineage_new", "Collection date")) %>% mutate(N=ifelse(is.na(N), 0, N))
  data_fitting$gisaid$data_GISAID_origin <- data_fitting$gisaid$data_GISAID_origin %>% group_by(sequencing_country, parent_country, lineage_new) %>% mutate(N = map_dbl(seq_along(N), ~mean(N[max(1, .x - window_size + 1):.x])))

  # plot_gisaid(data_fitting, outprefix="GISAID_mov_avg_")
  # plot_in_out_flow(data_fitting)

  return(data_fitting)
}

data_fitting_rds_path <- paste0("results/model_data/data_fitting_", model_name, ".rds")
if(check_redo_data_fitting){file.remove(data_fitting_rds_path)}
if(!file.exists(data_fitting_rds_path)){
  data_fitting <- preprocessing()
  saveRDS(data_fitting, data_fitting_rds_path)
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

# fixed parameters
source("scripts/model_fitting/helper/fixed_parameters.R")

## This step do Model building ##
source(here::here(paste0("scripts/model_building/spatPomp_", model_name, ".R")))

if (!exists("task_id")) {
  task_id <- 1
}
if (!exists("rebuild")) {
  rebuild <- T
}
if (task_id == 1 && rebuild == T){
  model_omicron20 <- build_global_spatPomp_Omicron20(
    global_data=data_fitting,
    dt_years = 1 / 365, 
    variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others"),
    fixed_shared_params_values = fixed_shared_params_values,
    fixed_unit_specific_params_values = fixed_unit_specific_params_values,
    N_introdution = 2000, # the number of Omicron BA.1 cases introduced into E and I compartments in South Africa, 1000 correspond to 3% of the infected population on day 0
    for_IBPF = FALSE
  )
  
  # params for testing
  if(any(coef(model_omicron20)==0)){
    coef(model_omicron20)[coef(model_omicron20)==0] <- 0.5
  }

  # # For testing
  # time_used <- system.time(print(paste("bpfilter logLik for", model_name, "model:", logLik(bpfilter(model_omicron20, Np=1, block_size=1)))))
  # print(paste0("Time used for logLik bpfilter:"))
  # print(time_used)
  
  saveRDS(model_omicron20, paste0("results/model_data/model_", model_name, ".rds"))
}
