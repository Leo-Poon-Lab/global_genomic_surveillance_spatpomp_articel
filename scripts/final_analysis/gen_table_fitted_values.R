library(tidyverse)

df_1 <- read_csv("results/model_data/profiling_HKU/Omicron20/profiling_unit_10/df_mcap_transformed_high.csv")
df_2 <- read_csv("results/model_data/profiling_HKU/Omicron20/profiling_shared_11/df_mcap_transformed_high.csv")

df_combined <- bind_rows(df_1, df_2) %>% filter(!is.na(profile_para))

df_combined$text_note <- paste0(round(df_combined$mle_transformed, 2), " (", round(df_combined$CI_low_transformed, 2), ", ", round(df_combined$CI_high_transformed, 2), ")")

df_combined %>%
  select(profile_para, text_note, everything()) %>%
  write_csv("results/model_data/profiling_HKU/Omicron20/df_fitted.csv")

source("scripts/model_fitting/helper/tranform_log_params.R")
# for beta_naught_29
k_values["beta_naught"] * 0.1926539 + b_values["beta_naught"]

# translate the Inverse of the incubation period to the incubation period
365/83.14
365/(k_values["epsilon_four"]*c(0, 1) + b_values["epsilon_four"])

