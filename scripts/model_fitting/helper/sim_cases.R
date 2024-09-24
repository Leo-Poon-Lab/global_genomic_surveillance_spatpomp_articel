suppressPackageStartupMessages(library(doParallel))
source("scripts/model_simulation/helper/get_dates.R")

# multi-core
simulate_step_mc <- function(
  model,
  model_nsim,
  n_cores,
  if_aggregate_sequenced_cases=TRUE,
  if_return_EDT=FALSE,
  variants_for_EDT=c("Omicron BA.1", "Omicron BA.2"),
  df_vline=NA
) {
  stopifnot(model_nsim > 0)
  if(if_return_EDT){
    stopifnot(!is.na(df_vline))
  }

  doParallel::registerDoParallel(n_cores)
  doRNG::registerDoRNG(20231118)

  model_sims <- foreach(i = 1:model_nsim, .combine = rbind, .packages = 'spatPomp') %dopar% {
    # Simulate from the fitted version of Model
    model_sims <- simulate(model, nsim = 1, format = 'data.frame')
    model_sims$`.id` <- as.character(model_sims$`.id`)
    model_sims$`.id` <- as.character(i)
    # Aggregate the sequenced cases
    if(if_aggregate_sequenced_cases){
      model_sims <- aggregating_sequenced_cases(model_sims)
    }
    if(if_return_EDT){
      model_sims <- get_EDT(model_sims, df_vline, variants_for_EDT)
    }
    model_sims
  }
  doParallel::stopImplicitCluster()

  return(model_sims)
}

# single core
simulate_step_sc <- function(
  model,
  model_nsim=2000,
  if_aggregate_sequenced_cases=TRUE
) {
  # Simulate from the fitted version of Model
  model_sims <- simulate(model, nsim = model_nsim, format = 'data.frame')
  if(if_aggregate_sequenced_cases){
    model_sims <- aggregating_sequenced_cases(model_sims)
  } else {
    model_sims
  }
}

get_quants_from_sim <- function(
  data_sim,
  events=c("daily_cases", "C_1_all", "C_2_all", "C_3_all", "C_4_all"),
  quants=c(0.05, 0.5, 0.95)
){
  # get quants and reframe
  data_quants <- data_sim %>% select(date_decimal, unitname, all_of(events)) %>%
    group_by(unitname, date_decimal) %>%
    reframe(across(all_of(events), \(x) quantile(x, probs = quants, na.rm = TRUE))) %>%
    ungroup()
  # add a column specifying the exact quantile for each row
  data_quants %>% mutate(quantile_level = rep(quants, nrow(data_quants) / length(quants)))
}

aggregating_sequenced_cases <- function(model_sims, events=c("C_1_sequenced", "C_2_sequenced", "C_3_sequenced", "C_4_sequenced"), patterns = c("C_1_._sequenced_.+", "C_2_._sequenced_.+", "C_3_._sequenced_.+", "C_4_._sequenced_.+")){
  if("C_1_all_sequenced_new" %in% colnames(model_sims)){
    for(i in seq_along(events)){
      event <- events[i]
      column_name_event <- gsub("_sequenced", "_all_sequenced_new", event, fixed = T)
      model_sims[[event]] <- model_sims[[column_name_event]]
    }
  } else {
    for(i in seq_along(events)){
      event <- events[i]
      pattern_event <- patterns[i]
      model_sims[[event]] <- model_sims %>% select(tidyselect::matches(pattern_event)) %>% apply(1, sum)
    }
  }
  # Aggregate the sequenced cases
  model_sims
}

aggregating_measurements <- function(df_meas){
  if("C_1_all" %in% colnames(df_meas)){
    df_meas$Delta <- df_meas$C_1_all
    df_meas$Omicron_BA_one <- df_meas$C_2_all
    df_meas$Omicron_BA_two <- df_meas$C_3_all
    df_meas$Others <- df_meas$C_4_all
  } else {
    df_meas$Delta <- df_meas %>% select(tidyselect::matches("C_1_.+")) %>% apply(1, sum)
    df_meas$Omicron_BA_one <- df_meas %>% select(tidyselect::matches("C_2_.+")) %>% apply(1, sum)
    df_meas$Omicron_BA_two <- df_meas %>% select(tidyselect::matches("C_3_.+")) %>% apply(1, sum)
    df_meas$Others <- df_meas %>% select(tidyselect::matches("C_4_.+")) %>% apply(1, sum)
  }
  df_meas$Cases <- df_meas$daily_cases
  df_meas$Deaths <- df_meas$daily_deaths
  df_meas
}
