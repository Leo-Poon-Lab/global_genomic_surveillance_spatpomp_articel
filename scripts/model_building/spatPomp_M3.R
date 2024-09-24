# This code builds the model M3: M0 with an adjustable sequencing propensity.
## The code is largely the same as M0, but differs in the rprocess_head and rprocess_tail parts.

library(spatPomp)
library(pomp)
library(tidyverse)

build_spatPomp_M3 <- function(
  global_data,
  dt_years = 1 / 365, 
  variants = c("Delta", "Omicron BA_one", "Omicron BA_two", "Others"),
  fixed_shared_params_values,
  fixed_unit_specific_params_values,
  N_introdution,
  for_IBPF = FALSE,
  change_IDR = NA, # can be a vector of length U, a single value, or NA
  change_DSR = NA, # can be a vector of length U, a single value, or NA
  sequencing_propensity = NA # can be a vector of length U, a single value, or NA
){
  stopifnot(!is.na(global_data))

  # check units
  units_all <- global_data$country_under_investigation
  U <- global_data$num_of_country_under_investigation
  ## make sure the units are in the same order as in the data (as implicitly assumed in the SpatPomp functions)
  stopifnot(all(units_all==unique(global_data$df_meas$code)))
  stopifnot(all(units_all==unique(global_data$df_covar$code)))

  # check variants
  strains_all <- sort(unique(global_data$gisaid$data_GISAID$lineage_new))
  stopifnot(all(strains_all==variants))
  stopifnot(variants[length(variants)]=="Others")
  num_of_strains <- length(strains_all)
  date_start_d <- global_data$dates$date_start_d

  # check M3 specific parameters
  if(length(change_IDR)==U){
    if(all(is.na(change_IDR))){
      change_IDR <- 1
    }
    stopifnot(all(change_IDR>=0))
  }
  if(length(change_DSR)==U){
    if(all(is.na(change_DSR))){
      change_DSR <- 1
    }
    stopifnot(all(change_DSR>=0))
  }
  if(length(sequencing_propensity)==U){
    if(all(is.na(sequencing_propensity))){
      sequencing_propensity <- 1
    }
    stopifnot(all(sequencing_propensity>=0 & sequencing_propensity<=1))
  }

  if(length(change_IDR)!=U){
    if(length(change_IDR)!=1){
      stopifnot(is.na(change_IDR))
    }
  }
  if(length(change_DSR)!=U){
    if(length(change_DSR)!=1){
      stopifnot(is.na(change_DSR))
    }
  }
  if(length(sequencing_propensity)!=U){
    if(length(sequencing_propensity)!=1){
      stopifnot(is.na(sequencing_propensity))
    }
  }

  # 3.3.1. Getting a measurements data.frame with columns for times, spatial units and measurements.
  # add origin-importation pairs for measurements
  # column names should be C_1_c, C_2_c, C_3_c, C_4_c, C_1_i_o{0:28}_, C_2_i_o{0:28}_, C_3_i_o{0:28}_, C_4_i_o{0:28}_
  data_measurements <- global_data$df_meas %>% filter(date_decimal>=date_start_d)
  data_measurements <- left_join(data_measurements, global_data$gisaid$data_GISAID_origin %>% ungroup() %>% mutate(date=`Collection date`, code=sequencing_country) %>% select(-date_decimal, -`Collection date`, -sequencing_country), c("date", "code"))
  data_measurements$code_num <- as.numeric(factor(data_measurements$code, levels=units_all))-1
  data_measurements$parent_country_num <- as.numeric(factor(data_measurements$parent_country, levels=units_all))-1
  data_measurements$lineage_new_num <- as.numeric(factor(data_measurements$lineage_new, levels=strains_all))

  data_measurements$case_origin <- mclapply(1:nrow(data_measurements), function(i){
    lineage_i <- data_measurements$lineage_new_num[i]
    code_num_i <- data_measurements$code_num[i]
    parent_country_num_i <- data_measurements$parent_country_num[i]

    if(parent_country_num_i==code_num_i){
      case_type_i <- "c"
    } else{
      case_type_i <- "i"
    }

    if(case_type_i=="c"){
      name_new <- paste0("C_", lineage_i, "_", case_type_i)
    } else{
      name_new <- paste0("C_", lineage_i, "_", case_type_i, "_o_", parent_country_num_i, "_")
    }
    return(name_new)
  }, mc.cores=8) %>% unlist()

  data_measurements <- data_measurements %>% select(-ends_with("_num")) 
  data_measurements <- data_measurements %>% select(-lineage_new, -parent_country)

  data_measurements$case_origin <- factor(data_measurements$case_origin, levels=c(paste0("C_", seq_along(strains_all), "_c"), paste0("C_", rep(seq_along(strains_all), each = U), "_i_o_", 0:(U-1), "_")), ordered = T)

  data_measurements <- data_measurements %>% pivot_wider(names_from = case_origin, values_from = N, values_fill = 0)

  data_measurements <- data_measurements %>% select(date, date_decimal, code, daily_cases, daily_deaths, all_of(c(paste0("C_", seq_along(strains_all), "_c"), paste0("C_", rep(seq_along(strains_all), each = U), "_i_o_", 0:(U-1), "_"))))

  if("NA" %in% names(data_measurements)){
    data_measurements <- data_measurements %>% select(-"NA")
  }
  data_measurements <- data_measurements %>% mutate_at(vars(daily_cases:`C_4_i_o_28_`), as.integer) # must be integer to avoid NAs in d_measure
  data_measurements <- data_measurements %>% select(-date)
  for(i in 1:num_of_strains){
    data_measurements[[paste0("C_", i, "_all")]] <- data_measurements %>% select(tidyselect::matches(paste0("C_", i, "_"))) %>% apply(1, sum)
  }

  # 3.3.2 Constructing model covariates.
  data_covariates <- global_data$df_covar %>% filter(date_decimal>=global_data$dates$date_start_d) %>% select(date_decimal:infection_fatality, cumulative_all_effectively_vaccinated) %>% select(-mandates_mean)

  ## add vaxx_rate (effective vaccination rate for all fully vaccinated individuals) (data from IHME "Fully vaccinated (one of one and†two†of two doses)")
  data_covariates_vaxx <- global_data$df_covar %>% filter(date_decimal>=global_data$dates$date_start_d) %>% group_by(code) %>% summarise(vacc_diff=max(cumulative_all_fully_vaccinated)-min(cumulative_all_fully_vaccinated), vaxx_num_per_year=vacc_diff/(max(date_decimal)-min(date_decimal)), vaxx_rate=vaxx_num_per_year/mean(population)) %>% select(code, vaxx_rate) %>% ungroup() 
  data_covariates <- left_join(data_covariates, data_covariates_vaxx)

  ## add total number of infection in the past 1-4 days and 5-10 days
  data_covariates_past_inf <- global_data$df_covar %>% filter(date_decimal<global_data$dates$date_start_d) %>% arrange(desc(date_decimal), code) %>% group_by(code) %>% summarise(inf_past_1_to_4=inf_cuml_mean[1]-inf_cuml_mean[5], inf_past_5_to_10=inf_cuml_mean[5]-inf_cuml_mean[10], population=mean(population)) %>% ungroup()
  stopifnot(!any(is.na(data_covariates_past_inf)))
  
  data_covariates_past_inf <- left_join(data_covariates_past_inf, global_data$travel$mov_mat %>% mutate(code=code_d) %>% group_by(code) %>% summarise(N_daily_arrival=mean(flow_all)), "code")
  data_covariates_past_inf$inf_past_1_to_4_imported <- data_covariates_past_inf$inf_past_1_to_4*(data_covariates_past_inf$N_daily_arrival*4/data_covariates_past_inf$population)
  data_covariates_past_inf$inf_past_1_to_4_community <- data_covariates_past_inf$inf_past_1_to_4-data_covariates_past_inf$inf_past_1_to_4_imported
  data_covariates_past_inf$inf_past_5_to_10_imported <- data_covariates_past_inf$inf_past_5_to_10*(data_covariates_past_inf$N_daily_arrival*6/data_covariates_past_inf$population)
  data_covariates_past_inf$inf_past_5_to_10_community <- data_covariates_past_inf$inf_past_5_to_10-data_covariates_past_inf$inf_past_5_to_10_imported

  ## The below part needs to be updated if studying in a different time period.
  data_covariates_variants <- global_data$gisaid$data_GISAID %>% filter(date_decimal<global_data$dates$date_start_d) %>% group_by(code, lineage_new) %>% summarise(N=sum(N)) %>% pivot_wider(names_from = lineage_new, values_from = N, values_fill = 0) %>% ungroup()
  if(any(units_all %in% data_covariates_variants$code)){
    missing_unit <- units_all[!units_all %in% data_covariates_variants$code]
    data_covariates_variants <- bind_rows(data_covariates_variants, tibble(code=missing_unit, Delta=1, `Omicron BA_one`=0, `Omicron BA_two`=0, Others=1))
  }
  data_covariates_variants$`Omicron BA_one` <- 0
  data_covariates_variants$`Omicron BA_two` <- 0

  data_covariates_past_inf <- left_join(data_covariates_past_inf, data_covariates_variants, "code")

  delta_prop <- mean(data_covariates_past_inf$Delta/data_covariates_past_inf$population, na.rm=T)
  others_prop <- mean(data_covariates_past_inf$Others/data_covariates_past_inf$population, na.rm=T)
  data_covariates_past_inf$Delta[is.na(data_covariates_past_inf$Delta)] <- data_covariates_past_inf$population[is.na(data_covariates_past_inf$Delta)]*delta_prop
  data_covariates_past_inf$Others[is.na(data_covariates_past_inf$Others)] <- data_covariates_past_inf$population[is.na(data_covariates_past_inf$Others)]*others_prop
  
  data_covariates_past_inf[,strains_all] <- t(apply(data_covariates_past_inf[,strains_all], 1, function(x){tmp=x/sum(x);if(is.na(sum(x))){return(x)}else{return(tmp)}}))
  ## The above part needs to be updated if studying in a different time period.

  ## for the number of infected in the past 10 days, the 1-4 days will be in the E compartment, and the 5-10 days will be in the I compartment.
  data_covariates_past_inf_sum <- data_covariates_past_inf %>% group_by(code) %>% transmute(
    E_1_i_init=inf_past_1_to_4_imported*Delta,
    E_2_i_init=inf_past_1_to_4_imported*`Omicron BA_one`,
    E_3_i_init=inf_past_1_to_4_imported*`Omicron BA_two`,
    E_4_i_init=inf_past_1_to_4_imported*Others,
    E_1_c_init=inf_past_1_to_4_community*Delta,
    E_2_c_init=inf_past_1_to_4_community*`Omicron BA_one`,
    E_3_c_init=inf_past_1_to_4_community*`Omicron BA_two`,
    E_4_c_init=inf_past_1_to_4_community*Others,
    I_1_i_init=inf_past_5_to_10_imported*Delta,
    I_2_i_init=inf_past_5_to_10_imported*`Omicron BA_one`,
    I_3_i_init=inf_past_5_to_10_imported*`Omicron BA_two`,
    I_4_i_init=inf_past_5_to_10_imported*Others,
    I_1_c_init=inf_past_5_to_10_community*Delta,
    I_2_c_init=inf_past_5_to_10_community*`Omicron BA_one`,
    I_3_c_init=inf_past_5_to_10_community*`Omicron BA_two`,
    I_4_c_init=inf_past_5_to_10_community*Others
  ) %>% ungroup()

  data_covariates <- left_join(data_covariates, data_covariates_past_inf_sum, "code")
  stopifnot(!any(is.na(data_covariates)))

  ## The fully vaccinated population and previous infected population.
  data_covariates <- data_covariates %>% select(-cumulative_all_effectively_vaccinated)
  data_covariates_immunised <- global_data$df_covar %>% filter(date_decimal<global_data$dates$date_start_d) %>% filter(date_decimal==min(date_decimal)) %>% select(code, population, cumulative_all_fully_vaccinated,  inf_cuml_mean_unvax)
  data_covariates_immunised <- data_covariates_immunised %>% mutate(V_init=cumulative_all_fully_vaccinated+inf_cuml_mean_unvax)
  data_covariates <- left_join(data_covariates, data_covariates_immunised %>% select(code, V_init), "code")

  ## S_init = population - all other compartments
  data_covariates$R_init <- 0
  data_covariates <- data_covariates %>% mutate(
    S_init = population - V_init - R_init - E_1_i_init - E_2_i_init - E_3_i_init - E_4_i_init - E_1_c_init - E_2_c_init - E_3_c_init - E_4_c_init - I_1_i_init - I_2_i_init - I_3_i_init - I_4_i_init - I_1_c_init - I_2_c_init - I_3_c_init - I_4_c_init 
  )
  check <- data_covariates$S_init < 0
  data_covariates$V_init[check] <- data_covariates$population[check] - (data_covariates[check,] %>% select(E_1_i_init:I_4_c_init) %>% apply(1, sum))
  data_covariates$S_init[check] <- 0
  stopifnot(all(abs(data_covariates$population-(data_covariates %>% select(contains("init")) %>% apply(1, sum)))<1))

  ## Add detection_sequencing ratio
  df_dec_seq <- data_measurements %>% group_by(date_decimal, code) %>% summarise(detection_sequencing=sum(c_across(C_1_c:C_4_i_o_28_),na.rm = T)/(daily_cases+daily_deaths))
  df_dec_seq$detection_sequencing[df_dec_seq$detection_sequencing>1] <- 1

  data_covariates <- left_join(data_covariates, df_dec_seq, c("date_decimal", "code"))

  if(any(is.infinite(data_covariates$detection_sequencing))){
    data_covariates$detection_sequencing[is.infinite(data_covariates$detection_sequencing)]<- 0
  }

  ## Add M3 specific covariates
  if(all(!is.na(change_IDR))){
    data_covariates <- left_join(data_covariates, tibble(code=units_all, change_IDR=change_IDR), "code")
  } else{
    warning("change_IDR is not complete, set to 1")
    data_covariates$change_IDR <- 1
  }
  if(all(!is.na(change_DSR))){
    data_covariates <- left_join(data_covariates, tibble(code=units_all, change_DSR=change_DSR), "code")
  } else{
    warning("change_DSR is not complete, set to 1")
    data_covariates$change_DSR <- 1
  }
  if(all(!is.na(sequencing_propensity))){
    data_covariates <- left_join(data_covariates, tibble(code=units_all, sequencing_propensity=sequencing_propensity), "code")
  } else{
    warning("sequencing_propensity is not complete, set to 1")
    data_covariates$sequencing_propensity <- 999 # magic number meaning that the sequencing propensity is not used
  }

  # 3.3.2 Constructing model components (latent state initializer, latent state transition simulator and measurement model).
  ## List all state-names in pomp object:
  source("scripts/model_building/helper/unit_state_names_M0.R") 
  unit_state_names <- get_unit_state_names(num_of_strains, U)# load unit_state_names

  ## Get params_fixed and params_est
  source("scripts/model_fitting/helper/parameter_names_est.R") # load shared_param_names_est, unit_specific_names_est
  params_RP_shared <- unique(c(
    names(fixed_shared_params_values),
    shared_param_names_est
  ))
  params_RP_diff <- c(
    names(fixed_unit_specific_params_values),
    unit_specific_names_est
  )
  params_fixed <- c(names(fixed_shared_params_values), names(fixed_unit_specific_params_values))
  stopifnot(all(params_fixed %in% c(params_RP_shared, params_RP_diff)))

  # Prepare the globals for the pomp object
  # globals: optional character or C snippet; arbitrary C code that will be hard-coded into the shared-object library created when C snippets are provided. If no C snippets are used, globals has no effect.
  # globals in our model should include: mov_mat

  set_unit_specific <- Csnippet(paste0("const int ", params_RP_diff,
    "_unit = 1;\n", collapse=" ")) # unit_specific parameters should refer to the corresponding (the n=n*1) element.

  if(for_IBPF){
    set_shared <- Csnippet(paste0("const int ", params_RP_shared,
    "_unit = 1;\n", collapse=" ")) # shared parameters should only refer to the corresponding (the n=n*1) element.
  } else {
    set_shared <- Csnippet(paste0("const int ", params_RP_shared,
    "_unit = 0;\n", collapse=" ")) # shared parameters should only refer to the first (the 0=n*0) element.
  }
  
  set_t_start_d <- Csnippet(paste0("const double t_start_d = ", date_start_d, ";\n"))

  model_globals <- Csnippet(
    paste(global_data$travel$v_by_g_C, set_unit_specific, set_shared, set_t_start_d, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  model_paramnames <- c(
    if(length(params_RP_shared)>0){
      paste0(rep(params_RP_shared, each=U), 1:U) # to be compatible with the use of ibpf
    },
    if(length(params_RP_diff)>0){
      paste0(rep(params_RP_diff, each=U), 1:U)
    }
  )

  # rinit
  source("scripts/model_building/helper/customize_rinit_M0.R")
  code_rinit_M0 <- get_code_rinit_M0(U=U) # load code_rinit_M0
  model_rinit <- spatPomp_Csnippet(
    unit_statenames = unit_state_names,
    unit_covarnames = c("S_init", "V_init", 
    "E_1_i_init", "E_2_i_init", "E_3_i_init", "E_4_i_init", "E_1_c_init", "E_2_c_init", "E_3_c_init", "E_4_c_init", "I_1_i_init", "I_2_i_init", "I_3_i_init", "I_4_i_init", "I_1_c_init", "I_2_c_init", "I_3_c_init", "I_4_c_init", "R_init"),
    code=code_rinit_M0
  )

  # rprocess 
  code_rprocess_head_M3 <- paste0(readLines("scripts/model_building/helper/code_rprocess_head_M3.C"), collapse = "\n") # load code_rprocess_head_M3, this will be changed in the model for simulation
  code_rprocess_main_M0 <- paste0(readLines("scripts/model_building/helper/code_rprocess_main_M0.C"), collapse = "\n") # load code_rprocess_main_M0, in the model for simulation, the code_rprocess_main_M0 will be reused as it is unchanged.
  code_rprocess_tail_M3 <- paste0(readLines("scripts/model_building/helper/code_rprocess_tail_M3.C"), collapse = "\n") # load code_rprocess_tail_M3, this will be changed in the model for simulation
  code_rprocess <- paste0(code_rprocess_head_M3, code_rprocess_main_M0, code_rprocess_tail_M3)

  source("scripts/model_building/helper/customize_rprocess_M3.R")
  code_rprocess <- customize_rprocess_M3(code_rprocess, U, units_all, num_of_strains, N_introdution)

  model_rprocess <- spatPomp_Csnippet(
    unit_statenames = unit_state_names,
    unit_covarnames = c('population', 'vaxx_rate', 'travel_control_level', 'infection_fatality', 'infection_detection', "detection_sequencing", "change_IDR", "change_DSR", "sequencing_propensity"),
    unit_paramnames = c(params_RP_shared[!grepl("theta_", params_RP_shared)], params_RP_diff), # all params except those in the measurement process 
    code=code_rprocess
  )

  # # dmeasure
  # source("scripts/model_building/helper/customize_dmeasure_M0.R") 
  # code_dmeasure_M0 <- get_code_dmeasure_M0(U, num_of_strains, num_of_units = U) # load code_dmeasure_M0
  # model_dmeasure <- spatPomp_Csnippet(
  #   unit_statenames = c(
  #     "C_reported_new", "D_new",
  #     unit_state_names[grepl("_sequenced.+new$", unit_state_names)]
  #     ),
  #   unit_obsnames = c(
  #     "daily_cases", "daily_deaths", 
  #     paste0("C_", seq_len(num_of_strains), "_c"),
  #     paste0("C_", seq_len(num_of_strains), "_i_o_", rep(seq_len(U)-1, each=num_of_strains), "_"),
  #     paste0("C_", seq_len(num_of_strains), "_all")
  #     ),
  #   unit_paramnames = c("tau_seq_all", "tau_seq_unit", "tau_cases", "tau_deaths"),
  #   code = code_dmeasure_M0
  # )

  # dunit_measure
  source("scripts/model_building/helper/customize_dunit_measure_M0.R") 
  code_dunit_measure_M0 <- get_code_dunit_measure_M0(U, num_of_strains, num_of_units = U)
  # load code_dunit_measure_M0
  model_dunit_measure <- spatPomp_Csnippet(
    unit_paramnames = c("tau_seq_all", "tau_seq_unit", "tau_cases", "tau_deaths"),
    code = code_dunit_measure_M0
  )

  # rmeasure
  source("scripts/model_building/helper/customize_rmeasure_M0.R")
  code_rmeasure_M0 <- get_code_rmeasure_M0(U, num_of_strains) # load code_rmeasure_M0
  model_rmeasure <- spatPomp_Csnippet(
    method='rmeasure',
    unit_statenames = c(
      "C_reported_new", "D_new",
      unit_state_names[grepl("_sequenced.+new$", unit_state_names)]
      ),
    unit_obsnames = c(
      "daily_cases", "daily_deaths", 
      paste0("C_", seq_len(num_of_strains), "_c"),
      paste0("C_", seq_len(num_of_strains), "_i_o_", rep(seq_len(U)-1, each=num_of_strains), "_"),
      paste0("C_", seq_len(num_of_strains), "_all")
      ),
    unit_paramnames = c("tau_seq_all", "tau_seq_unit", "tau_cases", "tau_deaths"),
    code = code_rmeasure_M0
  )

  ## vunit_measure, eunit_measure and skeleton not defined here, as ibpf does not require them?

  # parameter transformation to limit its value range
  basic_logit_names <- c(params_RP_shared, params_RP_diff)[!c(params_RP_shared, params_RP_diff) %in% params_fixed] # values from 0 to 1

  logit_names <- unlist(lapply(basic_logit_names, function(x,U) paste0(x,1:U),U))
  logit_names <- logit_names[logit_names %in% model_paramnames]

  model_partrans <- parameter_trans(logit=logit_names)

  # saveRDS(data_covariates, file = here::here("results/model_data/data_covariates_M3.rds"))

  # 3.4. Bringing all the data and model components together to form a spatPomp object via a call to spatPomp().
  model_spatpomp <- spatPomp(
    data = as.data.frame(data_measurements),
    units = "code",
    times = "date_decimal",
    covar = data_covariates,
    t0 = date_start_d,
    unit_statenames = unit_state_names,
    rprocess=euler(model_rprocess, delta.t=dt_years),
    unit_accumvars = unit_state_names[grepl("_new$", unit_state_names)],
    paramnames=model_paramnames,
    partrans=model_partrans,
    globals=model_globals,
    rinit=model_rinit,
    # dmeasure=model_dmeasure,
    rmeasure=model_rmeasure,
    dunit_measure=model_dunit_measure
  )

  model_params <- rep(0,length=length(model_paramnames))
  names(model_params) <- model_paramnames
  for(p in names(fixed_shared_params_values)) {model_params[paste0(p,seq_len(U))] <- as.numeric(fixed_shared_params_values[p])}
  for(p in names(fixed_unit_specific_params_values)) {model_params[paste0(p,seq_len(U))] <- as.numeric(fixed_unit_specific_params_values[[p]])}

  coef(model_spatpomp) <- model_params
  return(model_spatpomp)
}