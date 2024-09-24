get_EDT <- function( # earliest detection time (sequenced cases) for model_sims, for BA1 and BA2
  model_sims,
  df_vline,
  variants = c("Omicron BA.1", "Omicron BA.2")
){ 
  if(length(variants)==1 & all(variants=="Omicron BA.1")){
    df_sims <- model_sims %>% transmute(
      date = date(date_decimal(date_decimal)),
      unitname,
      id_sim = `.id`,
      C_2_c_sequenced_new,
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_sequenced_origin_")))
    )
  } else if (length(variants)==2 & all(variants %in% c("Omicron BA.1", "Omicron BA.2"))){
    df_sims <- model_sims %>% transmute(
      date = date(date_decimal(date_decimal)),
      unitname,
      id_sim = `.id`,
      C_2_c_sequenced_new,
      C_3_c_sequenced_new,
      C_2_i_sequenced_new = rowSums(select(., starts_with("C_2_i_sequenced_origin_"))),
      C_3_i_sequenced_new = rowSums(select(., starts_with("C_3_i_sequenced_origin_")))
    )
  }

  df_sims <- df_sims %>% pivot_longer(
    cols = starts_with("C_"),
    names_to = "event",
    values_to = "value"
  ) %>% filter(
    value >= 1
  ) %>% separate(
    col = event,
    into = c("C", "variant", "case_type", "sequenced", "new"),
    sep = "_",
    remove = FALSE
  ) %>% select(
    -C,
    -sequenced,
    -event,
    -new
  ) %>% arrange(
    date, desc(value)
  ) %>% group_by(
    variant, case_type, id_sim,
  ) %>% slice_head() %>% mutate(
    variant = factor(variant, levels=2:3, labels=c("Omicron BA.1", "Omicron BA.2")),
    case_type = ifelse(case_type == "c", "Community", "Imported")
  ) %>% ungroup() %>% left_join(
    df_vline %>% transmute(variant, date_intro=date),
    by="variant"
  ) %>% mutate(
    lag = as.numeric(date - date_intro)
  )
  df_sims
}

# # some of the units have no cases detected or sequenced, so the EDT is Inf, we need to calculate the proportion of units that have cases detected or sequenced
# summarise_EDT_props <- function(
#   data_EDT
# ){
#   data_EDT %>% summarise(across(.fns = ~ sum(!is.infinite(.x)) / n())) %>% select(-date_decimal)
# }

# summarise_EDT_quants <- function(
#   data_EDT,
#   quants_range = c(0.5)
# ){
#   data_EDT %>% select(-id_sim) %>% group_by(unitname) %>% reframe(across(.fns = ~ quantile(.x, probs = quants_range))) %>% select(-date_decimal)
# }

