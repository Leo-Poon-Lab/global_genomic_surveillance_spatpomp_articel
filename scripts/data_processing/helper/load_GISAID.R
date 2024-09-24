library(tidyverse)

## function for loading GISAID data
load_GISAID <- function(date_start, date_end, country_under_investigation, df_new_code){
    data_meta <- read_csv(here::here("data/gisaid_data/df_meta_lin_sum.csv"))
    # df_top_lineange <- read_csv(here::here("data/gisaid_data/df_top_lineange.csv"))
    if("variant" %in% colnames(data_meta)){
        data_meta$lineage <- data_meta$variant
    }

    data_meta_study_period <- data_meta %>% filter(`Collection date`>=date_start & `Collection date`<=date_end)

    data_meta_study_period$parent_country <- factor(data_meta_study_period$parent_country, levels=df_new_code$code, labels=df_new_code$new_code)
    data_meta_study_period$sequencing_country <- factor(data_meta_study_period$sequencing_country, levels=df_new_code$code, labels=df_new_code$new_code)

    data_meta_study_period <- data_meta_study_period %>% ungroup() %>% group_by(`Collection date`, date_decimal, sequencing_country, parent_country, lineage) %>% summarise(N=sum(N)) %>% ungroup()

    data_meta_study_period <- data_meta_study_period %>% filter(sequencing_country %in% country_under_investigation) %>% filter(parent_country %in% country_under_investigation)

    data_meta_study_period_aggregated <- data_meta_study_period %>% mutate(code=sequencing_country) %>% group_by(`Collection date`, date_decimal, code, lineage) %>% summarise(N=sum(N)) %>% ungroup()

    list(
      data_GISAID=data_meta_study_period_aggregated,
      data_GISAID_origin=data_meta_study_period,
      lineages=sort(table(data_meta_study_period$lineage),decreasing = T)
    )

}
