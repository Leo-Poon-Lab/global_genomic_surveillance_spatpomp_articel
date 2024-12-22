get_mov_mat_aggregate_2019 <- function(){
  if(file.exists("data/estimated_movement_matrix/move_mat_2019_country_aggregate.rds")){
    mov_mat_aggregate <- readRDS("data/estimated_movement_matrix/move_mat_2019_country_aggregate.rds")
    return(mov_mat_aggregate)
  } else{
    ## Build specific movement matrix for M_4, using only 2019 moving data##
    year_start <- 2019 # use the 2019 data to mimic the non-pandemic situation
    year_end <- 2019 # use the 2019 data to mimic the non-pandemic situation
    global_scale <- "country"
    mov_mat <- lapply(unique(c(year_start, year_end)), function(this_year) {
      tmp <- readRDS(here::here(paste0(
        "data/estimated_movement_matrix/move_mat_",
        this_year,
        "_",
        global_scale,
        ".rds")))
      tmp
    })
    mov_mat <- bind_rows(mov_mat)
    mov_mat <- mov_mat %>% group_by(code_o, code_d) %>% summarise(flow_all=mean(flow_all)) %>% ungroup() # average over 2019

    df_new_code <- readRDS("results/model_data/data_fitting_Omicron20.rds") %>% .$df_new_code
    mov_mat_aggregate <- mov_mat %>% left_join(df_new_code %>% transmute(code_o=code, new_code_o=new_code), by= "code_o") %>% left_join(df_new_code %>% transmute(code_d=code, new_code_d=new_code), by= "code_d") %>% mutate(code_o=new_code_o, code_d=new_code_d) %>% select(code_o, code_d, flow_all) %>% group_by(code_o, code_d) %>% summarise(flow_all=sum(flow_all)) %>% ungroup() %>% arrange(code_o, code_d) %>% filter(code_o!=code_d)

    saveRDS(mov_mat_aggregate, "data/estimated_movement_matrix/move_mat_2019_country_aggregate.rds")
    return(mov_mat_aggregate)
  }
}