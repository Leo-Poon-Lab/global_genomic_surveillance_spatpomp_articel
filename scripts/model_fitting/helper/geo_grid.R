geo_grid_29 <- function() {
  cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% dplyr::filter(level==3)
  data_measurements <- readRDS(paste0("results/model_data/data_measurements_", "Omicron20", ".rds")) # for reported cases, imported/community sequenced cases
  regions_code_3 <- cross_check_table %>% filter(code %in% unique(data_measurements$code)) %>% .$country_iso3c
  other_codes <- unique(data_measurements$code)[grepl("_others", unique(data_measurements$code))]
  geo_grid_Omicron_20 <- data.frame(
    name = c("Canada", "United Kingdom", "EU_others", "Germany", "Japan", "United States of America", "France", "South Korea", "AS_others", "Mexico", "NA_others", "Spain", "Hong Kong SAR", "Turkey", "United Arab Emirates", "India", "Tunisia", "Brazil", "SA_others", "Chile", "Nigeria", "Papua New Guinea", "Morocco", "OC_others", "AF_others", "Argentina", "South Africa", "New Zealand", "Australia"),
    code = c("CAN", "GBR", "EU_others", "DEU", "JPN", "USA", "FRA", "KOR", "AS_others", "MEX", "NA_others", "ESP", "HKG", "TUR", "ARE", "IND", "TUN", "BRA", "SA_others", "CHL", "NGA", "PNG", "MAR", "OC_others", "AF_others", "ARG", "ZAF", "NZL", "AUS"),
    # code_country = c("124", "826", "NA", "276", "392", "840", "250", "410", "NA", "484", "NA", "724", "NA", "792", "784", "356", "504", "076", "NA", "152", "788", "598", "566", "NA", "NA", "032", "710", "554", "036"),
    code_iso_3166_2 = c("ISO 3166-2:CA", "ISO 3166-2:GB", "NA", "ISO 3166-2:DE", "ISO 3166-2:JP", "ISO 3166-2:US", "ISO 3166-2:FR", "ISO 3166-2:KR", "NA", "ISO 3166-2:MX", "NA", "ISO 3166-2:ES", "NA", "ISO 3166-2:TR", "ISO 3166-2:AE", "ISO 3166-2:IN", "ISO 3166-2:TN", "ISO 3166-2:BR", "NA", "ISO 3166-2:CL", "ISO 3166-2:NG", "ISO 3166-2:PG", "ISO 3166-2:MA", "NA", "NA", "ISO 3166-2:AR", "ISO 3166-2:ZA", "ISO 3166-2:NZ", "ISO 3166-2:AU"),
    col = c(1, 2, 3, 4, 6, 1, 3, 6, 5, 1, 2, 3, 6, 4, 5, 5, 4, 2, 2, 1, 4, 6, 3, 5, 3, 1, 4, 6, 5),
    row = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
    stringsAsFactors = FALSE
  )
  geo_grid_Omicron_20$code_3 <- as.character(geo_grid_Omicron_20$code)
  geo_grid_Omicron_20$code <- sapply(geo_grid_Omicron_20$code, function(x){
    if(x %in% regions_code_3){
      return(cross_check_table %>% filter(country_iso3c==x) %>% .$code)
    } else if(x %in% other_codes){
      return(x)
    } else{
      return(NA)
    }
  })
  return(geo_grid_Omicron_20)
}
