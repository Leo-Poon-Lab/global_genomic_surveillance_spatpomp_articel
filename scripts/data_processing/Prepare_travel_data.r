library(tidyverse)
library(lubridate)
library(geosphere)
# library(rgeos)
library(maps)
library(rnaturalearth)
library(WDI)
library(gravity)
library(ggplot2)
library(ggrepel)

setwd(here::here())

# 1. load covid-19 airport data (ADS-B, Automatic Dependent Surveillance-Broadcast) #

files_asdb_data <- list.files(here::here("data/covid19-flight_data/"), ".gz$", full.names = T)
df_asdb <- lapply(files_asdb_data, function(x){
	# x=files_data[1]
	df_tmp <- read_csv(x, col_types=cols(.default = "c"), na="")
	# df_tmp %>% filter(!is.na(registration))
	df_tmp %>% select(
    callsign, # the identifier of the flight displayed on ATC screens (usually the first three letters are reserved for an airline: AFR for Air France, DLH for Lufthansa, etc.)
    icao24, # the transponder unique identification number
    registration, # the aircraft tail number (when available)
    typecode, # the aircraft model type (when available)
    origin, # a four letter code for the origin airport of the flight (when available)
    destination, # a four letter code for the destination airport of the flight (when available)
    day # the UTC day of the last message received by the OpenSky Network
  )
})
df_asdb <- bind_rows(df_asdb) # nrow: 123,159,149
# format(object.size(df_asdb), units = "auto")

# 2. load aircraft data #
# In the context of ICAO (International Civil Aviation Organization) Aircraft Type Designators, "WTC" stands for Wake Turbulence Category. The Wake Turbulence Category is a classification system that categorizes aircraft based on their maximum certificated takeoff mass. This classification is important for air traffic controllers to manage aircraft separation, ensuring that smaller aircraft avoid the turbulent wake created by larger ones. #

# There are three main categories of WTC: #
# Light (L): Aircraft with a maximum takeoff mass of 34,000 lb (15,500 kg) or less #
# Small (S): Aircraft with a maximum takeoff mass between 34,000 lb and 41,000 lb (15,500 kg and 18,600 kg) #
# Medium (M): Aircraft with a maximum takeoff mass between 41,000 lb and 299,000 lb (18,600 kg and 136,000 kg) #
# Heavy (H): Aircraft with a maximum takeoff mass exceeding 300,000 lb (136,000 kg) #
# Super (J): Currently, the Airbus A380-800 is the only aircraft designated as Super category. #

# It is important to note that the Wake Turbulence Category (WTC) is primarily based on an aircraft's maximum certificated takeoff mass, not on the number of passengers it can carry. However, the categories can provide a rough idea of the average number of passengers for each group. #
# Light (L): Light aircraft typically have a seating capacity of 1 to 8 passengers. Examples of light aircraft include popular general aviation planes like the Cessna 172 Skyhawk and Beechcraft Bonanza. #
# Small (S): Small aircraft may accommodate up to 19 passengers. Examples in this category include turboprops like the Beechcraft King Air 350 and small regional jets like the Embraer EMB-120 Brasilia. #
# Medium (M): Medium aircraft can accommodate between 20 and 250 passengers. Examples of medium aircraft include regional jets like the Bombardier CRJ900 and smaller narrow-body commercial jets like the Airbus A220 or the Boeing 737. #
# Heavy (H): Heavy aircraft are typically larger and have a passenger capacity between 250 and 550. Examples include wide-body airliners like the Boeing 777 and the Airbus A330. #
# Super (J): The Airbus A380-800, the only aircraft currently designated as Super category, has a seating capacity of around 853 passengers (868 in a typical two-class configuration or up to 853 in high-density configurations). #
# Please note that these numbers are approximate averages, as the actual number of passengers can vary between different aircraft configurations and airlines. #

# df_aircraft <- read_csv(here::here("data/open_sky_data/aircraft-database-complete-2023-02.csv", col_types=cols(.default = "c"))
df_aircraft_type <- read_csv(here::here("data/open_sky_data/doc8643AircraftTypes.csv"), col_types=cols(.default = "c"), na="")
df_aircraft_type <- df_aircraft_type %>% select(Designator, WTC) %>% unique()
names(df_aircraft_type)[1] <- "typecode"
df_aircraft_type <- df_aircraft_type %>% group_by(typecode) %>% filter(row_number()==1)
# df_aircraft <- left_join(df_aircraft, df_aircraft_type, "typecode")

sum(is.na(df_asdb$typecode))/nrow(df_asdb) # about 37% of the ADS-B flights don't have the aircraft model type
tmp <- df_asdb$typecode[!is.na(df_asdb$typecode)]
sum(!tmp %in% df_aircraft_type$typecode)/length(tmp)*100 # about 1.5% of the ADS-B flights don't have WTC aircraft type data
# unique(tmp[!tmp %in% df_aircraft_type$typecode])
df_asdb <- left_join(df_asdb, df_aircraft_type, "typecode")

df_asdb$WTC[is.na(df_asdb$WTC)] <- "S" # all the flights with missing WTC type, will be considered small flights
round(table(df_asdb$WTC)/nrow(df_asdb),2)*100
df_asdb$estimated_num_of_passenger <- NA
# We will approximately estimate the number of passengers in each flight, as 80% of the max capacity #
# On average, the passenger load factor (PLF) for major airlines tends to be around 75% to 85%. But we should consider there is some cargo flights or other flight not intended for carrying passengers. Here we assume 80% of the flights are for carrying passengers.
frac_capacity=0.8
frac_passenger_aircraf=0.8
df_asdb$estimated_num_of_passenger[df_asdb$WTC=="H"] <- frac_passenger_aircraf*frac_capacity*550
df_asdb$estimated_num_of_passenger[df_asdb$WTC=="M"] <- frac_passenger_aircraf*frac_capacity*(250+550)/2
df_asdb$estimated_num_of_passenger[df_asdb$WTC=="L/M"] <- frac_passenger_aircraf*frac_capacity*250
df_asdb$estimated_num_of_passenger[df_asdb$WTC=="L"] <- frac_passenger_aircraf*frac_capacity*19
df_asdb$estimated_num_of_passenger[df_asdb$WTC=="S"] <- frac_passenger_aircraf*frac_capacity*8

# 3. load airports data and find country and city of the flights #
df_airports <- read_csv(here::here("data/our_airports_data/airports.csv"), col_types=cols(.default = "c"), na="")
df_countries <- read_csv(here::here("data/our_airports_data/countries.csv"), col_types=cols(.default = "c"), na="")
df_regions <- read_csv(here::here("data/our_airports_data/regions.csv"), col_types=cols(.default = "c"), na="")

df_airports <- df_airports %>% select(ident:name, iso_country:municipality)
df_airports <- left_join(df_airports, df_countries %>% mutate(iso_country=code, country=name) %>% select(iso_country, country, continent))
df_airports <- left_join(df_airports, df_regions %>% mutate(iso_region=code, region=name) %>% select(iso_region, region))

sum(is.na(df_asdb$destination))/nrow(df_asdb)*100 # about 18% of the ADS-B flights don't have valid destination
tmp <- df_asdb$destination[!is.na(df_asdb$destination)]
sum(!tmp %in% df_airports$ident)/length(tmp)*100 # about 0.3% of the destination airports are not identifiable.
df_asdb_origin <- left_join(df_asdb %>% select(origin), df_airports %>% mutate(origin=ident) %>% select(-ident), "origin")
names(df_asdb_origin) <- paste0("ori_", names(df_asdb_origin))
df_asdb <- bind_cols(df_asdb, df_asdb_origin[,-1])

df_asdb_dest <- left_join(df_asdb %>% select(destination), df_airports %>% mutate(destination=ident) %>% select(-ident), "destination")
names(df_asdb_dest) <- paste0("dest_", names(df_asdb_dest))
df_asdb <- bind_cols(df_asdb, df_asdb_dest[,-1])

# 4. data cleaning: filter flights to be used in the model #

# 4.1 we only include flights with clear origin and destinations #
(sum(!is.na(df_asdb$destination) & !is.na(df_asdb$origin))/nrow(df_asdb))*100 # about 62% of the ADS-B flights have complete origin and destination info.
df_asdb <- df_asdb %>% filter(!is.na(df_asdb$destination) & !is.na(df_asdb$origin))
(sum(df_asdb$origin != df_asdb$destination, na.rm=T)/nrow(df_asdb))*100 # about 67% of the ADS-B flights traveled between different airports.
(sum(df_asdb$ori_iso_region != df_asdb$dest_iso_region, na.rm=T)/nrow(df_asdb))*100 # about 67% of the ADS-B flights traveled between different origin and destination regions.
df_asdb <- df_asdb %>% filter(ori_iso_region != dest_iso_region)

# 5. summaries the movement matrix, by day, by region #
# 50,524,872 entries between different locations # 
df_asdb$date <- ymd(ymd_hms(df_asdb$day))
saveRDS(df_asdb, here::here("data/covid19-flight_data/df_asdb.rds"))
df_asdb <- readRDS(here::here("data/covid19-flight_data/df_asdb.rds"))

for (this_level in c(3, 2)){
  # this_level=3
  # 5.1 reformat the per-flight data to a summarized table which uses per-day per-region as entry #
  stopifnot(!any(df_asdb$ori_iso_region==df_asdb$dest_iso_region))

  ## per level ##
  if(this_level==2){ # per continent #
    df_asdb_per_day <- df_asdb %>% filter(ori_continent!=dest_continent) %>% group_by(date, ori_continent, dest_continent) %>% summarise(N_flight=n(), N_passenger=sum(estimated_num_of_passenger))
    this_level_name <- "continent"
  } else if(this_level==3){ # per country #
    df_asdb_per_day <- df_asdb %>% filter(ori_country!=dest_country) %>% group_by(date, ori_continent, ori_iso_country, dest_continent, dest_iso_country) %>% summarise(N_flight=n(), N_passenger=sum(estimated_num_of_passenger))
    this_level_name <- "country"
  }  else if(this_level==4){ # per region #
    df_asdb_per_day <- df_asdb %>% group_by(date, ori_continent, ori_iso_country, ori_iso_region, dest_continent, dest_iso_country, dest_iso_region) %>% summarise(N_flight=n(), N_passenger=sum(estimated_num_of_passenger))
    this_level_name <- "region"
  } else {
    stop()
  }
  
  # compute yearly flow to compare with the WTO and ICAO data #
  df_asdb_per_day$code <- df_asdb_per_day[[paste0("ori_iso_", this_level_name)]]
  df_asdb_per_day_yearly_outflow <- df_asdb_per_day %>% mutate(year=year(date)) %>% group_by(year, code) %>% summarise(N_outflow_flight=sum(N_flight), N_outflow_passenger=sum(N_passenger))
  df_asdb_per_day$code <- df_asdb_per_day[[paste0("dest_iso_", this_level_name)]]
  df_asdb_per_day_yearly_inflow <- df_asdb_per_day %>% mutate(year=year(date)) %>% group_by(year, code) %>% summarise(N_inflow_flight=sum(N_flight), N_inflow_passenger=sum(N_passenger))
  df_asdb_per_day <- df_asdb_per_day %>% select(-code)

  df_asdb_per_day_yearly_flow <- left_join(df_asdb_per_day_yearly_inflow, df_asdb_per_day_yearly_outflow)
  df_asdb_per_day_yearly_flow$N_flow_passenger <- df_asdb_per_day_yearly_flow$N_inflow_passenger + df_asdb_per_day_yearly_flow$N_outflow_passenger
  df_asdb_per_day_yearly_flow <- df_asdb_per_day_yearly_flow %>% filter(!is.na(N_outflow_passenger))

  if(this_level==3){ # compete the missing data at level 3
    # 5.2 adjust the ADS-B coverage bias by using the international travel data from 1. World Tourism Organization, Yearbook of Tourism Statistics, Compendium of Tourism Statistics; 2. International Civil Aviation Organization #
    cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx"))
    df_ihme <- readRDS(here::here("data/IHME/df_ihme.rds"))
    df_ihme_pop <- df_ihme %>% mutate(year=year(date)) %>% filter(level==this_level) %>% select(year,population, code, level) %>% unique() %>% group_by(year, code, level) %>% summarise(population_mean=mean(population)) %>% ungroup()
      
    # 5.2.0 Validate the correlation between the data from WTO and ICAO #
    metadata_country_WTO <- read_csv(here::here("data/The world bank data/World Tourism Organization data/International tourism_number of arrivals/Metadata_Country_API_ST.INT.ARVL_DS2_en_csv_v2_5728898.csv"))
    if(!all(cross_check_table$country_iso3c %in% metadata_country_WTO$`Country Code`)){
      cross_check_table %>% filter(!country_iso3c %in% metadata_country_WTO$`Country Code`) # no data from Taiwan, others are good
      countries_of_interest <- metadata_country_WTO$`Country Code`[metadata_country_WTO$`Country Code` %in% cross_check_table$country_iso3c]
    }

    df_WTO_arrival <- read_csv(here::here("data/The world bank data/World Tourism Organization data/International tourism_number of arrivals/API_ST.INT.ARVL_DS2_en_csv_v2_5728898.csv"))
    df_WTO_arrival <- df_WTO_arrival %>% filter(`Country Code` %in% countries_of_interest)
    df_WTO_departure <- read_csv(here::here("data/The world bank data/World Tourism Organization data/International tourism_number of departures/API_ST.INT.DPRT_DS2_en_csv_v2_5734520.csv"))
    df_WTO_departure <- df_WTO_departure %>% filter(`Country Code` %in% countries_of_interest)
    stopifnot(all(df_WTO_departure$`Country Code` == df_WTO_arrival$`Country Code`))

    df_ICAO_passenger <- read_csv(here::here("data/The world bank data/International Civil Aviation Organization data/Air transport_passengers carried/API_IS.AIR.PSGR_DS2_en_csv_v2_5728833.csv"))
    df_ICAO_passenger <- df_ICAO_passenger %>% filter(`Country Code` %in% countries_of_interest)
    stopifnot(all(df_ICAO_passenger$`Country Code` == df_WTO_arrival$`Country Code`))

    df_ICAO_carrier <- read_csv(here::here("data/The world bank data/International Civil Aviation Organization data/registered carrier departures worldwide/API_IS.AIR.DPRT_DS2_en_csv_v2_5735762.csv"))
    df_ICAO_carrier <- df_ICAO_carrier %>% filter(`Country Code` %in% countries_of_interest)
    stopifnot(all(df_ICAO_carrier$`Country Code` == df_WTO_arrival$`Country Code`))

    # take 2019 data as an example #
    df_WTO_arrival_2019 <- df_WTO_arrival %>% select(`Country Code`, `2019`)
    df_WTO_departure_2019 <- df_WTO_departure %>% select(`Country Code`, `2019`)
    df_WTO_total_2019 <- left_join(df_WTO_arrival_2019, df_WTO_departure_2019, "Country Code")
    df_WTO_total_2019$total_flow <- df_WTO_total_2019$`2019.x` + df_WTO_total_2019$`2019.y`

    check <- (!is.na(df_ICAO_passenger$`2019`) & !is.na(df_WTO_total_2019$total_flow))
    plot(df_ICAO_passenger$`2019`[check], df_WTO_total_2019$total_flow[check])
    cor(df_ICAO_passenger$`2019`[check], df_WTO_total_2019$total_flow[check]) # 0.7538952, strong positive correlation
    check <- (!is.na(df_ICAO_carrier$`2019`) & !is.na(df_WTO_total_2019$total_flow))
    plot(df_ICAO_carrier$`2019`[check], df_WTO_total_2019$total_flow[check])
    cor(df_ICAO_carrier$`2019`[check], df_WTO_total_2019$total_flow[check]) # 0.7226753, strong positive correlation

    # take 2018 data as an example #
    df_WTO_arrival_2018 <- df_WTO_arrival %>% select(`Country Code`, `2018`)
    df_WTO_departure_2018 <- df_WTO_departure %>% select(`Country Code`, `2018`)
    df_WTO_total_2018 <- left_join(df_WTO_arrival_2018, df_WTO_departure_2018, "Country Code")
    df_WTO_total_2018$total_flow <- df_WTO_total_2018$`2018.x` + df_WTO_total_2018$`2018.y`

    check <- (!is.na(df_ICAO_passenger$`2018`) & !is.na(df_WTO_total_2018$total_flow))
    plot(df_ICAO_passenger$`2018`[check], df_WTO_total_2018$total_flow[check])
    cor(df_ICAO_passenger$`2018`[check], df_WTO_total_2018$total_flow[check]) # 0.7421278, strong positive correlation
    check <- (!is.na(df_ICAO_carrier$`2018`) & !is.na(df_WTO_total_2018$total_flow))
    plot(df_ICAO_carrier$`2018`[check], df_WTO_total_2018$total_flow[check])
    cor(df_ICAO_carrier$`2018`[check], df_WTO_total_2018$total_flow[check]) # 0.712501, strong positive correlation

    # The ratio between carrier and passenger # 
    quantile(df_ICAO_passenger$`2019`/df_ICAO_carrier$`2019`, na.rm=T)
    #       0%        25%        50%        75%       100% 
    # 5.256278  66.191279  99.217525 126.387596 220.564600 
    quantile(df_WTO_total_2019$total_flow/df_ICAO_carrier$`2019`, na.rm=T)
    #        0%          25%          50%          75%            100% 
    # 20.98390     78.52497    165.66945    613.25608    354772.72727
    quantile(df_asdb_per_day$N_passenger/df_asdb_per_day$N_flight, na.rm = T)
    #       0%      25%      50%      75%     100% 
    #   5.1200 130.5600 233.8743 265.2800 352.0000 

    # 5.2.1 we try to link the movement statistics between the ADS-B and census data #
    # At the end, we chose to use the ICAO data (total flow with gravity model) rather than the WTO data (inflow and outflow with pairwise interaction model), because WTO data have too much missing data after 2019 #

    # 5.2.1.1 in the US and Europe, are the two datasets (ADS-B and yearly census) consistent? is there any offsets to be set? #
    # take 2019 data as an example #

    df_flow_comp_wto_2019 <- left_join(left_join(df_asdb_per_day_yearly_flow %>% filter(year==2019), cross_check_table %>% filter(level==3), c("code")) %>% select(country_iso3c, code, N_flow_passenger, continent), df_WTO_total_2019 %>% mutate(country_iso3c=`Country Code`, flow_wto=total_flow) %>% select(country_iso3c, flow_wto))
    ggpubr::ggscatter(df_flow_comp_wto_2019, x = "N_flow_passenger", y = "flow_wto", # the correlation is linear #
      add = "reg.line",  # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
      conf.int = TRUE # Add confidence interval
    )+
      xlab("Num of passengers (ADS-B)")+
      ylab("Num of passengers (WTO)")+
      geom_text_repel(aes(label=code), data=df_flow_comp_wto_2019 %>% filter(abs(flow_wto-N_flow_passenger)>8*10^7))+
      ggpubr::stat_cor(method = "pearson", label.x = 3, label.y = 350000000)+
      geom_abline(intercept = 0, slope = 1)
    ggsave(here::here("results/figs/fig_total_flow_adsb_wto_2019_all.pdf"), width=8, height=5)
    # The data between WTO and ADS-B is highly correlated #
    # CN ADS-B data should be replaced by the gravity model estimation #
    ggpubr::ggscatter(df_flow_comp_wto_2019 %>% filter(continent %in% c("EU", "NA")), x = "N_flow_passenger", y = "flow_wto", # the correlation is linear #
      add = "reg.line",  # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
      conf.int = TRUE # Add confidence interval
    )+
      xlab("Num of passengers (ADS-B)")+
      ylab("Num of passengers (WTO)")+
      geom_text_repel(aes(label=code), data=. %>% filter(abs(flow_wto-N_flow_passenger)>8*10^7))+
      ggpubr::stat_cor(method = "pearson", label.x = 3, label.y = 350000000)+
      geom_abline(intercept = 0, slope = 1)
    ggsave(here::here("results/figs/fig_total_flow_adsb_wto_2019_eu_na.pdf"), width=8, height=5)

    df_flow_comp_icao_2019 <- left_join(left_join(df_asdb_per_day_yearly_flow %>% filter(year==2019), cross_check_table %>% filter(level==3), c("code")) %>% select(country_iso3c, code, continent, N_flow_passenger), df_ICAO_passenger %>% select(`Country Code`, `2019`) %>% mutate(country_iso3c=`Country Code`, flow_icao=`2019`) %>% select(country_iso3c, flow_icao))
    ggpubr::ggscatter(df_flow_comp_icao_2019, x = "N_flow_passenger", y = "flow_icao", # the correlation is linear #
      add = "reg.line",  # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
      conf.int = TRUE # Add confidence interval
    )+
      xlab("Num of passengers (ADS-B)")+
      ylab("Num of passengers (ICAO)")+
      geom_text_repel(aes(label=code), data=. %>% filter(abs(flow_icao-N_flow_passenger)>2*10^8))+
      ggpubr::stat_cor(method = "pearson", label.x = 3, label.y = 300000000)+
      geom_abline(intercept = 0, slope = 1)
    ggsave(here::here("results/figs/fig_total_flow_adsb_icao_2019_all.pdf"), width=8, height=5)
    
    ggpubr::ggscatter(df_flow_comp_icao_2019 %>% filter(continent %in% c("EU", "NA")), x = "N_flow_passenger", y = "flow_icao", # the correlation is linear #
      add = "reg.line",  # Add regressin line
      add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
      conf.int = TRUE # Add confidence interval
    )+
      xlab("Num of passengers (ADS-B)")+
      ylab("Num of passengers (ICAO)")+
      geom_text_repel(aes(label=code), data=. %>% filter(abs(flow_icao-N_flow_passenger)>2*10^8))+
      ggpubr::stat_cor(method = "pearson", label.x = 3, label.y = 300000000)+
      geom_abline(intercept = 0, slope = 1)
    ggsave(here::here("results/figs/fig_total_flow_adsb_icao_2019_eu_na.pdf"), width=8, height=5)
    # Note: China is the most significant outlier. USA also have higher ADS-B estimates, this is partially due to US have high number of domestic flights (ICAO includes both international and domestic flights) #
    # in the following analysis, we should use gravity model data for China #
    
    # As the international flights approximately contribute to 36.3% of the total flights, we should adjust the ICAO number for US to 36.3% #
    df_ICAO_passenger_2019 <- df_ICAO_passenger %>% select(`Country Code`, `2019`)
    left_join(df_asdb_per_day_yearly_flow, cross_check_table %>% filter(level==3), c("code")) %>% filter(year==2019) %>% filter(country_iso3c=="USA") %>% t()
    df_ICAO_passenger_2019 %>% filter(`Country Code`=="USA")
    df_WTO_total_2019 %>% filter(`Country Code`=="USA")
    (usa_adjust_ratio <- df_WTO_total_2019$`total_flow`[df_WTO_total_2019$`Country Code`=="USA"]/df_ICAO_passenger_2019$`2019`[df_ICAO_passenger_2019$`Country Code`=="USA"])

    left_join(df_asdb_per_day_yearly_flow, cross_check_table %>% filter(level==3), c("code")) %>% filter(year==2019) %>% filter(country_iso3c=="ZAF") %>% t()
    df_ICAO_passenger_2019 %>% filter(`Country Code`=="ZAF")
    df_WTO_total_2019 %>% filter(`Country Code`=="ZAF")

    # 5.2.1.2 We use the ICAO data to predict the actual flow (represented by ADS-B data, excluding outlier China), using a gravity model, the gravity model also considered the GDP, urban population of each country, and distance between countries #
    # In the gravity model, please refer to notes_gravity_model.md #
    years_all <- 2019:2022 
    sapply(years_all, function(this_year) {
      # this_year=2022
      # Step 0: prepare data #
      # flow data #
      
      # ICAO data for 2022 is missing
      ## use 2021 data for 2022 onwards
      if (this_year==2022){
        this_year_ICAO <- 2021
      } else{
        this_year_ICAO <- this_year
      }

      df_ICAO_passenger_i <- df_ICAO_passenger %>% select(any_of(c("Country Code", this_year_ICAO)))
      names(df_ICAO_passenger_i)[names(df_ICAO_passenger_i)==this_year_ICAO] <- "ICAO_flow"
      
      df_flow_comp_icao_i <- left_join(left_join(df_asdb_per_day_yearly_flow %>% filter(year==this_year), cross_check_table %>% filter(level==3), c("code")) %>% select(country_iso3c, code, continent, N_inflow_passenger, N_outflow_passenger,N_flow_passenger), df_ICAO_passenger_i %>% mutate(country_iso3c=`Country Code`))
      df_flow_comp_icao_i <- df_flow_comp_icao_i %>% filter(country_iso3c != "CHN") %>% filter(!is.na(ICAO_flow)) %>% filter(!is.na(N_flow_passenger))
      num_of_model_countries_i <- nrow(df_flow_comp_icao_i)

      panel_data <- expand.grid(origin = 1:num_of_model_countries_i,
                    destination = 1:num_of_model_countries_i,
                    stringsAsFactors = FALSE)
      panel_data$code_o <- df_flow_comp_icao_i$code[panel_data$origin]
      panel_data$code_d <- df_flow_comp_icao_i$code[panel_data$destination]
      panel_data$code_3_o <- df_flow_comp_icao_i$country_iso3c[panel_data$origin]
      panel_data$code_3_d <- df_flow_comp_icao_i$country_iso3c[panel_data$destination]
      panel_data <- panel_data %>% filter(origin!=destination)

      # Fetch GDP data (NY.GDP.MKTP.CD indicates GDP per capita in current USD
      gdp_pc_data <- WDI(indicator = "NY.GDP.PCAP.CD", start = this_year, end = this_year)
      # Fetch urban population data (https://data.worldbank.org/indicator/SP.URB.TOTL?view=chart
      urban_pop_data <- WDI(indicator = "SP.URB.TOTL", start = this_year, end = this_year)
      # location data #
      spdf_world <- ne_countries()
      # spdf_world <- spdf_world@data

      # ADS-B data
      df_asdb_pairflow_this_year <- df_asdb_per_day %>% ungroup() %>% mutate(year=year(date)) %>% filter(year==this_year) 
      df_asdb_pairflow_year <- df_asdb_pairflow_this_year %>% group_by(year, ori_iso_country, dest_iso_country) %>% summarise(N_flight=sum(N_flight), N_passenger=sum(N_passenger))

      panel_data_extra <- apply(panel_data, 1, function(x){
        # x=pairs_i[,1]
        code_o <- x[3]
        code_d <- x[4]
        code_o_3 <- x[5]
        code_d_3 <- x[6]
        
        # Get geospatial data
        point1 <- spdf_world[which(spdf_world$adm0_a3 == code_o_3), ]
        point2 <- spdf_world[which(spdf_world$adm0_a3 == code_d_3), ]

        # Extract latitude and longitude
        point1 <- c(point1$label_x[1], point1$label_y[1])
        point2 <- c(point2$label_x[1], point2$label_y[1])

        # Calculate distance (in meters)
        # distance_meters <- distGeo(point1, point2)
        dist <- distVincentySphere(point1, point2)

        # GDP
        gdp_data_o <- gdp_pc_data$NY.GDP.PCAP.CD[gdp_pc_data$iso3c==code_o_3]
        gdp_data_d <- gdp_pc_data$NY.GDP.PCAP.CD[gdp_pc_data$iso3c==code_d_3]

        # Urban population
        pop_data_o <- urban_pop_data$SP.URB.TOTL[urban_pop_data$iso3c==code_o_3]
        pop_data_d <- urban_pop_data$SP.URB.TOTL[urban_pop_data$iso3c==code_d_3]

        # ICAO
        icao_flow_o <- df_flow_comp_icao_i$ICAO_flow[df_flow_comp_icao_i$code==code_o]
        icao_flow_d <- df_flow_comp_icao_i$ICAO_flow[df_flow_comp_icao_i$code==code_d]
        if(code_o_3=="USA"){icao_flow_o <- icao_flow_o*usa_adjust_ratio}
        if(code_d_3=="USA"){icao_flow_d <- icao_flow_o*usa_adjust_ratio}

        # ADS-B flow
        adsb_flow <- df_asdb_pairflow_year$N_passenger[df_asdb_pairflow_year$ori_iso_country==code_o & df_asdb_pairflow_year$dest_iso_country==code_d][1]

        tibble(code_o=code_o, code_d=code_d, code_3_o=code_o_3, code_3_d=code_d_3, dist=dist, gdp_o=gdp_data_o, gdp_d=gdp_data_d, upop_o=pop_data_o, upop_d=pop_data_d, icao_flow_o=icao_flow_o, icao_flow_d=icao_flow_d, adsb_flow=adsb_flow)
      })
      panel_data_extra <- bind_rows(panel_data_extra)
      panel_data_completed <- as_tibble(left_join(panel_data, panel_data_extra))

      # Fit the gravity model using Poisson regression #
      gravity_model_flow <- glm(adsb_flow ~ log(icao_flow_o) + log(icao_flow_d) + log(gdp_o) + log(gdp_d) + log(upop_o) + log(upop_d) + log(dist), data = panel_data_completed, family = poisson(link = "log"))
      summary(gravity_model_flow)
      qplot(panel_data_completed$adsb_flow, predict(gravity_model_flow, panel_data_completed, type = "response"))
      (rmse_flow <- sqrt(mean((panel_data_completed$adsb_flow - exp(predict(gravity_model_flow, panel_data_completed)))^2, na.rm=T)))
      
      # gravity_model_2_inflow <- glm(adsb_inflow ~ log(gdp_o) + log(gdp_d) + log(upop_o) + log(upop_d) + log(dist), data = panel_data_completed, family = poisson(link = "log"))
      

      # 5.2.2 we calculate the expected annual movement between nodes using the gravity (GLM) model, predicting the flow for China and other countries with no ADS-B data #
      # we check which countries lacks the ADS-B data # 
      countries_with_asdb_data <- df_asdb_per_day_yearly_flow %>% filter(year==this_year) %>% filter(!is.na(N_flow_passenger)) %>% .$code
      countries_to_fill <- c("CN", unique(df_ihme_pop$code)[!unique(df_ihme_pop$code) %in% countries_with_asdb_data])
      countries_to_fill_code_3 <- sapply(countries_to_fill, function(x) {cross_check_table %>% filter(level==3) %>% filter(code==x) %>% .$country_iso3c})
      check <- countries_to_fill_code_3 %in% (df_ICAO_passenger_i %>% filter(!is.na(ICAO_flow)) %>% .$`Country Code`)
      countries_to_fill_can_fill <- countries_to_fill[check]
      countries_to_fill_can_fill_3 <- sapply(countries_to_fill_can_fill, function(x) {cross_check_table %>% filter(level==3) %>% filter(code==x) %>% .$country_iso3c})
      countries_to_fill_cant_fill <- countries_to_fill[!check]
      (num_of_countries_can_fill <- length(countries_to_fill_can_fill)) # 41 countries

      countries_no_need_to_fill <- unique(df_ihme_pop$code)[!unique(df_ihme_pop$code) %in% countries_to_fill]
      countries_no_need_to_fill_3 <- sapply(countries_no_need_to_fill, function(x) {cross_check_table %>% filter(level==3) %>% filter(code==x) %>% .$country_iso3c})
      check <- countries_no_need_to_fill_3 %in% (df_ICAO_passenger_i %>% filter(is.na(ICAO_flow)) %>% .$`Country Code`)
      countries_no_need_to_fill <- countries_no_need_to_fill[!check]
      countries_no_need_to_fill_3 <- countries_no_need_to_fill_3[!check]
      (num_of_countries_no_need_to_fill <- length(countries_no_need_to_fill)) # 97 countries

      cross_check_table %>% filter(level==3) %>% filter(code %in% countries_to_fill_can_fill) %>% .$continent %>% table()
      # AF AS EU NA OC SA 
      # 22  7  2  5  2  3 
      # Most are from Africa 
      cross_check_table %>% filter(level==3) %>% filter(code %in% countries_to_fill_cant_fill) %>% .$loc_name
      # we lack the data for these countries, so we don't include them in the current study

      countries_analysis <- c(countries_to_fill_can_fill, countries_no_need_to_fill)
      countries_analysis_3 <- sapply(countries_analysis, function(x) {cross_check_table %>% filter(level==3) %>% filter(code==x) %>% .$country_iso3c})
      (num_of_countries_analysis <- length(countries_analysis)) # we predict 138 countries

      # prediction #
      predict_data <- expand.grid(origin = 1:num_of_countries_analysis,
                    destination = 1:num_of_countries_analysis,
                    stringsAsFactors = FALSE) %>% as.tibble()
      predict_data$code_o <- countries_analysis[predict_data$origin]
      predict_data$code_d <- countries_analysis[predict_data$destination]
      predict_data$code_3_o <- countries_analysis_3[predict_data$origin]
      predict_data$code_3_d <- countries_analysis_3[predict_data$destination]
      predict_data <- predict_data %>% filter(origin!=destination)
      predict_data <- predict_data %>% filter(!((code_o %in% countries_no_need_to_fill) & (code_d %in% countries_no_need_to_fill)))

      predicted_flow_yearly <- apply(predict_data, 1, function(x){
        # x=unlist(predict_data[9472,])
        code_o <- x[3]
        code_d <- x[4]
        code_o_3 <- x[5]
        code_d_3 <- x[6]
        
        # Get geospatial data
        point1 <- spdf_world[which(spdf_world$adm0_a3 == code_o_3), ]
        point2 <- spdf_world[which(spdf_world$adm0_a3 == code_d_3), ]

        # Extract latitude and longitude
        point1 <- c(point1$label_x[1], point1$label_y[1])
        point2 <- c(point2$label_x[1], point2$label_y[1])

        # Calculate distance (in meters)
        # distance_meters <- distGeo(point1, point2)
        dist <- distVincentySphere(point1, point2)

        # GDP
        gdp_data_o <- gdp_pc_data$NY.GDP.PCAP.CD[gdp_pc_data$iso3c==code_o_3]
        gdp_data_d <- gdp_pc_data$NY.GDP.PCAP.CD[gdp_pc_data$iso3c==code_d_3]

        # Urban population
        pop_data_o <- urban_pop_data$SP.URB.TOTL[urban_pop_data$iso3c==code_o_3]
        pop_data_d <- urban_pop_data$SP.URB.TOTL[urban_pop_data$iso3c==code_d_3]

          # ICAO
        icao_flow_o <- df_ICAO_passenger_i$ICAO_flow[df_ICAO_passenger_i$`Country Code`==code_o_3]
        icao_flow_d <- df_ICAO_passenger_i$ICAO_flow[df_ICAO_passenger_i$`Country Code`==code_d_3]
        if(code_o_3=="USA"){icao_flow_o <- icao_flow_o*usa_adjust_ratio}
        if(code_d_3=="USA"){icao_flow_d <- icao_flow_o*usa_adjust_ratio}

        data_tmp <- tibble(code_o=code_o, code_d=code_d, code_3_o=code_o_3, code_3_d=code_d_3, dist=dist, gdp_o=gdp_data_o, gdp_d=gdp_data_d, upop_o=pop_data_o, upop_d=pop_data_d, icao_flow_o=icao_flow_o, icao_flow_d=icao_flow_d)
        flow <- predict(gravity_model_flow, data_tmp, type = "response")
        # print(data_tmp)
        if(length(flow)==0){return(NA)}else{return(flow)}
        if(is.na(flow)){return(NA)}else{return(flow)}
      })
      predict_data$predicted_flow_yearly <- predicted_flow_yearly
      predict_data <- as.tibble(predict_data)
      predict_data %>% filter(is.na(predicted_flow_yearly))

      # manually adjust the predicted flow for China (both inbound and outbound) to 10% as strict travel ban implemented in China #
      # https://www.mckinsey.com/industries/travel-logistics-and-infrastructure/our-insights/chinas-uneven-travel-recovery-long-road-to-international-travel-furthers-domestic-opportunities
      if(this_year>2019){
        predict_data$predicted_flow_yearly[predict_data$code_o=="CN" | predict_data$code_d=="CN"] <- predict_data$predicted_flow_yearly[predict_data$code_o=="CN" | predict_data$code_d=="CN"]*0.1
      }
      
      # construct the completed daily movement matrix
      # 5.2.3 we supplement the ADS-B data with gravity model estimations and complete the travel matrix #
      full_data <- expand.grid(origin = 1:num_of_countries_analysis,
                    destination = 1:num_of_countries_analysis,
                    date=sort(unique(df_asdb_pairflow_this_year$date)),
                    stringsAsFactors = FALSE) %>% as.tibble()
      full_data$code_o <- countries_analysis[full_data$origin]
      full_data$code_d <- countries_analysis[full_data$destination]
      full_data$code_3_o <- countries_analysis_3[full_data$origin]
      full_data$code_3_d <- countries_analysis_3[full_data$destination]
      full_data <- full_data %>% filter(origin!=destination)
      
      full_data <- left_join(full_data, df_asdb_pairflow_this_year %>% mutate(code_o=ori_iso_country, code_d=dest_iso_country, flow_adsb=N_passenger) %>% select(date, code_o, code_d, flow_adsb), c("date", "code_o", "code_d"))
      
      full_data <- left_join(full_data, predict_data %>% mutate(flow_gravity=predicted_flow_yearly/365.25) %>% select(code_o, code_d, flow_gravity), c("code_o", "code_d"))

      full_data$flow_all <- full_data$flow_adsb
      full_data$data_source <- "adsb"
      full_data$data_source[is.na(full_data$flow_all)] <- "gravity"
      full_data$data_source[is.na(full_data$flow_all) & is.na(full_data$flow_gravity)] <- NA
      full_data$flow_all[is.na(full_data$flow_all)] <- full_data$flow_gravity[is.na(full_data$flow_all)]

      move_mat <- full_data %>% filter(!is.na(flow_all))
      saveRDS(move_mat, here::here(paste0("data/estimated_movement_matrix/move_mat_", this_year, "_", this_level_name, ".rds")))
    })

  }

  if(level==4){ # compete the missing data at level 4
      ## TODO
      ## Not to do now because this study is not focusing on the level 4 data yet

      # 6. Other fixes #
      # Note that Jammu & Kashmir and Ladakh are combined together in the IHME data, only applicable when analysing at level 4 #

  }

  if(level==2){ # compete the missing data at level 4
      ## TODO
      ## Not to do now because this study is not focusing on the level 2 data yet
  }

}
