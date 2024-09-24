library(tidyverse)
library(geofacet)
library(ggrepel)
source("scripts/data_processing/helper/color_scheme.R")

# Fig 2
dir_rst_fig2 <- "results/figs/model_simulation/Omicron20/"

# Fig 2a
## The peaks
df_fig_2a <- readRDS("results/figs/model_simulation/Omicron20/df_fig_2a_source.rds")
(df_fig_2a_peak <- df_fig_2a %>% group_by(cases_level, date) %>% summarise(total_value=sum(value)) %>% group_by(cases_level) %>% slice_max(total_value))

# Fig 2b
## Why did not HK diagnostic BA.1 imported at very early stage?
## infected cases are being detected at around 2021-11-10
data_covar <- readRDS("results/model_data/data_covariates_Omicron20.rds")
data_covar %>% filter(code=="HK") %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date, infection_detection, detection_sequencing) %>% filter(date>ymd("2021-11-01"), date<ymd("2021-12-01")) %>% print(n=200)
### This is because the infection-diagnostic ratio is very low in HK. Despite the diagnostic-sequencing ratio is high (~50% to 100%), the infection-diagnostic ratio is 1% and the sequencing rate is ~1%.

## plot the IDR and DSR for different countries
source("scripts/model_fitting/helper/geo_grid.R")
layout_geo_grid_29 <- geo_grid_29()
layout_geo_grid_29$code_intro <- layout_geo_grid_29$code
layout_geo_grid_29$name_intro <- layout_geo_grid_29$name

df_plot_IDR_DSR <- data_covar %>% mutate(date=as_date(date_decimal(date_decimal))) %>% select(date_decimal, date, code, infection_detection, detection_sequencing, travel_control_level) %>% mutate(infection_detection=infection_detection*100, detection_sequencing=detection_sequencing*100, travel_control_level=travel_control_level*25) %>% pivot_longer(cols = c("infection_detection", "detection_sequencing", "travel_control_level"), names_to = "Type", values_to = "Value") %>% mutate(Type=factor(Type, levels=c("infection_detection", "detection_sequencing", "travel_control_level"), labels=c("IDR", "DSR", "Travel control level")))

df_plot_IDR_DSR %>% filter(Type=="IDR") %>% summarise(mean_value=mean(Value), sd_value=sd(Value))
df_plot_IDR_DSR %>% filter(Type=="DSR") %>% summarise(mean_value=mean(Value), sd_value=sd(Value))

p_IDR_DSR <-ggplot(data = df_plot_IDR_DSR) +
  geom_line(aes(x = date, y = Value, color = Type), linewidth = 0.7, linetype="solid", alpha=0.9) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(name = "Covariate",values = colors_IDR_DSR_TraLevel) +
  xlab("Date") +
  facet_geo(vars(code), grid = layout_geo_grid_29, label="name") +
  theme_minimal()+
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(angle=30),
    legend.position = "bottom"
    )+
  guides(size="none")+
  NULL
ggsave(p_IDR_DSR, filename=paste0(dir_rst_fig2, "fig_IDR_DSR_TraLevel.jpg"), width=12, height=10, dpi=300)
ggsave(p_IDR_DSR, filename=paste0(dir_rst_fig2, "fig_IDR_DSR_TraLevel.pdf"), width=12, height=10)

## compare the medians in fig 2b, between community and imported cases
df_fig_2b <- readRDS("results/figs/model_simulation/Omicron20/df_fig_2b_source.rds")
### Infected
df_fig_2b %>% group_by(unitname, cases_level, variants, lineage_type) %>% summarise(median_EDT=median(earliest_date)) %>% ungroup() %>% filter(cases_level=="Infected") %>% group_by(unitname, variants) %>% summarise(is_imported_earlier=median_EDT[which(lineage_type=="Imported")]<median_EDT[which(lineage_type=="Community")]) %>% arrange(variants, is_imported_earlier) %>% print(n=200)
### Diagnostic
df_fig_2b %>% group_by(unitname, cases_level, variants, lineage_type) %>% summarise(median_EDT=median(earliest_date)) %>% ungroup() %>% filter(cases_level=="Diagnostic") %>% group_by(unitname, variants) %>% summarise(is_imported_earlier=median_EDT[which(lineage_type=="Imported")]<median_EDT[which(lineage_type=="Community")]) %>% arrange(variants, is_imported_earlier) %>% print(n=200)
### Sequenced
df_fig_2b %>% group_by(unitname, cases_level, variants, lineage_type) %>% summarise(median_EDT=median(earliest_date)) %>% ungroup() %>% filter(cases_level=="Sequenced") %>% group_by(unitname, variants) %>% summarise(is_imported_earlier=median_EDT[which(lineage_type=="Imported")]<median_EDT[which(lineage_type=="Community")]) %>% arrange(variants, is_imported_earlier) %>% print(n=200)

## How long does it take for a variant to arrive? (emergence-arrival lag)
lags_intro <- readxl::read_excel("results/figs/model_simulation/Omicron20/params_transformed_Omicron20.xlsx") %>% select(day_BA_one, day_BA_two) %>% unique() %>% unlist()
df_vline <- tibble(lag=lags_intro)
df_vline$date_decimal_intro <- floor(df_vline$lag)/365.25 + min(model_sims$date_decimal)
df_vline$date_intro <- date_decimal(df_vline$date_decimal_intro)
df_vline$variants <- c("Omicron BA.1", "Omicron BA.2")

df_fig_2b_median <- df_fig_2b %>% group_by(unitname, cases_level, variants, lineage_type) %>% summarise(median_EDT=median(earliest_date)) 
df_fig_2b_median <- left_join(df_fig_2b_median, df_vline %>% select(variants, date_intro), by="variants")

df_fig_2b_median_earliest <- df_fig_2b_median %>% group_by(unitname, cases_level, variants) %>% slice_min(median_EDT)

## The emergence-arrival lag, arrival-diagnostic and diagnostic-sequencing lags
df_fig_2b_median_earliest %>% select(-lag) %>% pivot_wider(names_from = cases_level, values_from = median_EDT) %>% mutate(infection_arrival_lag=as.numeric(difftime(Infected, date_intro, units="days")), arrival_diagnostic_lag=as.numeric(difftime(Diagnostic, Infected, units="days")), diagnostic_sequencing_lag=as.numeric(difftime(Sequenced, Diagnostic, units="days"))) %>% group_by(variants) %>% summarise(mean_infection_arrival_lag=mean(infection_arrival_lag, na.rm=T), mean_arrival_diagnostic_lag=mean(arrival_diagnostic_lag, na.rm=T), mean_diagnostic_sequencing_lag=mean(diagnostic_sequencing_lag, na.rm=T))

## Does travel hubs have shorter arrival-diagnostic and diagnostic-sequencing lags?
df_travel_hubs <- read_csv("results/figs/model_simulation/Omicron20/travel_hub_rank.csv")
names(df_travel_hubs)[1] <- "unitname"
df_plot <- df_fig_2b_median_earliest %>% select(-lag) %>% pivot_wider(names_from = cases_level, values_from = median_EDT) %>% mutate(infection_arrival_lag=as.numeric(difftime(Infected, date_intro, units="days")), arrival_diagnostic_lag=as.numeric(difftime(Diagnostic, Infected, units="days")), diagnostic_sequencing_lag=as.numeric(difftime(Sequenced, Diagnostic, units="days"))) %>% select(unitname, variants, lineage_type, ends_with("lag")) %>% pivot_longer(cols = c("infection_arrival_lag", "arrival_diagnostic_lag", "diagnostic_sequencing_lag"), names_to = "Type", values_to = "Value") %>% filter(!is.na(Value))
df_plot <- left_join(df_plot, df_travel_hubs, by="unitname")
df_plot$group <- paste0(df_plot$variants, " (", df_plot$lineage_type, ")")
df_plot$Type <- factor(df_plot$Type, levels=c("infection_arrival_lag", "arrival_diagnostic_lag", "diagnostic_sequencing_lag"), labels=c("Infection-arrival lag", "Arrival-diagnostic lag", "Diagnostic-sequencing lag"))

p_lag_travel <- ggplot(df_plot, aes(x=travel_per_capita, y=Value, color=variants))+
  geom_point(alpha=0.7)+
  geom_text_repel(aes(x=travel_per_capita, y=Value, label=unitname), size=2, nudge_x=0.1, nudge_y=0.1)+
  scale_color_manual(values=colors_lineage)+
  facet_grid(Type~group)+
  ylab("Lag (days)")+
  xlab("Traveler per capita")+
  theme_minimal()+
  theme(
    legend.position = "bottom"
    )+
  geom_smooth(method='lm')+
  ylim(0, 80)+
  NULL
ggsave(p_lag_travel, filename=paste0(dir_rst_fig2, "fig_lag_travel.jpg"), width=10, height=10, dpi=300)
ggsave(p_lag_travel, filename=paste0(dir_rst_fig2, "fig_lag_travel.pdf"), width=10, height=10)

# Fig 2c
df_plot_2c <- readRDS("results/figs/model_simulation/Omicron20/df_plot_2c_source.rds")
## The peak of the daily travel-related infections for BA.1 and BA.2
(df_plot_2c_peak <- df_plot_2c %>% group_by(variants, date) %>% summarise(total_value=sum(value)) %>% group_by(variants) %>% slice_max(total_value))

df_plot_2c_peak$date[df_plot_2c_peak$variants=="Omicron BA.1"] - df_fig_2a_peak$date[df_fig_2a_peak$cases_level=="Infected"]
df_plot_2c_peak %>% left_join(df_vline %>% transmute(variants=gsub(" emergence", "", v_line_label), date_intro=date), by="variants") %>% mutate(lag=as.numeric(difftime(date, date_intro, units="days")))

df_plot_2c %>% filter(date==df_plot_2c_peak$date[df_plot_2c_peak$variants=="Omicron BA.2"]) %>% group_by(variants) %>% summarise(total_value=sum(value)) # On 2022-02-26

## The time for the seeding period (during seeding period the total daily travel-related infection should not exceed 100)
(df_seeding_period <- df_plot_2c %>% group_by(variants, date) %>% summarise(total_value=sum(value)) %>% filter(total_value<500) %>% group_by(variants) %>% summarise(date=max(date)) %>% left_join(df_vline %>% transmute(variants=gsub(" emergence", "", v_line_label), date_intro=date), by="variants") %>% mutate(seeding_period=as.numeric(difftime(date, date_intro, units="days"))))

## The proportion for cases of diffent origins, among all the travel-related infections, in the seeding period and out of the seeding period.
df_plot_2c %>% left_join(df_seeding_period %>% transmute(date_seeding=date, variants)) %>% mutate(period=ifelse(date<date_seeding, "Seeding", "After seeding")) %>% group_by(variants, period, origin_id) %>% summarise(origin_period_value=sum(value)) %>% ungroup() %>% group_by(variants, period) %>% mutate(total_period_value=sum(origin_period_value)) %>% mutate(proportion=origin_period_value/total_period_value) %>% group_by(variants, period) %>% arrange(variants, desc(proportion)) %>% slice_head(n=3) %>% print(n=50)

df_plot_2c %>% left_join(df_seeding_period %>% transmute(date_seeding=date, variants)) %>% mutate(period=ifelse(date<date_seeding, "Seeding", "After seeding")) %>% group_by(variants, period, origin_id) %>% summarise(origin_period_value=sum(value)) %>% ungroup() %>% group_by(variants, period) %>% mutate(total_period_value=sum(origin_period_value)) %>% mutate(proportion=origin_period_value/total_period_value) %>% filter(origin_id=="ZA")
