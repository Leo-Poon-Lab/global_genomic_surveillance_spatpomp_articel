library(tidyverse)
library(ggpubr)
source("scripts/data_processing/helper/color_scheme.R")

# Fig 4
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

## Fig 4a
### count how many regions prefer community detection, how many regions prefer imported detection
df_fig_4a_label <- read_csv("results/figs/model_simulation/M4/Fig4a_label_source.csv")
df_detection_type_sim <- df_fig_4a_label %>% group_by(code_intro) %>% count(case_type) %>% slice_max(n) %>% select(-n) %>% arrange(case_type, code_intro) %>% print(n=30)
(df_detection_type_sim <- df_detection_type_sim %>% left_join(cross_check_table %>% transmute(code_intro=code, continent, loc_name)) %>% arrange(case_type, continent, code_intro) %>% print(n=30)) # For supplementing Table S2

# reorder df_detection_type_sim by Table S2
code_travel_order <- c("AE", "HK", "ES", "GB", "DE", "FR", "CA", "EU_others", "NZ", "US", "TR", "AU", "KR", "NA_others", "MA", "JP", "OC_others", "CL", "MX", "AR", "ZA", "SA_others", "BR", "AS_others", "NG", "PG", "AF_others", "IN", "TN")
df_detection_type_sim[match(code_travel_order, df_detection_type_sim$code_intro),] %>% write_csv("results/figs/model_simulation/M4/Table_S2_supplement.csv")

df_detection_type <- df_fig_4a_label %>% group_by(code_intro,parval_1,parval_2) %>% count(case_type) %>% arrange(case_type, code_intro) %>% select(-n)

### association between IDR, DSR, travel weight and travel hubs
df_covar_m4 <- readRDS("results/model_data/data_covariates_M4.rds")

df_detection_type <- left_join(df_detection_type, df_covar_m4 %>% select(code_intro=code, infection_detection, detection_sequencing) %>% unique())
df_detection_type$travel_weight <- df_detection_type$parval_1
df_detection_type$travel_hubs <- as.numeric(df_detection_type$parval_2)
df_detection_type$travel_hubs[is.na(df_detection_type$travel_hubs)] <- 29

df_plot <- df_detection_type %>% ungroup() %>% select(-parval_1, -parval_2) %>% pivot_longer(cols=c(infection_detection, detection_sequencing, travel_weight, travel_hubs), names_to="parameter", values_to="value")
df_plot$parameter <- factor(df_plot$parameter, levels=c("infection_detection", "detection_sequencing", "travel_weight", "travel_hubs"), labels = c("IDR", "DSR", "Traveler weight", "Travel hubs"))
df_plot$case_type <- factor(df_plot$case_type, levels=c("Community", "Imported"), labels=c("Community-detected", "Traveler-detected"))

p_fig_s5 <- ggplot(df_plot, aes(x=case_type, y=value, fill=case_type, group=case_type))+
  geom_boxplot(position="dodge", alpha=0.8)+
  geom_point(aes(x=case_type, y=value, color=code_intro, group=case_type), position=position_jitterdodge(), alpha=0.5, show.legend = FALSE)+
  stat_compare_means(size=3, vjust=-2)+
  scale_color_manual(values=colors_spatial_units)+
  scale_fill_manual(values=c("darkblue", "darkred"))+
  facet_wrap(~parameter, nrow = 1, scales="free")+
  scale_y_continuous(expand = c(0.1,0))+
  theme_minimal()+
  xlab("")+
  ylab("Value")+
  theme(
    axis.text.x = element_text(angle = 15),
    legend.position = "bottom")+
  NULL
ggsave("results/figs/model_simulation/M4/Fig_S5.jpg", p_fig_s5, width=8, height=6, dpi=300)
ggsave("results/figs/model_simulation/M4/Fig_S5.pdf", p_fig_s5, width=8, height=6)

## Elaborating the recommended strategy, how many fold increase in what regions? What exactly IDR and DSR are required for these regions? What continents are covered?
dir_model_data <- paste0("results/model_data/model_simulation/", "M4", "/")
df_all_values <- list.files(dir_model_data, "df_all_values_M4_\\d{1,2}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

df_s2_IDR_DSR <- tibble(code=sort(code_travel_order), change_IDR=df_all_values %>% filter(scenario=="TH_SP_ID_9_0.01_0") %>% pull(change_IDR) %>% unlist(), change_DSR=df_all_values %>% filter(scenario=="TH_SP_ID_9_0.01_0") %>% pull(change_DSR) %>% unlist())

data_covar_m4 <- readRDS("results/model_data/data_covariates_M4.rds")
df_s2_IDR_DSR <- data_covar_m4 %>% transmute(code, IDR_before=infection_detection, DSR_brfore=detection_sequencing) %>% unique() %>% left_join(df_s2_IDR_DSR)

df_s2_IDR_DSR <- df_s2_IDR_DSR %>% mutate(IDR_after=IDR_before*change_IDR, DSR_after=DSR_brfore*change_DSR) 
df_s2_IDR_DSR <- df_s2_IDR_DSR[match(code_travel_order, df_s2_IDR_DSR$code),]
write_csv(df_s2_IDR_DSR, "results/figs/model_simulation/M4/Table_S2_IDR_DSR.csv")

## Fig 4c
### how many days are saved by using the recommended strategy, compared to baseline
df_baseline <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="0.01", parval_2=="Max")

df_fig_4c_label <- readRDS("results/figs/model_simulation/M4/Fig4c_label_source.rds")

lags_recommended_strategy <- df_fig_4c_label %>% filter(parval_1=="1", parval_2=="9") %>% pull(lag)
lags_baseline <- df_baseline$lag
t.test(lags_recommended_strategy, lags_baseline)


## fig 4f,g
list_p_fig_4efg <- readRDS("results/figs/model_simulation/M4/list_p_fig_4efg.rds")

df_fig_4f <- list_p_fig_4efg$p_fig_4fg[[1]]$data
lags_IDR_001 <- df_fig_4f %>% filter(parval==0.01) %>% pull(lag) 
t.test(lags_IDR_001, df_baseline$lag)
mean(df_baseline$lag)-mean(lags_IDR_001)

df_fig_4g <- list_p_fig_4efg$p_fig_4fg[[2]]$data
lags_DSR_0.6 <- df_fig_4g %>% filter(parval==0.6) %>% pull(lag)
t.test(lags_DSR_0.6, df_baseline$lag)

lags_DSR_1 <- df_fig_4g %>% filter(parval==1) %>% pull(lag)
t.test(lags_DSR_1, df_baseline$lag)
mean(df_baseline$lag)-mean(lags_DSR_1)

set.seed(2024)
(sample(lags_DSR_1, 10000, replace=TRUE)-sample(df_baseline$lag, 10000, replace=TRUE)) %>% quantile(seq(0,1,0.01))
