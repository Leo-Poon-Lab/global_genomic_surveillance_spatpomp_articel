library(tidyverse)
library(ggpubr)
source("scripts/data_processing/helper/color_scheme.R")

# Fig 4
cross_check_table <- readxl::read_excel(here::here("data//our_airports_data/cross_check_table_ihme_input_completed.xlsx")) %>% filter(level==3)

## Fig 4a
### count how many regions prefer community detection, how many regions prefer imported detection
df_fig_4a_label <- read_csv("results/figs/model_simulation/M4/Fig4a_label_source.csv")
df_detection_type_sim <- df_fig_4a_label %>% group_by(code_intro) %>% count(case_type) %>% slice_max(n) %>% select(-n) %>% arrange(case_type, code_intro)

df_s2 <- read_csv("results/figs/model_simulation/M3/Table_S2.csv")
df_s2 <- df_s2 %>% left_join(df_detection_type_sim %>% transmute(code=code_intro, earliest_detection_type=case_type))
# For supplementing Table S2
write_csv(df_s2, "results/figs/model_simulation/M4/Table_S2_supp.csv")

### association between IDR, DSR, travel weight and travel hubs
df_fig_4a_label_resampled <- df_fig_4a_label %>% group_by(code_intro) %>% slice_sample(n = 1024, replace = TRUE)
df_detection_type <- df_fig_4a_label_resampled

df_covar_m4 <- readRDS("results/model_data/data_covariates_M4.rds")

df_detection_type <- left_join(df_detection_type, df_covar_m4 %>% select(code_intro=code, infection_detection, diagnostic_capacity, detection_sequencing, sequencing_capacity) %>% unique(), by = "code_intro")
df_detection_type$travel_weight <- df_detection_type$parval_1
df_detection_type$travel_weight <- gsub("%", "", df_detection_type$travel_weight)
df_detection_type$travel_weight <- as.numeric(df_detection_type$travel_weight)/100

df_s5 <- read_csv("results/figs/model_simulation/M3/Table_S4.csv")
df_s5$num_travel_hubs <- as.numeric(gsub(":.+", "", df_s5$info_travelhubs))
df_detection_type <- left_join(df_detection_type, df_s5 %>% transmute(parval_2=st_abbr, travel_hubs=num_travel_hubs) %>% unique(), by="parval_2")

df_plot <- df_detection_type %>% ungroup() %>% 
  select(code_intro, case_type, infection_detection, diagnostic_capacity, detection_sequencing, sequencing_capacity, travel_weight) %>%
  pivot_longer(cols=c(infection_detection, diagnostic_capacity, detection_sequencing, sequencing_capacity, travel_weight), names_to="parameter", values_to="value")
df_plot$parameter <- factor(df_plot$parameter, levels=c("infection_detection", "diagnostic_capacity", "detection_sequencing", "sequencing_capacity", "travel_weight"), labels = c("IDR", "Diagnostic capacity", "DSR", "Sequencing capacity", "Traveler weight"))
df_plot$case_type <- factor(df_plot$case_type, levels=c("Community", "Imported"), labels=c("Community-detected", "Traveler-detected"))

p_fig_s5 <- ggplot(df_plot, aes(x=case_type, y=value, fill=case_type, group=case_type))+
  geom_boxplot(position="dodge", alpha=0.8)+
  # geom_point(aes(x=case_type, y=value, color=code_intro, group=case_type), position=position_jitterdodge(), alpha=0.5, show.legend = FALSE)+
  geom_point(aes(x=case_type, y=value, color=code_intro, group=case_type), position=position_jitterdodge(), alpha=0.05, size=0.5)+
  stat_compare_means(size=3, vjust=-2)+
  scale_color_manual(name="Regions", values=colors_spatial_units)+
  scale_fill_manual(name="Earliest detection type", values=c("darkblue", "darkred"))+
  facet_wrap(~parameter, nrow = 1, scales="free")+
  scale_y_continuous(expand = c(0.1,0))+
  theme_minimal()+
  xlab("")+
  ylab("")+
  theme(
    axis.text.x = element_text(angle = 15),
    legend.box = "vertical",
    legend.box.spacing = unit(0.1, "lines"),
    legend.key.height = unit(0.4, "lines"),
    legend.key.width = unit(1, "lines"),
    legend.text = element_text(size=8),
    legend.position = "bottom")+
  guides(color = guide_legend(nrow = 4, byrow = TRUE, override.aes = list(alpha = 1, size = 3)))+
  NULL
ggsave("results/figs/model_simulation/M4/Fig_S6.jpg", p_fig_s5, width=10, height=10, dpi=300)
ggsave("results/figs/model_simulation/M4/Fig_S6.pdf", p_fig_s5, width=10, height=10)

## Fig 4c
### how many days are saved by using the recommended strategy, compared to baseline
list_p_fig_4abcd <- readRDS("results/figs/model_simulation/M4/list_p_fig_4abcd.rds")
(best_param_set <- list_p_fig_4abcd$df_EDT_fig_4bcd_percentile %>%
  filter((parval_2=="P2-H" & parval_1=="0.1%") | (parval_2=="T2-R" & parval_1=="70%")))

lags_P_selected <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="0.1%", parval_2=="P2-H") %>% pull(lag)
lags_T_selected <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="70%", parval_2=="T2-R") %>% pull(lag)
lags_baseline <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="1%", parval_2=="Max") %>% pull(lag)

t.test(lags_P_selected, lags_baseline)
t.test(lags_T_selected, lags_baseline)
