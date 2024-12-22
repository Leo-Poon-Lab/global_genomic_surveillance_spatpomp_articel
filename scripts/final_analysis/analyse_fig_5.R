library(tidyverse)

## recap the settings for Figure 5
dir_rst <- paste0("results/figs/model_simulation/", "M4", "/")
list_p_fig_4abcd <- readRDS(paste0(dir_rst, "list_p_fig_4abcd.rds"))
list_fig_5abc <- readRDS(paste0(dir_rst, "list_p_fig_5abc.rds"))

(best_param_set <- list_p_fig_4abcd$df_EDT_fig_4bcd_percentile %>%
  filter((parval_2=="P2-H" & parval_1=="0.1%") | (parval_2=="T2-R" & parval_1=="70%")))

list_fig_5abc$change_IDR_min_better
list_fig_5abc$change_DSR_min_better

## fig 5b,c
list_p_fig_5abc <- readRDS("results/figs/model_simulation/M4/list_p_fig_5abc.rds")
df_baseline <- readRDS("results/figs/model_simulation/M4/df_EDT_fig_4bcd.rds") %>% filter(parval_1=="1%", parval_2=="Max")

names(list_p_fig_5abc)

df_fig_5b <- list_p_fig_5abc$p_fig_5bc[[1]]$data
lags_P2H_010 <- df_fig_5b %>% filter(parval=="10%", strategy=="P2-H") %>% pull(lag) 
t.test(lags_P2H_010, df_baseline$lag)
round(mean(df_baseline$lag)-mean(lags_P2H_010),2)

lags_T2R_001 <- df_fig_5b %>% filter(parval=="1%", strategy=="T2-R") %>% pull(lag) 
t.test(lags_T2R_001, df_baseline$lag)
round(mean(df_baseline$lag)-mean(lags_T2R_001),2)

df_fig_5c <- list_p_fig_5abc$p_fig_5bc[[2]]$data
lags_P2H_040 <- df_fig_5c %>% filter(parval=="40%", strategy=="P2-H") %>% pull(lag)
t.test(lags_P2H_040, df_baseline$lag)
round(mean(df_baseline$lag)-mean(lags_P2H_040),2)

lags_T2R_060 <- df_fig_5c %>% filter(parval=="60%", strategy=="T2-R") %>% pull(lag)
t.test(lags_T2R_060, df_baseline$lag)
round(mean(df_baseline$lag)-mean(lags_T2R_060),2)

# In discussion, we want to calculate the actual number needed for diagnostics and sequencing in travels

## T2-R with 70% traveler weight under the half resources scenario (Figure 5d) would require the United States ...
### diagnosing N travelers
745944 * # the diagnosting capacity
0.7 * # the traveler weight
0.5 # the half resources scenario

### sequencing N travelers
14675 * # the sequencing capacity
0.7 * # the traveler weight
0.5 # the half resources scenario

## P2-H with 0.1% traveler weight under the half resources scenario (Figure 5d) would require the Hong Kong ...
### diagnosing N travelers
#### HK
6675/100 * # the diagnosing capacity
0.001 * # the traveler weight
0.5 * # the half resources scenario
39046 # the number of diagnostic capacity
#### UAE
37020/100 * # the diagnosing capacity
0.001 * # the traveler weight
0.5 * # the half resources scenario
2801 # the number of diagnostic capacity

### sequencing N times
#### HK
36049/100 * # the sequencing capacity
0.001 * # the traveler weight
0.5 * # the half resources scenario
77 # the number of sequencing capacity
#### UAE
60706/100 * # the sequencing capacity
0.001 * # the traveler weight
0.5 * # the half resources scenario
20 # the number of sequencing capacity