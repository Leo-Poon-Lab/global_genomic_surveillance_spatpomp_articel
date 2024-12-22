library(tidyverse)
library(patchwork)

# Now that we have generated some simulations from `Simulate_M4_4.R`.

## Run on local machine
model_name <- "M4"
base_model_name <- "Omicron20"
n_cores <- 64
dir_rst <- paste0("results/figs/model_simulation/", model_name, "/")
dir.create(dir_rst, showWarnings = FALSE, recursive = TRUE)
dir_model_data <- paste0("results/model_data/model_simulation/", model_name, "/")

# ## Load data
data_fitting_rds_path <- paste0("results/model_data/data_fitting_", base_model_name, ".rds")
if(!file.exists(data_fitting_rds_path)){
  stop("Please run the data processing script in the M_0 model first.")
} else{
  data_fitting <- readRDS(data_fitting_rds_path)
}

df_all_values <- list.files(dir_model_data, "df_all_values_M4_2\\d{3}.rds$", full.names = TRUE) %>% naturalsort::naturalsort() %>% lapply(readRDS) %>% bind_rows()

## combine all fig 4 sub plots
source("scripts/model_simulation/helper/plot_fig_5d.R")
## Plotting Fig. 5
list_fig_5d <- plot_fig_5d(
  dir_rst = dir_rst,
  dir_model_data = dir_model_data,
  df_all_values = df_all_values,
  data_fitting = data_fitting,
  bootstrap_each_param=1024*10,
  n_cores = n_cores
)

list_fig_5abc <- readRDS(paste0(dir_rst, "list_p_fig_5abc.rds"))
list_fig_5d <- readRDS(paste0(dir_rst, "list_p_fig_5d.rds"))
p_fig_5bcd <- (list_fig_5abc$p_fig_5bc / (list_fig_5d$p_fig_5d + ggtitle("d"))) + plot_layout(heights = c(1, 2))
ggsave(p_fig_5bcd, filename = paste0(dir_rst, "fig_5bcd.jpg"), width = 10, height = 10, units = "in", dpi = 300)
ggsave(p_fig_5bcd, filename = paste0(dir_rst, "fig_5bcd.pdf"), width = 10, height = 10, units = "in")

p_fig_5 <- (list_fig_5abc$p_fig_5a + ggtitle("a")) + p_fig_5bcd + 
  plot_layout(nrow = 1, widths = c(3,2))
ggsave(p_fig_5, filename = paste0(dir_rst, "fig_5.jpg"), width = 20, height = 20*0.618, units = "in", dpi = 300)
ggsave(p_fig_5, filename = paste0(dir_rst, "fig_5.pdf"), width = 20, height = 20*0.618, units = "in")
