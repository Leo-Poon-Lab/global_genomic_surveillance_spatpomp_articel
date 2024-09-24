#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("input/output dir must be supplied.n", call.=FALSE)
} else {
  input_dir <- args[1]
  output_dir <- args[2]
}
# input_dir <- "results/model_data/profiling_HKU/profiling_unit_1/"
# output_dir <- "results/figs/profiling_HKU/profiling_unit_1/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files_rds <- list.files(input_dir, full.names = TRUE, pattern = "profile_.+\\.rds$")

# ggsave may not run on servers lackzing X11 support
suppressPackageStartupMessages(library(tidyverse))
source("scripts/model_fitting/helper/parameter_names_est.R")

if(all(grepl("_low.rds$", files_rds))){
  search_level <- "low"
} else if(all(grepl("_mid.rds$", files_rds))){
  search_level <- "mid"
} else if(all(grepl("_high.rds$", files_rds))){
  search_level <- "high"
} else{
  stop("The files in the input directory should be all low, mid or high.")
}

for (file_rds_i in files_rds) {
  # file_rds_i <- files_rds[1]
  profile_para <- gsub(paste0("_", search_level, ".rds"), "", gsub("profile_", "", basename(file_rds_i)), fixed=T)
  
  profile_data <- readRDS(file_rds_i)
  profile_data <- cbind(profile_data$logLiks,profile_data$params)
  if(any(names(profile_data)=="ll")){
    names(profile_data)[names(profile_data)=="ll"] <- "logLik"
  }
  if(profile_para %in% shared_param_names_est){
    profile_para1 <- paste0(profile_para,"1")
  }
  if(profile_para %in% paste0(rep(unit_specific_names_est, each = 29), seq_len(29))){
    profile_para1 <- profile_para
  }
  prof_res <- profile_data %>%
    group_by(.data[[profile_para1]]) %>%
    summarize(logLik = max(logLik))

  mcap_results <- pomp::mcap(prof_res$logLik, prof_res[[profile_para1]])

  profile_plot <- ggplot() +
    geom_point(data = prof_res, aes(x = .data[[profile_para1]], y = logLik)) +
    geom_line(data = mcap_results$fit, aes(x = parameter, y = smoothed), col = 'blue') +
    geom_vline(xintercept = mcap_results$ci[1], linetype = 'dashed') +
    geom_vline(xintercept = mcap_results$ci[2], linetype = 'dashed') +
    geom_vline(xintercept = mcap_results$mle, col = 'blue') +
    labs(x = profile_para, y = 'Log Likelihood') +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  # This ggsave may not run on servers lacking X11 support
  ggsave(paste0(output_dir, "/", profile_para,".png"), plot=profile_plot, width=10, height=7, dpi=300)
}
