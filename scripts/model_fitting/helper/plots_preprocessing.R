plot_gisaid <- function(data_fitting, outprefix="GISAID_raw_"){
  p_GISAID <- data_fitting$gisaid$data_GISAID %>% ggplot() +
    geom_line(aes(x=`Collection date`, y=N+1, group=lineage_new, color=lineage_new)) +
    facet_wrap(~code, ncol=3)+
    scale_color_manual(values=colors_lineage)+
    scale_y_log10(  # Plot on log-scale
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  if (!exists("support_ggsave")) {
    support_ggsave <- TRUE
  }
  if (support_ggsave) {
    ggsave(here::here(paste0("results/figs/", outprefix, model_name, ".jpg")), plot=p_GISAID, width=12*0.618, height=12, dpi = 450)
  }
}

plot_in_out_flow <- function(data_fitting){
  ## Data inspection #
  ### confirm that for each country, the inflow and outflow are relatively stable over time. ###
  data_plot_outflow <- data_fitting$travel$mov_mat %>% mutate(code=code_o) %>% group_by(date, code) %>% summarise(outflow=sum(flow_all)) %>% ungroup()
  data_plot_inflow <- data_fitting$travel$mov_mat %>% mutate(code=code_d) %>% group_by(date, code) %>% summarise(inflow=sum(flow_all)) %>% ungroup()

  data_plot_flow <- left_join(data_plot_outflow, data_plot_inflow, by=c("date", "code")) %>% mutate(ratio=outflow/inflow)
  data_plot_flow <- left_join(data_plot_flow, cross_check_table %>% select(code, loc_name))
  data_plot_flow$loc_name[is.na(data_plot_flow$loc_name)] <- data_plot_flow$code[is.na(data_plot_flow$loc_name)]

  df_meas <- data_fitting$df_meas
  df_meas <- left_join(df_meas, cross_check_table %>% select(code, loc_name))
  df_meas$loc_name[is.na(df_meas$loc_name)] <- df_meas$code[is.na(df_meas$loc_name)]

  data_plot_flow <- data_plot_flow %>% pivot_longer(cols=c("inflow", "outflow"), names_to="flow_type", values_to="flow")
  df_meas$flow_type <- "cases"
  df_meas$flow <- df_meas$daily_cases

  p_in_out_flow_log <- bind_rows(data_plot_flow, df_meas) %>% group_by(date, loc_name, flow_type) %>% summarise(flow=mean(flow)) %>% ggplot()+
    geom_line(aes(x=date, y=log10(flow+1), group=flow_type, color=flow_type))+
    facet_wrap(~loc_name, ncol=5)+
    scale_color_manual(values=c("black", "red", "blue"))+
    NULL
  if (!exists("support_ggsave")) {
    support_ggsave <- TRUE
  }
  if (support_ggsave) {
    ggsave(here::here(paste0("results/figs/in_out_flow_log_", model_name, ".jpg")), plot=p_in_out_flow_log, width=12, height=15, dpi = 450)
  }

  p_in_out_flow <- bind_rows(data_plot_flow, df_meas) %>% group_by(date, loc_name, flow_type) %>% summarise(flow=mean(flow)) %>% ggplot()+
    geom_line(aes(x=date, y=flow, group=flow_type, color=flow_type))+
    facet_wrap(~loc_name, ncol=5)+
    scale_color_manual(values=c("black", "red", "blue"))+
    NULL
  if (!exists("support_ggsave")) {
    support_ggsave <- TRUE
  }
  if (support_ggsave) {
    ggsave(here::here(paste0("results/figs/in_out_flow_", model_name, ".jpg")), plot=p_in_out_flow, width=12, height=15, dpi = 450)
  }
}
