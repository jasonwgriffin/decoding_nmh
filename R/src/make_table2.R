make_table2 <- function(fpath, d, iq_control){
  
  d = d_t1
  fpath=fig_path
  
  df_list <- list(
    "All Images" =           tidyr::pivot_wider(d$all_images_metrics,      names_from = metric,values_from = value) %>% mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% left_join(phenotype %>% select(id,visit_participant_group,diagnosis_summary_full_scale_best_iq.t1), by = "id") %>% rename(iq=diagnosis_summary_full_scale_best_iq.t1),
    "Faces vs. Houses" =     tidyr::pivot_wider(d$faces_houses_metrics,    names_from = metric,values_from = value) %>% mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% left_join(phenotype %>% select(id,visit_participant_group,diagnosis_summary_full_scale_best_iq.t1), by = "id") %>% rename(iq=diagnosis_summary_full_scale_best_iq.t1),
    "Face Identity" =        tidyr::pivot_wider(d$face_identity_metrics,   names_from = metric,values_from = value) %>% mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% left_join(phenotype %>% select(id,visit_participant_group,diagnosis_summary_full_scale_best_iq.t1), by = "id") %>% rename(iq=diagnosis_summary_full_scale_best_iq.t1),
    "Upright vs. Inverted" = tidyr::pivot_wider(d$upright_inverted_metrics,names_from = metric,values_from = value) %>% mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% left_join(phenotype %>% select(id,visit_participant_group,diagnosis_summary_full_scale_best_iq.t1), by = "id") %>% rename(iq=diagnosis_summary_full_scale_best_iq.t1)
    )

  # List of metrics
  metrics <- c("decoding_onset", "decoding_sus", "decoding_peak", "decoding_latency")
  # Mapping metrics to simpler names for table columns
  metric_names <- c(decoding_onset = "onset", decoding_sus = "sustainability", decoding_peak = "peak",decoding_latency = "latency")
  
  # Initialize empty list to store results
  tab_list <- list()
  
  # Loop over datasets
  for(type_name in names(df_list)){
    df <- df_list[[type_name]]
    
    # Loop over metrics
    for(metric in metrics){
      
      # Run ANOVA
      if (iq_control == "yes") {
        a <- anova(lm(as.formula(paste0(metric, "~ group * age_cat + iq")), data = df))
      } else {
        a <- anova(lm(as.formula(paste0(metric, "~ group * age_cat")), data = df))
      }
      
      a_tidy <- broom::tidy(a)
      
      # Residual df for formatting
      df_res <- ifelse(iq_control =="yes", a_tidy$df[5],a_tidy$df[4])
      
      # Clean up table
      a_tidy <- a_tidy %>%
        filter(term != "Residuals") %>%
        mutate(
          type = type_name,
          df = paste0(df, ",", df_res),
          f = statistic,
          p = p.value,
          metric = metric_names[[metric]]
        ) %>%
        select(type, term, df, f, p, metric)
      
      # Append to list
      tab_list[[paste(type_name, metric, sep = "_")]] <- a_tidy
    }
  }
  
  # Combine all results
  tab2 <- bind_rows(tab_list)
  
  # Pivot to wide format
  tab2 <- tab2 %>%
    pivot_wider(names_from = metric, values_from = c(df, f, p)) %>%
    select(1,2,3,
           f_onset, p_onset,
           f_sustainability, p_sustainability,
           f_peak, p_peak,
           f_latency, p_latency) %>%
    mutate_if(is.numeric, round, 3)
  
  if (iq_control == "yes") {
    write.csv(tab2, file = paste0(fpath, "table2_iq_control.csv"), row.names = F)
  } else {
    write.csv(tab2, file = paste0(fpath, "table2.csv"), row.names = F)
  }
  
  
   
}