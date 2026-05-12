## this is for the CNR measure
tv_files <- list.files("") # director of voltage files

tv_files %>%
  set_names(basename(.) %>% str_remove("_trial_voltages\\.csv$")) %>% 
  imap(~ {
    
    out <- read.csv(.x) %>%
      filter(Electrode < 121) %>%
      group_by(Trial, Electrode, Bin) %>%
      summarize(volt = mean(Voltage), .groups = "drop")
    
    write_csv(out, paste0("../../decoding_manuscript/Submission/Revision/trial_data_averages/",.y, "_avg.csv"))
    
  }, .progress = TRUE)

compute_cnr <- function(df) {
  
  # Run ANOVA
  model <- aov(volt ~ Bin * Electrode, data = df)
  anova_table <- anova(model)
  
  # Extract SS
  SS_side        <- anova_table["Bin", "Sum Sq"]
  SS_electrode   <- anova_table["Electrode", "Sum Sq"]
  SS_interaction <- anova_table["Bin:Electrode", "Sum Sq"]
  SS_residual    <- anova_table["Residuals", "Sum Sq"]
  
  # Extract df
  df_side        <- anova_table["Bin", "Df"]
  df_electrode   <- anova_table["Electrode", "Df"]
  df_interaction <- anova_table["Bin:Electrode", "Df"]
  df_residual    <- anova_table["Residuals", "Df"]
  
  # RMS values
  RMS_side        <- sqrt(SS_side / df_side)
  RMS_electrode   <- sqrt(SS_electrode / df_electrode)
  RMS_interaction <- sqrt(SS_interaction / df_interaction)
  RMS_noise       <- sqrt(SS_residual / df_residual)
  
  # CNR
  CNR <- RMS_interaction / RMS_noise
  
  return(CNR)
}


tva_files <- list.files("../../decoding_manuscript/Submission/Revision/trial_data_averages", full.names = T)

tva_list <- tva_files %>%
  set_names(basename(.) %>% str_remove("_avg\\.csv$")) %>% 
  map(~ {
    read.csv(.x)
  },
  .progress = TRUE
  )

cnr_results <- 
  tibble(
    id = names(tva_list),
    CNR = map_dbl(tva_list, compute_cnr) 
  ) %>% 
  left_join(d %>% select(id,visit_participant_group), by = "id") %>%  rename(group=visit_participant_group) %>% 
  mutate(CNR = round(CNR,2))
write.csv(cnr_results,"../../decoding_manuscript/Submission/Revision/CNR_results.csv",row.names = F)

# Incorporate with decoding data
cnr_results <- read_csv("../../decoding_manuscript/Submission/Revision/CNR_results.csv")

ai_ind = left_join(d_t1$all_images_age, cnr_results %>% select(id,CNR), by = "id")
fh_ind = left_join(d_t1$faces_houses_age, cnr_results %>% select(id,CNR), by = "id")
ui_ind = left_join(d_t1$upright_inverted_age, cnr_results %>% select(id,CNR), by = "id")
fi_ind = left_join(d_t1$face_identity_age, cnr_results %>% select(id,CNR), by = "id")

ai_avg <- ai_ind %>% group_by(id,group,CNR,age_cat) %>% summarize(mean_acc=mean(acc)) %>% ungroup() %>% mutate(cat="Image Selectivity")
fh_avg <- fh_ind %>% group_by(id,group,CNR,age_cat) %>% summarize(mean_acc=mean(acc)) %>% ungroup() %>% mutate(cat="Face Selectivity")
ui_avg <- ui_ind %>% group_by(id,group,CNR,age_cat) %>% summarize(mean_acc=mean(acc)) %>% ungroup() %>% mutate(cat="Orientation Selectivity")
fi_avg <- fi_ind %>% group_by(id,group,CNR,age_cat) %>% summarize(mean_acc=mean(acc)) %>% ungroup() %>% mutate(cat="Identity Selectivity")

# CNR differences by group and age
bind_rows(
  broom::tidy(anova(lm(CNR~group, data = ai_avg))) %>% mutate(type="all_images"),
  broom::tidy(anova(lm(CNR~group, data = fh_avg))) %>% mutate(type="faces_houses"),
  broom::tidy(anova(lm(CNR~group, data = ui_avg))) %>% mutate(type="upright_inverted"),
  broom::tidy(anova(lm(CNR~group, data = fi_avg))) %>% mutate(type="face_identity")
) %>% select(type, term,p.value) %>% 
  filter(term=="group")
bind_rows(
  broom::tidy(anova(lm(CNR~group+age_cat, data = ai_avg))) %>% mutate(type="all_images"),
  broom::tidy(anova(lm(CNR~group+age_cat, data = fh_avg))) %>% mutate(type="faces_houses"),
  broom::tidy(anova(lm(CNR~group+age_cat, data = ui_avg))) %>% mutate(type="upright_inverted"),
  broom::tidy(anova(lm(CNR~group+age_cat, data = fi_avg))) %>% mutate(type="face_identity")
) %>% select(type, term,p.value) %>% 
  filter(term=="age_cat")

# The impact of CNR on decoding accuracy
bind_rows(
  broom::tidy(anova(lm(mean_acc~CNR, data = ai_avg))) %>% mutate(type="all_images"),
  broom::tidy(anova(lm(mean_acc~CNR, data = fh_avg))) %>% mutate(type="faces_houses"),
  broom::tidy(anova(lm(mean_acc~CNR, data = ui_avg))) %>% mutate(type="upright_inverted"),
  broom::tidy(anova(lm(mean_acc~CNR, data = fi_avg))) %>% mutate(type="face_identity")
) %>% select(type,term,p.value) %>% 
  filter(term=="CNR")

bind_rows(
  broom::tidy(anova(lm(mean_acc~CNR*age_cat, data = ai_avg))) %>% mutate(type="all_images"),
  broom::tidy(anova(lm(mean_acc~CNR*age_cat, data = fh_avg))) %>% mutate(type="faces_houses"),
  broom::tidy(anova(lm(mean_acc~CNR*age_cat, data = ui_avg))) %>% mutate(type="upright_inverted"),
  broom::tidy(anova(lm(mean_acc~CNR*age_cat, data = fi_avg))) %>% mutate(type="face_identity")
) %>% select(type,term,p.value) %>% 
  filter(term=="CNR:age_cat")


all_avg = rbind(ai_avg,fh_avg, ui_avg,fi_avg) %>% 
  rename(`Age Group (years)` = age_cat)

all_avg$cat <- factor(all_avg$cat, levels = c("Image Selectivity", "Face Selectivity", "Orientation Selectivity", "Identity Selectivity")) 



p1e <- all_avg %>%
  ggplot(aes(x=CNR, y=mean_acc,color=`Age Group (years)`)) +
  facet_wrap(~cat, scales = "free",ncol=4)+
  geom_point(alpha=.6) +
  geom_smooth(method="lm") +
  scale_color_bmj()+
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white")) +
  xlab("Contrast-to-Noise Ratio (CNR)") +
  ylab("Decoding Accuracy") +
  scale_y_continuous(labels = scales::percent_format())
p1e
ggsave(p1e, file = "../../decoding_manuscript/Submission/Revision/_Figures/CNR_Supplement.png", dpi=300, height=3,width=12)

## analysis regarding IQ

iq_table <- bind_rows(
  broom.mixed::tidy(anova(lmer(acc~group * iq + (1|id), data = fh_ind))) %>% mutate(type="faces_houses"),
  broom.mixed::tidy(anova(lmer(acc~group * iq + (1|id), data = ui_ind))) %>% mutate(type="upright_inverted"),
  broom.mixed::tidy(anova(lmer(acc~group * iq + (1|id), data = fi_ind))) %>% mutate(type="face_identity"),
  broom.mixed::tidy(anova(lmer(acc~group * iq + (1|id), data = ai_ind))) %>% mutate(type="all_images")
  ) %>% 
  select(type,term,NumDF,DenDF,statistic, p.value) %>% 
  filter(term == "iq" | term =="group:iq") %>% 
  mutate_if(is.numeric, round,2)

write.csv(iq_table,file = "../../decoding_manuscript/Submission/Revision/_Figures/iq_on_acc.csv", row.names = F)


#create datasets with clinical variables
clin_vars <- c("diagnosis_summary_full_scale_best_iq.t1", "imported_srs_t_score.t1", "visit_ados2_comparison_score_ss.t1", "imported_vineland_2_vi3_soc_ss.t1", "imported_vineland_2_vi3_com_ss.t1", "@vineland_2_vi3_2scale_ss.t1", "imported_nepsy_calc_Mf_Scaled_Score.t1")
fh_df <- tidyr::pivot_wider(d_t1$faces_houses_metrics,names_from = metric,values_from = value) %>% 
  mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% 
  left_join(phenotype %>% select(id,clin_vars), by = "id") %>% 
  rename(iq = diagnosis_summary_full_scale_best_iq.t1,
         ados=visit_ados2_comparison_score_ss.t1,
         srs=imported_srs_t_score.t1,
         vab_abc = `@vineland_2_vi3_2scale_ss.t1`,
         vab_soc = imported_vineland_2_vi3_soc_ss.t1,
         vab_com = imported_vineland_2_vi3_com_ss.t1,
         face_memory = imported_nepsy_calc_Mf_Scaled_Score.t1)
ui_df <- tidyr::pivot_wider(d_t1$upright_inverted_metrics,names_from = metric,values_from = value) %>% 
  mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% 
  left_join(phenotype %>% select(id,clin_vars), by = "id") %>% 
  rename(iq = diagnosis_summary_full_scale_best_iq.t1,
         ados=visit_ados2_comparison_score_ss.t1,
         srs=imported_srs_t_score.t1,
         vab_abc = `@vineland_2_vi3_2scale_ss.t1`,
         vab_soc = imported_vineland_2_vi3_soc_ss.t1,
         vab_com = imported_vineland_2_vi3_com_ss.t1,
         face_memory = imported_nepsy_calc_Mf_Scaled_Score.t1)
fi_df <- tidyr::pivot_wider(d_t1$face_identity_metrics,names_from = metric,values_from = value) %>% 
  mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% 
  left_join(phenotype %>% select(id,clin_vars), by = "id") %>% 
  rename(iq = diagnosis_summary_full_scale_best_iq.t1,
         ados=visit_ados2_comparison_score_ss.t1,
         srs=imported_srs_t_score.t1,
         vab_abc = `@vineland_2_vi3_2scale_ss.t1`,
         vab_soc = imported_vineland_2_vi3_soc_ss.t1,
         vab_com = imported_vineland_2_vi3_com_ss.t1,
         face_memory = imported_nepsy_calc_Mf_Scaled_Score.t1)
ai_df <- tidyr::pivot_wider(d_t1$all_images_metrics,names_from = metric,values_from = value) %>% 
  mutate(decoding_sus = tidyr::replace_na(decoding_sus,0)) %>% 
  left_join(phenotype %>% select(id,clin_vars), by = "id") %>% 
  rename(iq = diagnosis_summary_full_scale_best_iq.t1,
         ados=visit_ados2_comparison_score_ss.t1,
         srs=imported_srs_t_score.t1,
         vab_abc = `@vineland_2_vi3_2scale_ss.t1`,
         vab_soc = imported_vineland_2_vi3_soc_ss.t1,
         vab_com = imported_vineland_2_vi3_com_ss.t1,
         face_memory = imported_nepsy_calc_Mf_Scaled_Score.t1)

iq_metrics <- 
  bind_rows(
    broom::tidy(summary(lm(decoding_onset~iq, data = fh_df))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
    broom::tidy(summary(lm(decoding_onset~iq, data = ui_df))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
    broom::tidy(summary(lm(decoding_onset~iq, data = fi_df))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
    broom::tidy(summary(lm(decoding_onset~iq, data = ai_df))) %>% mutate(metric= "decoding_onset", type = "all_images"),
    broom::tidy(summary(lm(decoding_sus ~ iq, data = fh_df))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
    broom::tidy(summary(lm(decoding_sus ~ iq, data = ui_df))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
    broom::tidy(summary(lm(decoding_sus ~ iq, data = fi_df))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
    broom::tidy(summary(lm(decoding_sus ~ iq, data = ai_df))) %>% mutate(metric= "decoding_sus", type = "all_images"),
    broom::tidy(summary(lm(decoding_peak ~ iq, data = fh_df))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
    broom::tidy(summary(lm(decoding_peak ~ iq, data = ui_df))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
    broom::tidy(summary(lm(decoding_peak ~ iq, data = fi_df))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
    broom::tidy(summary(lm(decoding_peak ~ iq, data = ai_df))) %>% mutate(metric= "decoding_peak", type = "all_images"),
    broom::tidy(summary(lm(decoding_latency ~ iq, data = fh_df))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
    broom::tidy(summary(lm(decoding_latency ~ iq, data = ui_df))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
    broom::tidy(summary(lm(decoding_latency ~ iq, data = fi_df))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
    broom::tidy(summary(lm(decoding_latency ~ iq, data = ai_df))) %>% mutate(metric= "decoding_latency", type = "all_images")
) %>%
  select(type,metric,term,statistic,p.value) %>% 
  filter(term=="iq") %>% 
  mutate_if(is.numeric, round,3)

write.csv(iq_metrics,file = "../../decoding_manuscript/Submission/Revision/_Figures/iq_on_metrics.csv", row.names = F) 

## evaluate associations with clinical phenotype.
tab3 <- bind_rows(
  broom::tidy(summary(lm(decoding_onset~ados,        data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_onset~srs,         data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_onset~vab_abc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_onset~vab_soc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_onset~vab_com,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_onset~face_memory, data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "faces_houses"),
  
  broom::tidy(summary(lm(decoding_onset~ados,        data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_onset~srs,         data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_onset~vab_abc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_onset~vab_soc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_onset~vab_com,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_onset~face_memory, data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "upright_inverted"),
  
  broom::tidy(summary(lm(decoding_onset~ados,        data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  broom::tidy(summary(lm(decoding_onset~srs,         data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  broom::tidy(summary(lm(decoding_onset~vab_abc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  broom::tidy(summary(lm(decoding_onset~vab_soc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  broom::tidy(summary(lm(decoding_onset~vab_com,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  broom::tidy(summary(lm(decoding_onset~face_memory, data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "face_identity"),
  
  broom::tidy(summary(lm(decoding_onset~ados,        data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  broom::tidy(summary(lm(decoding_onset~srs,         data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  broom::tidy(summary(lm(decoding_onset~vab_abc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  broom::tidy(summary(lm(decoding_onset~vab_soc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  broom::tidy(summary(lm(decoding_onset~vab_com,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  broom::tidy(summary(lm(decoding_onset~face_memory, data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_onset", type = "all_images"),
  ## decoding sus
  broom::tidy(summary(lm(decoding_sus~ados,        data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_sus~srs,         data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_sus~vab_abc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_sus~vab_soc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_sus~vab_com,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_sus~face_memory, data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "faces_houses"),
  
  broom::tidy(summary(lm(decoding_sus~ados,        data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_sus~srs,         data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_sus~vab_abc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_sus~vab_soc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_sus~vab_com,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_sus~face_memory, data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "upright_inverted"),
  
  broom::tidy(summary(lm(decoding_sus~ados,        data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  broom::tidy(summary(lm(decoding_sus~srs,         data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  broom::tidy(summary(lm(decoding_sus~vab_abc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  broom::tidy(summary(lm(decoding_sus~vab_soc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  broom::tidy(summary(lm(decoding_sus~vab_com,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  broom::tidy(summary(lm(decoding_sus~face_memory, data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "face_identity"),
  
  broom::tidy(summary(lm(decoding_sus~ados,        data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  broom::tidy(summary(lm(decoding_sus~srs,         data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  broom::tidy(summary(lm(decoding_sus~vab_abc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  broom::tidy(summary(lm(decoding_sus~vab_soc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  broom::tidy(summary(lm(decoding_sus~vab_com,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  broom::tidy(summary(lm(decoding_sus~face_memory, data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_sus", type = "all_images"),
  
  ## decoding peak
  broom::tidy(summary(lm(decoding_peak~ados,        data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_peak~srs,         data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_peak~vab_abc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_peak~vab_soc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_peak~vab_com,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_peak~face_memory, data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "faces_houses"),
  
  broom::tidy(summary(lm(decoding_peak~ados,        data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_peak~srs,         data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_peak~vab_abc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_peak~vab_soc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_peak~vab_com,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_peak~face_memory, data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "upright_inverted"),
  
  broom::tidy(summary(lm(decoding_peak~ados,        data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  broom::tidy(summary(lm(decoding_peak~srs,         data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  broom::tidy(summary(lm(decoding_peak~vab_abc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  broom::tidy(summary(lm(decoding_peak~vab_soc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  broom::tidy(summary(lm(decoding_peak~vab_com,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  broom::tidy(summary(lm(decoding_peak~face_memory, data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "face_identity"),
  
  broom::tidy(summary(lm(decoding_peak~ados,        data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  broom::tidy(summary(lm(decoding_peak~srs,         data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  broom::tidy(summary(lm(decoding_peak~vab_abc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  broom::tidy(summary(lm(decoding_peak~vab_soc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  broom::tidy(summary(lm(decoding_peak~vab_com,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  broom::tidy(summary(lm(decoding_peak~face_memory, data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_peak", type = "all_images"),
  
  ## decoding latency
  broom::tidy(summary(lm(decoding_latency~ados,        data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_latency~srs,         data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_latency~vab_abc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_latency~vab_soc,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_latency~vab_com,     data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  broom::tidy(summary(lm(decoding_latency~face_memory, data = fh_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "faces_houses"),
  
  broom::tidy(summary(lm(decoding_latency~ados,        data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_latency~srs,         data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_latency~vab_abc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_latency~vab_soc,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_latency~vab_com,     data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  broom::tidy(summary(lm(decoding_latency~face_memory, data = ui_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "upright_inverted"),
  
  broom::tidy(summary(lm(decoding_latency~ados,        data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  broom::tidy(summary(lm(decoding_latency~srs,         data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  broom::tidy(summary(lm(decoding_latency~vab_abc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  broom::tidy(summary(lm(decoding_latency~vab_soc,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  broom::tidy(summary(lm(decoding_latency~vab_com,     data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  broom::tidy(summary(lm(decoding_latency~face_memory, data = fi_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "face_identity"),
  
  broom::tidy(summary(lm(decoding_latency~ados,        data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images"),
  broom::tidy(summary(lm(decoding_latency~srs,         data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images"),
  broom::tidy(summary(lm(decoding_latency~vab_abc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images"),
  broom::tidy(summary(lm(decoding_latency~vab_soc,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images"),
  broom::tidy(summary(lm(decoding_latency~vab_com,     data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images"),
  broom::tidy(summary(lm(decoding_latency~face_memory, data = ai_df %>% filter(group =="asd")))) %>% mutate(metric= "decoding_latency", type = "all_images")
) %>%
  filter(term != "(Intercept)") %>% 
  select(type,metric, term, estimate, std.error, p.value) %>%  
  mutate_if(is.numeric,round,5)
tab3$type <- factor(tab3$type, levels = c("all_images", "faces_houses","face_identity","upright_inverted"))

col_order <- c(
  "type",
  "term",
  "decoding_onset_estimate","decoding_onset_std.error","decoding_onset_p.value",
  "decoding_sus_estimate","decoding_sus_std.error","decoding_sus_p.value",
  "decoding_peak_estimate","decoding_peak_std.error","decoding_peak_p.value",
  "decoding_latency_estimate","decoding_latency_std.error","decoding_latency_p.value")

tab3_indexed <- tab3 %>%
  group_by(type, metric) %>%
  mutate(row_id = row_number()) %>%
  ungroup()

tab3_wide <- tab3_indexed %>%
  pivot_wider(
    id_cols = c(type, term, row_id),
    names_from = metric,
    values_from = c(estimate, std.error, p.value),
    names_glue = "{metric}_{.value}"
  ) %>% select(all_of(col_order)) %>% 
  arrange(type)

write.csv(tab3_wide,file = "../../decoding_manuscript/Submission/Revision/_Figures/Table3.csv", row.names = F) 

##
fh_df %>% filter(group=="asd") %>% ggplot(aes(x=face_memory, y = decoding_sus,color=group)) + geom_point() +geom_smooth(method="lm")
fi_df %>% filter(group=="asd") %>% ggplot(aes(x=ados, y = decoding_sus,color=group)) + geom_point() +geom_smooth(method="lm")



