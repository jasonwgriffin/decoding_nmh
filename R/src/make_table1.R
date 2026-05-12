make_table1 <- function(fpath, d){
  
  fpath = fig_path
  d=d
  jb = c('visit_participant_group','age_cat')
  d = d %>% 
    filter(id %in% t1_ids) %>% 
    mutate(age = floor(abcct_demographic_form_age_at_evaluation/365)) %>% 
    mutate(age_cat = case_when(age == 6 ~ "a6-7",
                               age == 7 ~ "a6-7",
                               age == 8 ~ "b8-9",
                               age == 9 ~ "b8-9",
                               age == 10 ~ "c10-11",
                               age == 11 ~ "c10-11"))
  
  #N
  t1 <- d %>% group_by(visit_participant_group,age_cat) %>% summarize(n = n(), .groups = "drop") %>% ungroup()
  
  #sex
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group,visit_sex,age_cat) %>% summarize(male_count = n(), .groups = "drop") %>% filter(visit_sex=="male") %>% ungroup(), by = jb) %>% 
    mutate(male_percentage = round(male_count/n*100,2))
  #age
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group,age_cat) %>% summarize(mean_age = round(mean(abcct_demographic_form_age_at_evaluation)/365,2), sd_age = round(sd(abcct_demographic_form_age_at_evaluation)/365,2), .groups = "drop") %>% ungroup(), by = jb)
  
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, abcct_demographic_form_Participant_Ethnicity,age_cat) %>% summarize(count = n(), .groups = "drop") %>% ungroup() %>% 
                           pivot_wider(
                             names_from = abcct_demographic_form_Participant_Ethnicity,
                             values_from = count,
                             values_fill = 0   # fills missing race counts with 0
                           ), by = jb) %>%
    mutate(hispanic_cell = paste0(`hispanic-or-latino`,' (', round(`hispanic-or-latino`/n*100,2),"%)"),
           nonhispanic_cell = paste0(`not-hispanic-or-latino`,' (', round(`not-hispanic-or-latino`/n*100,2),"%)")
    )
  
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, abcct_demographic_form_Participant_Race,age_cat) %>% summarize(count = n(), .groups = "drop") %>% ungroup() %>% 
    pivot_wider(
      names_from = abcct_demographic_form_Participant_Race,
      values_from = count,
      values_fill = 0   # fills missing race counts with 0
    ), by = jb) %>%
    mutate(black_cell = paste0(`african-american-or-black`,' (', round(`african-american-or-black`/n*100,2),"%)"),
           white_cell = paste0(white,' (', round(white/n*100,2),"%)"),
           #alaska_cell = paste0(`american-indian-or-alaska-native`,' (', round(`american-indian-or-alaska-native`/n*100,2),"%)"),
           asian_cell = paste0(asian,' (', round(asian/n*100,2),"%)"),
           other_cell = paste0(other,' (', round(other/n*100,2),"%)"),
           mixed_cell = paste0(`mixed-race`,' (', round(`mixed-race`/n*100,2),"%)")
    )
  

  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, age_cat) %>% summarize(mean_iq = round(mean(diagnosis_summary_full_scale_best_iq.t1),2),    sd_iq = round(sd(diagnosis_summary_full_scale_best_iq.t1),2), .groups = "drop"), by=jb)
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, age_cat) %>% summarize(mean_ados = round(mean(visit_ados2_comparison_score_ss.t1),2),       sd_ados = round(sd(visit_ados2_comparison_score_ss.t1),2), .groups = "drop"), by=jb)
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, age_cat) %>% summarize(mean_srs = round(mean(imported_srs_t_score.t1,na.rm=T),2),           sd_srs = round(sd(imported_srs_t_score.t1,na.rm=T),2), .groups = "drop"), by=jb)
  t1 <- t1 %>% left_join(d %>% group_by(visit_participant_group, age_cat) %>% summarize(mean_vab = round(mean(imported_vineland_2_vi3_abc_ss.t1,na.rm=T),2), sd_vab = round(sd(imported_vineland_2_vi3_abc_ss.t1,na.rm=T),2), .groups = "drop"), by=jb)
  
  t1 <- t1 %>% 
    mutate(iq_cell = paste0(mean_iq, ' (', sd_iq, ')')) %>% 
    mutate(age_cell = paste0(mean_age, ' (', sd_age, ')')) %>% 
    mutate(sex_cell = paste0(male_count, ' (', male_percentage, '%)')) %>%
    mutate(ados_cell = paste0(mean_ados, ' (', sd_ados, ')')) %>%
    mutate(srs_cell = paste0(mean_srs, ' (', sd_srs, ')')) %>% 
    mutate(vab_cell = paste0(mean_vab, ' (', sd_vab, ')'))
  
  
  
  t1 <- t1 %>% select(age_cat,visit_participant_group,n,sex_cell, age_cell, white_cell,mixed_cell,black_cell,asian_cell,other_cell, nonhispanic_cell, hispanic_cell, iq_cell,ados_cell,srs_cell,vab_cell)
  t1 <- t1 %>% arrange(age_cat,visit_participant_group)
  table1 <- t(t1)
  
  
  write.csv(table1, file = paste0(fpath,'/table1.csv'), row.names = T)
  
}
