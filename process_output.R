process_output <- function(timepoint){
  #timepoint = "T1"
  demo_key <- read_csv("demos.csv",show_col_types = FALSE) %>% select(1,3,7,9) %>% rename(id = visit_individual, sex = visit_sex, group = visit_participant_group, age = abcct_demographic_form_age_at_evaluation) %>%
  mutate(id = str_remove(id,"^.")) %>% 
  mutate(id = str_replace(id,"^(..).", "\\1.")) %>%
    mutate(age = floor(age/365))


  # Section 1 (Category Bins): Face Selectivity: Faces vs Houses
  files = list.files(path = paste0("../1. MVPA/", timepoint,"/mvpc_output/category_bins/faces vs houses"), pattern = "*.csv", full.names = T)
  data_list <- map(files, ~ read_csv(., show_col_types = F) %>% mutate(id = tools::file_path_sans_ext(basename(.x))), .progress = T)
  ll <- length(seq(100,250,10))
  faces_houses <- bind_rows(data_list) %>% 
    mutate(id = str_extract(id, "^[^_]+")) %>%
    mutate(id = str_remove(id,"^.")) %>% 
    mutate(id = str_replace(id,"^(..).", "\\1.")) %>% 
    select(id,everything()) %>% 
    left_join(demo_key, by = "id") %>% 
    mutate(chance_level = .50) %>% 
    mutate(individual_chance = ifelse(acc-(se*2) > chance_level,"yes","no")) %>% 
    mutate(group = factor(group, levels = c("td", "asd")))
  #01.9182
  faces_houses_metrics <- bind_rows(
    faces_houses %>% filter(individual_chance == "yes" & time > 0) %>% 
      group_by(id,group,sex,age) %>%
      summarize(value = min(time), .groups = "drop") %>% 
      mutate(metric="decoding_onset"),
    faces_houses %>% filter(individual_chance == "yes" & time >=100 & time <=250) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = n()/ll) %>% 
      mutate(metric = "decoding_sus"),
    faces_houses %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = max(acc),.groups = "drop") %>% 
      mutate(metric = "decoding_peak"),
    faces_houses %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>%
      slice(which.max(acc)) %>% 
      ungroup() %>% 
      mutate(value = time) %>%
      select(id,group,sex,age,value) %>% 
      mutate(metric = "decoding_latency")
  ) %>% mutate(age_cat = case_when(age == 6 ~ "6-7",age == 7 ~ "6-7",age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) 

  avg_faces_houses <- faces_houses %>%
    group_by(group,time) %>% 
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .50)
  
  avg_faces_houses_age <- faces_houses %>%
    mutate(age_cat = case_when(age == 6 ~ "6-7",
                               age == 7 ~ "6-7",
                               age == 8 ~ "8-9",
                               age == 9 ~ "8-9",
                               age == 10 ~ "10-11",
                               age == 11 ~ "10-11")) %>% 
    group_by(age_cat, group,time) %>% 
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .50) %>%
    mutate(age_cat=factor(age_cat, levels = c("10-11", "8-9", "6-7"))) 
  
  
  # Section 2 (Category Bins): Orientation Selectivity: Upright vs. Inverted
  files = list.files(path = paste0("../1. MVPA/", timepoint, "/mvpc_output/category_bins/upright vs inverted"), pattern = "*.csv", full.names = T)
  #as_tibble(t(read_delim(files[1], col_names=F))) %>% rename(time=V1,acc=V2)
  data_list <- map(files, ~ read_csv(., show_col_types = F) %>% mutate(id = tools::file_path_sans_ext(basename(.x))))
  upright_inverted <- bind_rows(data_list) %>%
    mutate(id = str_extract(id, "^[^_]+")) %>%
    mutate(id = str_remove(id,"^.")) %>% 
    mutate(id = str_replace(id,"^(..).", "\\1.")) %>% 
    select(id,everything()) %>% 
    left_join(demo_key, by = "id") %>% 
    mutate(chance_level = .50) %>% 
    mutate(individual_chance = ifelse(acc-(se*2) > chance_level,"yes","no")) %>% 
    mutate(group = factor(group, levels = c("td", "asd")))
  
  upright_inverted_metrics <- bind_rows(
    upright_inverted %>% filter(individual_chance == "yes" & time > 0) %>% 
      group_by(id,group,sex,age) %>%
      summarize(value = min(time), .groups = "drop") %>% 
      mutate(metric="decoding_onset"),
    upright_inverted %>% filter(individual_chance == "yes" & time >=100 & time <=250) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = n()/ll,.groups = "drop") %>% 
      mutate(metric = "decoding_sus"),
    upright_inverted %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = max(acc),.groups = "drop") %>% 
      mutate(metric = "decoding_peak"),
    upright_inverted %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>%
      slice(which.max(acc)) %>% 
      ungroup() %>% 
      mutate(value = time) %>%
      select(id,group,sex,age,value) %>% 
      mutate(metric = "decoding_latency")
  ) %>% mutate(age_cat = case_when(age == 6 ~ "6-7",age == 7 ~ "6-7",age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) 
  
  avg_upright_inverted <- upright_inverted %>%
    group_by(group,time) %>% 
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .50)
  
  avg_upright_inverted_age <- upright_inverted %>%
    mutate(age_cat = case_when(age == 6 ~ "6-7",
                               age == 7 ~ "6-7",
                               age == 8 ~ "8-9",
                               age == 9 ~ "8-9",
                               age == 10 ~ "10-11",
                               age == 11 ~ "10-11")) %>% 
    group_by(age_cat, group,time) %>%  
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .50) %>%
    mutate(age_cat=factor(age_cat, levels = c("10-11", "8-9", "6-7"))) 
  
  
  
  # Section 3 (Image Bins): Identity: CauU vs. afrU vs AsiU
  files = list.files(path = paste0("../1. MVPA/", timepoint, "/mvpc_output/image_bins/face_identity"), pattern = "*.csv", full.names = T)
  #as_tibble(t(read_delim(files[1], col_names=F))) %>% rename(time=V1,acc=V2)
  data_list <- map(files, ~ read_csv(., show_col_types = F) %>% mutate(id = tools::file_path_sans_ext(basename(.x))))
  
  face_identity <- bind_rows(data_list) %>%
    mutate(id = str_remove(id,"_mvpc")) %>% 
    select(id,everything()) %>% 
    left_join(demo_key, by = "id") %>% 
    mutate(chance_level = .33) %>% 
    mutate(individual_chance = ifelse(acc-(se*2) > chance_level,"yes","no")) %>% 
    mutate(group = factor(group, levels = c("td", "asd")))
  
  face_identity_metrics <- bind_rows(
    face_identity %>% filter(individual_chance == "yes" & time > 0) %>% 
      group_by(id,group,sex,age) %>%
      summarize(value = min(time), .groups = "drop") %>% 
      mutate(metric="decoding_onset"),
    face_identity %>% filter(individual_chance == "yes" & time >=100 & time <=250) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = n()/ll,.groups = "drop") %>% 
      mutate(metric = "decoding_sus"),
    face_identity %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = max(acc),.groups = "drop") %>% 
      mutate(metric = "decoding_peak"),
    face_identity %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>%
      slice(which.max(acc)) %>% 
      ungroup() %>% 
      mutate(value = time) %>%
      select(id,group,sex,age,value) %>% 
      mutate(metric = "decoding_latency")
  ) %>% mutate(age_cat = case_when(age == 6 ~ "6-7",age == 7 ~ "6-7",age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) 
  
  avg_face_identity <- face_identity %>%
    group_by(group,time) %>% 
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .33)
  
  avg_face_identity_age <- face_identity %>%
    mutate(age_cat = case_when(age == 6 ~ "6-7",
                               age == 7 ~ "6-7",
                               age == 8 ~ "8-9",
                               age == 9 ~ "8-9",
                               age == 10 ~ "10-11",
                               age == 11 ~ "10-11")) %>% 
    group_by(age_cat, group,time) %>% 
    summarise(mean_acc = mean(acc),
              sd = sd(acc),
              se = sd/sqrt(n()),
              .groups = "drop") %>% 
    mutate(chance_level = .33) %>%
    mutate(age_cat=factor(age_cat, levels = c("10-11", "8-9", "6-7"))) 
  
  
  # Section 4 (Image Bins): All images
  files = list.files(path = paste0("../1. MVPA/", timepoint, "/mvpc_output/image_bins/all_bins"), pattern = "*.csv", full.names = T)
  #as_tibble(t(read_delim(files[1], col_names=F))) %>% rename(time=V1,acc=V2)
  data_list <- map(files, ~ read_csv(., show_col_types = F) %>% mutate(id = tools::file_path_sans_ext(basename(.x))))
 
  all_images <- bind_rows(data_list) %>%
   mutate(id = str_remove(id,"_mvpc")) %>% 
   select(id,everything()) %>% 
   left_join(demo_key, by = "id") %>% 
   mutate(chance_level = .11) %>% 
   mutate(individual_chance = ifelse(acc-(se*2) > chance_level,"yes","no")) %>% 
    mutate(group = factor(group, levels = c("td", "asd")))
 
  all_images_metrics <- bind_rows(
    all_images %>% filter(individual_chance == "yes" & time > 0) %>% 
      group_by(id,group,sex,age) %>%
      summarize(value = min(time), .groups = "drop") %>% 
      mutate(metric="decoding_onset"),
    all_images %>% filter(individual_chance == "yes" & time >=100 & time <=250) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = n()/ll,.groups = "drop") %>% 
      mutate(metric = "decoding_sus"),
    all_images %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>% 
      summarize(value = max(acc),.groups = "drop") %>% 
      mutate(metric = "decoding_peak"),
    all_images %>% filter(individual_chance == "yes" & time > 0) %>%
      group_by(id,group,sex,age) %>%
      slice(which.max(acc)) %>% 
      ungroup() %>% 
      mutate(value = time) %>%
      select(id,group,sex,age,value) %>% 
      mutate(metric = "decoding_latency")
  ) %>% mutate(age_cat = case_when(age == 6 ~ "6-7",age == 7 ~ "6-7",age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) 
  
 avg_all_images <- all_images %>%
   group_by(group,time) %>% 
   summarise(mean_acc = mean(acc),
             sd = sd(acc),
             se = sd/sqrt(n()),
             .groups = "drop") %>% 
   mutate(chance_level = .11)
 
 avg_all_images_age <- all_images %>%
   mutate(age_cat = case_when(age == 6 ~ "6-7",
                              age == 7 ~ "6-7",
                              age == 8 ~ "8-9",
                              age == 9 ~ "8-9",
                              age == 10 ~ "10-11",
                              age == 11 ~ "10-11")) %>% 
   group_by(age_cat, group,time) %>% 
   summarise(mean_acc = mean(acc),
             sd = sd(acc),
             se = sd/sqrt(n()),
             .groups = "drop") %>% 
   mutate(chance_level = .11) %>% 
   mutate(age_cat=factor(age_cat, levels = c("10-11", "8-9", "6-7"))) 
 
 ##prepare developmental age group data
 
 faces_houses_age <- faces_houses         %>% mutate(age_cat = case_when(age == 6 ~ "6-7", age == 7 ~ "6-7", age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) %>% mutate(age_cat=factor(age_cat, levels = c("6-7","8-9","10-11"))) %>% mutate(group=factor(group, levels = c("td","asd")))
 upright_inverted_age <- upright_inverted %>% mutate(age_cat = case_when(age == 6 ~ "6-7", age == 7 ~ "6-7", age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) %>% mutate(age_cat=factor(age_cat, levels = c("6-7","8-9","10-11"))) %>% mutate(group=factor(group, levels = c("td","asd")))
 face_identity_age <- face_identity       %>% mutate(age_cat = case_when(age == 6 ~ "6-7", age == 7 ~ "6-7", age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) %>% mutate(age_cat=factor(age_cat, levels = c("6-7","8-9","10-11"))) %>% mutate(group=factor(group, levels = c("td","asd")))
 all_images_age <- all_images             %>% mutate(age_cat = case_when(age == 6 ~ "6-7", age == 7 ~ "6-7", age == 8 ~ "8-9",age == 9 ~ "8-9",age == 10 ~ "10-11",age == 11 ~ "10-11")) %>% mutate(age_cat=factor(age_cat, levels = c("6-7","8-9","10-11"))) %>% mutate(group=factor(group, levels = c("td","asd")))

 

 ## function for running tests
 input=all_images_age
 run_tests <- function(input){
   time_values <- unique(input$time);
   model_summaries = list();
   model_summaries_cat = list();
   
   for(time_val in time_values){
     subset_data = input %>% filter(time == time_val)
     
     if("age_cat" %in% colnames(input)){
       lm_model = lm(acc~group*age,data = subset_data)
       lm_model_cat = lm(acc~group*age_cat, data = subset_data)
       model_summary_cat = broom::tidy(anova(lm_model_cat))
       model_summaries_cat[[as.character(time_val)]] =  model_summary_cat
     } else {
       lm_model = lm(acc~group,data = subset_data)
     }
     
     model_summary = broom::tidy(lm_model)
     model_summaries[[as.character(time_val)]] =  model_summary
     
     
   }
   # save the results
   if("age_cat" %in% colnames(input)){
     tests = bind_rows(model_summaries, .id = "time") %>% dplyr::filter(term == "groupasd:age")
   } else {
     tests = bind_rows(model_summaries, .id = "time") %>% dplyr::filter(term == "groupasd")
   }
   
   tests$adj_p.value <- p.adjust(tests$p.value, method = "fdr") 
   tests$time <- as.numeric(tests$time)
   
   #return_df = left_join(input, tests, by = "time")
   return(tests)
 }
 
 # run_tests(faces_houses) %>%
 #   group_by(time) %>% 
 #   filter(row_number() == 1) %>% 
 #   ungroup() %>% 
 #   mutate(type = "faces_houses") %>% select(-group,-id,-sex),
  

 ## statistical test data
 test_data <- bind_rows(
   run_tests(faces_houses) %>%
     mutate(type = "faces_houses"),
   run_tests(upright_inverted) %>%
     mutate(type = "upright_inverted"),
   run_tests(face_identity) %>%
     mutate(type = "face_identity"),
   run_tests(all_images) %>%
     mutate(type = "all_images")
   )
 
 
 test_data_age <- 
   bind_rows(
     run_tests(faces_houses_age) %>%
       mutate(type = "faces_houses"),
     run_tests(upright_inverted_age) %>%
       mutate(type = "upright_inverted"),
     run_tests(face_identity_age) %>%
       mutate(type = "face_identity"),
     run_tests(all_images_age) %>%
       mutate(type = "all_images")
     )
 
 

 # create one average data for list
 
 avg_data <- bind_rows(avg_faces_houses %>% mutate(type = "faces_houses"),
                       avg_upright_inverted %>% mutate(type = "upright_inverted"),
                       avg_face_identity %>% mutate(type = "face_identity"),
                       avg_all_images %>% mutate(type = "all_images")) %>% 
   mutate(chance = ifelse(mean_acc-(se*2) > chance_level,"yes","no"))
 
 
 avg_data_age <- bind_rows(avg_faces_houses_age %>% mutate(type = "faces_houses"),
                           avg_upright_inverted_age %>% mutate(type = "upright_inverted"),
                           avg_face_identity_age %>% mutate(type = "face_identity"),
                           avg_all_images_age %>% mutate(type = "all_images")) %>%
   mutate(chance = ifelse(mean_acc-(se*2) > chance_level,"yes","no"))
 
  return_list = list(faces_houses = faces_houses,
                     upright_inverted = upright_inverted,
                     face_identity = face_identity,
                     all_images = all_images,
                     
                     faces_houses_age = faces_houses_age,
                     upright_inverted_age = upright_inverted_age,
                     face_identity_age = face_identity_age,
                     all_images_age = all_images_age,
                     
                     faces_houses_metrics = faces_houses_metrics,
                     upright_inverted_metrics = upright_inverted_metrics,
                     face_identity_metrics = face_identity_metrics,
                     all_images_metrics = all_images_metrics,
                     
                     test_data = test_data,
                     test_data_age = test_data_age,
                     
                     avg_data = avg_data,
                     avg_data_age =  avg_data_age)
  return(return_list)
  
  
  
  
}

  
