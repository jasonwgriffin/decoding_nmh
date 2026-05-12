make_figure3 <- function(fpath){
  #fpath = "~/Dropbox/Decoding/decoding_manuscript/Submission/Figures/"
  # mats
  mats = R.matlab::readMat("../2. RSA/RDMs/RSA_all_bins.mat")
  td = mats$mats[[1]][[1]]
  asd = mats$mats[[2]][[1]]
  # do you want this to be a dissimilairty matrix?
  td = 1-abs(td)
  asd = 1-abs(asd)
  #remove diagnoal
  diag(td) <- NA
  diag(asd) <- NA
  # Convert correlation matrix to long format
  df <- bind_rows(melt(td) %>% mutate(group="NT"), melt(asd) %>% mutate(group="Autistic")) %>%
    mutate(Var1 = case_when(Var1 == 1 ~ "House 1",
                          Var1 == 2 ~ "House 2",
                          Var1 == 3 ~ "House 3",
                          Var1 == 4 ~ "Caucasian U",
                          Var1 == 5 ~ "African U",
                          Var1 == 6 ~ "Asian U",
                          Var1 == 7 ~ "Caucasian I",
                          Var1 == 8 ~ "African I",
                          Var1 == 9 ~ "Asian I"),
         Var2 = case_when(Var2 == 1 ~ "House 1",
                          Var2 == 2 ~ "House 2",
                          Var2 == 3 ~ "House 3",
                          Var2 == 4 ~ "Caucasian U",
                          Var2 == 5 ~ "African U",
                          Var2 == 6 ~ "Asian U",
                          Var2 == 7 ~ "Caucasian I",
                          Var2 == 8 ~ "African I",
                          Var2 == 9 ~ "Asian I"))
  df$Var1 <- factor(df$Var1, levels = c("House 1", "House 2", "House 3", "African I", "Asian I", "Caucasian I", "African U", "Asian U", "Caucasian U"))
  df$Var2 <- factor(df$Var2, levels = c("House 1", "House 2", "House 3", "African I", "Asian I", "Caucasian I", "African U", "Asian U", "Caucasian U"))

# Rename columns
  gg_aes = list(
    geom_tile(color = "white"),
    
    #scale_fill_steps2(na.value = 'white', show.limits = T, limits = c(.05,1)),
    #scale_fill_gradient2(low = "yellow", high = "red", na.value = NA, limits = c(0,.3)),
    
    #scale_fill_viridis_c(option = "F"),
    #scale_fill_gradientn(
    #  colours = colorRampPalette(c('#709AE1', 'white', '#FD7446'))(100),
    #  breaks=seq(0,1,.1),
    #  limits=c(0, 1)
    #),
    #scale_fill_gradient2(low = "blue", mid = "white", high = "red",na.value = NA),
    scale_fill_distiller(palette = "BrBG", direction = 1),
    scale_x_discrete(limits = rev(levels(df$Var1)), position = "top"),
    scale_y_discrete(limits = rev(levels(df$Var2)), position = "left"),
    labs(fill="1-r"),
    theme_classic(),
    theme(plot.background = element_rect(fill="white",color="white"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(.5,.1,.1,.1, "cm"),
          plot.title = element_text(hjust = 1),
          aspect.ratio = 1)
    )
  

  p_td <- ggplot(df %>% filter(group=="NT"), aes(Var1, Var2, fill = value)) + gg_aes
  p_td
  legg <- cowplot::get_legend(p_td)
  #legg = cowplot::get_plot_component(p_td, 'guide-box-top', return_all = FALSE)
  
  p_td <- p_td +theme(legend.position = "none")
  pimage_x <- axis_canvas(p_td, axis = "x") + draw_image("../src/stim_images/production_cropped/array1.png",x = 4.5, scale=8.95)
  pimage_y <- axis_canvas(p_td, axis = "y") + draw_image("../src/stim_images/production_cropped/array1v.png",y = 4.5, scale=8.95)
  td_1 <- insert_xaxis_grob(p_td,pimage_x, position = "center")
  td_2 <- insert_yaxis_grob(td_1,pimage_y, position = "center")   
  td_p <- ggdraw(td_2)

  p_asd <- ggplot(df %>% filter(group=="Autistic"), aes(Var1, Var2, fill = value)) + gg_aes 
  p_asd <- p_asd +theme(legend.position = "none")
  pimage_x2 <- axis_canvas(p_asd, axis = "x") + draw_image("../src/stim_images/production_cropped/array1.png",x = 4.5, scale=8.95)
  pimage_y2 <- axis_canvas(p_asd, axis = "y") + draw_image("../src/stim_images/production_cropped/array1v.png",y = 4.5, scale=8.95)
  asd_1 <- insert_xaxis_grob(p_asd,pimage_x2, position = "center")
  asd_2 <- insert_yaxis_grob(asd_1,pimage_y2, position = "center")   
  asd_p <- ggdraw(asd_2)
  p <- plot_grid(td_p,
               asd_p,
               legg,
               ncol=3,
               rel_widths = c(1,1,.2))
  p
  plot_grid(p,legg,ncol=2,rel_widths = c(1,.1))
  
  
  ## MDS
 
  image_cols = c( 
                 "#6699ff", # House 1
                 "#7f66ff", # House 2
                 "#961aff", # House 3
                 "#ffcc33", # Caucasian u
                 "#ff6633", # African U
                 "#ff9333", # Asian U
                 "#ff0066", # Caucasian I
                 "#b40dbf", # African I
                 "#e2079f" # Asian I
                 )
  image_cols <- pal_npg(alpha = 1)(9)
  #create dissimilarity matrix, aka, 1-correlation
  
  #dm_td <- as.dist(1-td)
  #dm_asd <- as.dist(1-asd)
  dm_td <- as.dist(td)
  dm_asd <- as.dist(asd)
  
  # Perform Multidimensional Scaling (MDS)
  mds_td <- cmdscale(dm_td, k = 2)
  mds_asd <- cmdscale(dm_asd, k = 2)
  
  
  mds_df <- 
    bind_rows(as_tibble(mds_td, .name_repair = "unique") %>% rename(d1 = "...1", d2 = "...2") %>% mutate(group="NT"),
              as_tibble(mds_asd, .name_repair = "unique") %>% rename(d1 = "...1", d2 = "...2") %>% mutate(group="Autistic")
    )
  
  mds_df <- mds_df %>% 
    group_by(group) %>% 
    mutate(img = row_number()) %>% 
    ungroup() %>% 
    mutate(stim = case_when(img == 1 ~ "House 1",
                            img == 2 ~ "House 2",
                            img == 3 ~ "House 3",
                            img == 4 ~ "Caucasian U",
                            img == 5 ~ "African U",
                            img == 6 ~ "Asian U",
                            img == 7 ~ "Caucasian I",
                            img == 8 ~ "African I",
                            img == 9 ~ "Asian I")) 
  
  write_csv(mds_df,paste0("../2. RSA/RDMs/mds_df.csv"))
  
  mds_df <- mds_df %>%
    mutate(path = paste0("../src/stim_images/production_rsa_faces/",stim, ".jpg"))
  
  mds_df$stim <- factor(mds_df$stim, levels = c("House 1", "House 2", "House 3", "African I", "Asian I", "Caucasian I", "African U", "Asian U", "Caucasian U"))
  
  mds_df <- mds_df %>% 
    mutate(bin = case_when(str_detect(stim, "House") ~ "Houses",
                           str_detect(stim, " U") ~ "Upright Faces",
                           str_detect(stim, " I") ~ "Inverted Faces"),
           bin_color = case_when(str_detect(stim, "House") ~ "red",
                                 str_detect(stim, " U") ~ "blue",
                                 str_detect(stim, " I") ~ "orange")
           )
  
  mds_df$bin <- factor(mds_df$bin, levels = c("Houses", "Upright Faces", "Inverted Faces"))
  #mds_df <- mds_df %>% mutate(path = "../stim_images/array1.png")

  library(ggtext)
  library(ggimage)
  axis_df <- tibble(
    label = c("***PC2***", "***PC1***"),
    x=c(.02,.6),
    y=c(.6,-.03),
    fill = c("grey", "red"),
    angle=c(90,0)
  )
  bin_cols = c("red", "blue", "#440d92")
  bin_cols = pal_bmj("default")(3)
  
  
  mds_df <- mds_df %>% mutate(cropped_imgs = circle_crop(path, border_size = 50, border_colour = c(rep(bin_cols, each=3))))
  #library(cropcircles)
  
  mds_aes = list(
    geom_hline(yintercept = 0,linewidth=1, color="gray93"),
    geom_vline(xintercept = 0,linewidth=1, color="gray93"),
    geom_polygon(aes(group=bin,color = bin), fill = NA, linewidth=2,show.legend = FALSE),
    geom_point(aes(color=bin)),
    
    geom_image(size=.10),
    scale_color_bmj(),
    #scale_color_manual(values = bin_cols),
    theme_half_open(12),
    facet_wrap(~group),
    #facet_grid(bin~group),
    coord_cartesian(xlim = c(-.6,.6), ylim = c(-.6,.6)),
    geom_richtext(inherit.aes = F,data = axis_df, aes(x=x, y=y, label = label, fill = fill, angle=angle),
                  fill=NA,label.color=NA,size=2.5),
    guides(color = guide_legend(override.aes = list(size = 6, shape = 16))),
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text =  
        element_textbox(
        size = 12,
        color = "black", fill = "grey", box.color = "grey",
        halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
        padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)
        )
    )
  )

  ggdraw(legg)
  #ggplot(mds_df %>% filter(group == "asd", bin == "Face"), aes(x=d1,y=d2, image=path, color = bins)) + mds_aes
  asd_mds_p <- ggplot(mds_df %>% filter(group == "Autistic"), aes(x=d1,y=d2,label=stim,image=cropped_imgs)) + mds_aes + 
    annotation_custom(
      grob = legg,
      xmin = -0.72,
      xmax = -0.82,
      ymin = 0.4,
      ymax = 0.5
    )
  td_mds_p <- ggplot(mds_df %>% filter(group == "NT") %>% mutate(group = str_replace_all(group, "NT", "Neurotypical")), aes(x=d1,y=d2,label=stim,image=cropped_imgs)) + mds_aes
  
  mds_p = plot_grid(asd_mds_p, td_mds_p)

  asd_3 = plot_grid(asd_p, asd_mds_p, rel_widths = c(1,1.2))
  td_3 = plot_grid(td_p, td_mds_p, rel_widths = c(1,1.2))
  
  p_final <- plot_grid(asd_3,td_3,ncol=1,labels = c("A", "B"))
  
  
  
  save_plot(p_final, file=paste0(fpath, "Figure 3.png"),dpi=300,base_height=8,base_width=11,bg = "white")
  

}
