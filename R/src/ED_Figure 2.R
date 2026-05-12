make_figure5 <- function(fpath,td_matrix, asd_matrix){
  #fpath = "~/Dropbox/Decoding/decoding_manuscript/Submission/Figures/"
  
  dfs <- c("faces_houses", "upright_inverted", "face_identity", "all_bins", "category_bins")
  

  
  pull_tg_matrix <- function(fp,mp){
    df = read.csv(fp, header = F)
    rownames(df) = 1:nrow(df)
    colnames(df) = 1:ncol(df)
    df = df %>% mutate(row = rownames(df)) %>% 
      tidyr::pivot_longer(cols = -row, names_to = "col", values_to = "value") %>% 
      mutate(group = ifelse(str_detect(fp,'asd'),'asd','td'))
    
    mask = read.csv(mp, header = F)
    rownames(mask) = 1:nrow(mask)
    colnames(mask) = 1:ncol(mask)
    mask = mask %>% mutate(row = rownames(mask)) %>% 
      tidyr::pivot_longer(cols = -row, names_to = "col", values_to = "significant") %>% 
      mutate(group = ifelse(str_detect(fp,'asd'),'asd','td'))
    
    df_mask <- full_join(df,mask, by=c("col","row","group"))
    
  }
  
  
  i=0
  plot_list <- vector("list", 5)
  for(i in 1:length(dfs)){
    td_matrix = paste0("src/tg_dataframes/", dfs[i], "_td.csv")
    asd_matrix = paste0("src/tg_dataframes/", dfs[i], "_asd.csv")
    mask_matrix = paste0("src/tg_dataframes/", dfs[i], "_mask.csv")
    
    df = bind_rows(pull_tg_matrix(td_matrix,mask_matrix),pull_tg_matrix(asd_matrix,mask_matrix))
    
    df$row <- as.numeric(df$row)*10-200
    df$col <- as.numeric(df$col)*10-200
    
    df <- df %>% mutate(group = case_when(group =="td"~ "NT",
                                    group == "asd"~"Autistic"))
    
    #display.brewer.pal(n = 8, name = 'BrBG')
    #brewer.pal(n = 8, name = "BrBG")
    #"#8C510A" "#BF812D" "#DFC27D" "#F6E8C3" "#C7EAE5" "#80CDC1" "#35978F" "#01665E"
    #library(viridisLite)
    #magma_cols <- viridis(8, option = "magma")
    #magma_cols"#000004FF" "#231151FF" "#5F187FFF" "#982D80FF" "#D3436EFF" "#F8765CFF" "#FEBA80FF" "#FCFDBFFF"
    mask <- df %>%
      mutate(label = "Difference") %>% 
      ggplot(aes(x=row,y=col,fill=as.factor(significant))) +
      geom_tile() +
      geom_hex(color=scales::alpha("white", 0.1)) +
      geom_vline(xintercept = 0,color="white", linewidth=2) +
      geom_hline(yintercept = 0,color="white", linewidth=2) +
      
      #scale_fill_viridis_d(option="A") +
      
      #scale_fill_distiller(palette = "BrBG") + 
      scale_fill_manual(values = c("#231151FF", "#FEBA80FF"))+
      
      scale_x_continuous(expand=c(0,0), breaks = seq(-200,500,100)) +
      scale_y_continuous(expand=c(0,0), breaks = seq(-200,500,100), position = "right") +
      coord_cartesian(xlim=c(5,500), ylim = c(5,500)) +
      theme_classic() +
      ylab(NULL) +
      xlab("Generalization Time") +
      facet_wrap(~label) + 
      theme(legend.position = "none",
            #axis.text.y = element_blank(),
            plot.margin = margin(1.5,0.5,0,0, "cm"),
            strip.background = element_blank(),
            strip.text = element_textbox(
              size = 12,
              color = "black", fill = "grey", box.color = "grey",
              halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
              padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)))
    mask  
    p <- df %>% 
      ggplot(aes(x=row,y=col,fill=value)) + 
      geom_tile()+
      #geom_raster(aes(fill = value), interpolate = TRUE) +
      geom_vline(xintercept = 0,color="white", linewidth=2) +
      geom_hline(yintercept = 0,color="white", linewidth=2) +
      scale_fill_viridis_c(option="A") +
      #scale_fill_distiller(palette = "BrBG", direction = 1) + 
      scale_color_viridis_c(option="A") +
      #scale_fill_gsea() +
      #scale_fill_bmj() +
      #scale_y_continuous(labels = ) +
      #scale_x_continuous(labels = number_format(accuracy = 0.01)) +
      scale_x_continuous(expand=c(0,0), breaks = seq(-200,500,100)) +
      scale_y_continuous(expand=c(0,0), breaks = seq(-200,500,100)) +
      #scale_fill_continuous(value = label_number(accuracy = 0.01, trim = TRUE)) +
      coord_cartesian(xlim=c(5,500), ylim = c(5,500)) +
      ylab("Testing Time") +
      xlab("Generalization Time") +
      facet_wrap(~group) +
      theme_classic() +
      theme(legend.position = "right",
            legend.title = element_blank(),
            legend.direction = "vertical",
            legend.key.height = unit(1.25, 'cm'),  # Increase legend key height
            legend.key.width = unit(.2,'cm'),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margin around legend box
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),  
            #legend.spacing.y = unit(0.9, 'cm'),  # Add ver
            legend.text = element_text(size=10),
            plot.margin = margin(1.5,.5,0,.5, "cm"),
            plot.background = element_rect(fill = "white"),
            strip.background = element_blank(),
            strip.text = element_textbox(
              size = 12,
              color = "black", fill = "grey", box.color = "grey",
              halign = 0.5, linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
              padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3)
            ))
    
    p
    # get legend for grid
    p_leg = get_legend(p)
    # remove legend from panel
    
    
    
    pgrid = plot_grid(p, mask,ncol=2,rel_widths = c(1,.45))
    pgrid
    
    plot_list[[i]] = pgrid
  }
    
  
  a1 = "../src/stim_images/production_cropped/array1.png"
  a2 = "../src/stim_images/production_cropped/array2.png"
  a3 = "../src/stim_images/production_cropped/array3.png"
  a4 = "../src/stim_images/production_cropped/array4.png"
  p1 = ggdraw() + draw_plot(plot_list[[4]]) + draw_image(a1,x = .08, y = 1, hjust = 0, vjust = .16, halign =  0,valign = 0, width = .5)
  p2 = ggdraw() + draw_plot(plot_list[[1]]) + draw_image(a2,x = .08, y = 1, hjust = 0, vjust = .15, halign =  0,valign = 0, width = .5)
  p3 = ggdraw() + draw_plot(plot_list[[3]]) + draw_image(a3,x = .08, y = 1, hjust = 0, vjust = .15, halign =  0,valign = 0, width = .5)
  p4 = ggdraw() + draw_plot(plot_list[[2]]) + draw_image(a4,x = .08, y = 1, hjust = 0, vjust = .15, halign =  0,valign = 0, width = .5)
  
  p_final = plot_grid(p1 +  theme(plot.margin = margin(t = 20, r = 0, b = 0, l = 0),
                                  plot.background = element_rect(fill = "white")),
                      p2,
                      p3,
                      p4,nrow=4,
                      labels = c("A", "B", "C", "D"))
  
  #ggsave(p_final, file = paste0(fpath, "Figure 5.png"),dpi=300,height=15,width=10)
  ggsave(p_final, file = paste0(fpath, "Figure 5.pdf"),dpi=300,height=15,width=10)
  
  
}

