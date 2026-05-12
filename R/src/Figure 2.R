make_figure3 <- function(fpath,h,w,run_label){
  
  df = d_t1$avg_data %>% mutate(timepoint = "T1") %>% rename(Time=time)
  diff_df <- d_t1$test_data %>% rename(Time=time) %>% mutate(sig_group_diff = ifelse(p.value < .051,"yes","no")) %>%
    filter(term=="group")
  df1 <- left_join(df,diff_df %>% select(Time,sig_group_diff, type), by = c("Time", "type"))
  
  types = unique(df$type)
  chance_levels = c(.50,.50,.33,.11)
  y_coord_min = c(.4,.4,.31,.1)
  #y_coord_max = c(.8,.8,.42,.15)
  titles = c("Face Selectivity", "Orientation Selectivity", "Identity Selectivity", "Image Selectivity")
  ybreaks = data.frame(faces = c(.4,.5,.6,.7,.8),
                       orientation = c(.4,.5,.6,.7,.8),
                       identity = c(.31,.33,.35,.37,.39),
                       image = c(.1,.11,.12,.13,.14))
  
  
  label_scaling = c(.95,.95,.975,.97)
  
  y_coord_max = c(ybreaks[5,1] + ybreaks[5,1]-ybreaks[4,1],
                  ybreaks[5,2] + ybreaks[5,2]-ybreaks[4,2],
                  ybreaks[5,3] + ybreaks[5,3]-ybreaks[4,3],
                  ybreaks[5,4] + ybreaks[5,4]-ybreaks[4,4])
  p_grid = list();
  #group.labs <- c("Neurotypical", "Autism")
  #names(group.labs) <- c("td", "asd")
  
  df$group <- str_replace(df$group, "td", "NT")
  df$group <- str_replace(df$group, "asd", "Autistic")
  
  
  #chance_rects = c(.4,.4,.30+.0025,.10)
  chance_rects = c(.4,.4,.31,.10)
  test_diff_rects = c(ybreaks[5,1],ybreaks[5,2],ybreaks[5,3],ybreaks[5,4])
  #chance_rects_widths = c(.02,.02,.01,.005)
  chance_rects_widths = c((ybreaks[2,1]-ybreaks[1,1])*.25,
                          (ybreaks[2,2]-ybreaks[1,2])*.25,
                          (ybreaks[2,3]-ybreaks[1,3])*.25*.5,
                          (ybreaks[2,4]-ybreaks[1,4])*.25)

  #xmin = c(200,200,200,200)
  #xmax = c(300,300,300,400)
  
  
  i=0
  for(i in 1:length(types)){
    rs=chance_rects[i]
    diff_rs=test_diff_rects[i]
    
    chance_diff = df %>% filter(type ==types[i]) %>% filter(chance == "yes",Time > 0)
    chance_rect <- data.frame(xmin = c(as.numeric(chance_diff$Time)-5), xmax = c(as.numeric(chance_diff$Time)), ymin = c(rs), ymax = c(rs + chance_rects_widths[i]), timepoint = chance_diff$timepoint, group = chance_diff$group)
    
    test_diff = diff_df %>% filter(type ==types[i]) %>% filter(sig_group_diff == "yes",Time > 50)
    test_diff_rect = data.frame(xmin=c(as.numeric(test_diff$Time)-5), xmax = c(as.numeric(test_diff$Time)), ymin = c(diff_rs), ymax = c(diff_rs + chance_rects_widths[i]*.8))


    p = 
      ggplot(
      df %>% filter(type == types[i]) %>%
        mutate(group=factor(group, levels = c("NT","Autistic")))
    ) +
      geom_hline(yintercept = chance_levels[i], color = "grey") +
      #geom_label(colour = "white", fontface = "bold", label = "Chance Level Accuracy", x = 375, y = chance_levels[i], size = 2, fill="black", label.padding = unit(0.15, "lines")) +
      geom_ribbon(aes(x=Time,y=mean_acc,color=group,fill=group,ymin=mean_acc-se, ymax=mean_acc+se), alpha=.4,color=NA) +
      geom_line(inherit.aes=F, show.legend = T,aes(x=Time,y=mean_acc,color=group)) +
      geom_rect(inherit.aes = F,data = chance_rect, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=group), position = ggstance::position_dodgev(height=chance_rects_widths[i]), alpha=.7) +
      #colors
      #scale_color_iterm() +
      #scale_fill_iterm() +
      scale_color_manual(values = c("#FF666DB2", "#00CDE8B2")) +
      scale_fill_manual(values = c("#FF666DB2", "#00CDE8B2")) +
      
      scale_x_continuous(expand = c(0,0),breaks=c(seq(-200,500,100))) +
      scale_y_continuous(breaks = ybreaks[[i]],
                         labels = scales::percent_format())+
      
      coord_cartesian(ylim = c(y_coord_min[i], y_coord_max[i]),
                      xlim = c(-200,500))+
      # category top right text
      annotate(
        "label", 
        x = -190,  # Slightly less than max x for padding
        y = y_coord_max[i]*label_scaling[i],  # Slightly less than max y for padding
        #label = bquote(italic(.(titles[i]))), 
        label = titles[i], 
        hjust = 0,  # Align text to the right
        vjust = 0,  # Align text to the top
        size = 2.5,
        color = "black"
      ) +
      

      #ggtitle("Hello") +
      
      theme_cowplot() + 
      theme(plot.background = element_rect(fill = "white"), 
            axis.title = element_text(size=10),
            axis.title.x = element_text(size=8, vjust=5),
            axis.text.x = element_text(size=6, angle=45),
            #legend.position = c(.02,.88),
            legend.position = "inside",
            legend.position.inside = c(.02, .81),
            legend.direction = "vertical",
            legend.title = element_blank(),
            plot.margin = margin(1,.5,0,.5, "cm"),
            #plot.title = element_text(size=15, hjust = 0, vjust = 0),
            legend.text = element_text(margin = margin(l=0), size=6),
            legend.key.width = unit(.3,"cm"),
            legend.key.size=unit(.2,"lines")) +
      
      #facet_wrap(~timepoint) +
                 #labeller = labeller(group = group.labs)) +
      
      guides(color = guide_legend(override.aes = list(linewidth=2.5)),
             fill = "none") +
      ylab("Decoding Accuracy")
    
    #add brackets if needed
    if(types[i] != "face_identity"){
      p = p + geom_bracket(xmin = test_diff %>% slice_min(Time) %>% pull(Time)-15, xmax = test_diff %>% slice_max(Time) %>% pull(Time)+10, y.position = ybreaks[5,i] + chance_rects_widths[i]*1.5, label = "Autistic < NT", label.size = 2,tip.length = c(.02,0)) +
        geom_rect(inherit.aes = F,data = test_diff_rect, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=.8,fill = "grey")  
    }
    
    p_grid[[i]] = p
  }

  a1 = "../src/stim_images/production_cropped/array1.png"
  a2 = "../src/stim_images/production_cropped/array2.png"
  a3 = "../src/stim_images/production_cropped/array3.png"
  a4 = "../src/stim_images/production_cropped/array4.png"

  p1 = ggdraw() + draw_plot(p_grid[[4]]) +                                                      draw_image(a1,x = 1, y = 1, hjust = 1.10, vjust = .72, halign =  0,valign = .6, width = .65)
  p2 = ggdraw() + draw_plot(p_grid[[1]]) +  draw_image(a2,x = 1, y = 1, hjust = 1.10, vjust = .72, halign =  0,valign = .6, width = .65)
  p3 = ggdraw() + draw_plot(p_grid[[3]]) +                                                      draw_image(a3,x = 1, y = 1, hjust = 1.10, vjust = .72, halign =  0,valign = .6, width = .65)
  p4 = ggdraw() + draw_plot(p_grid[[2]]) +  draw_image(a4,x = 1, y = 1, hjust = 1.10, vjust = .72, halign =  0,valign = .6, width = .65)
  p = plot_grid(p1,p2,p3,p4,ncol=2,
                labels = c("A", "B", "C", "D"))
  
  p
  
  #ggsave(p1,file=paste0(fpath,"fig2a.png"),dpi=300,height=3,width=4)
  #ggsave(p2,file=paste0(fpath,"fig2b.png"),dpi=300,height=3,width=4)
  #ggsave(p3,file=paste0(fpath,"fig2c.png"),dpi=300,height=3,width=4)
  #ggsave(p4,file=paste0(fpath,"fig2d.png"),dpi=300,height=3,width=4.5)
  
  ggsave(p, file = paste0("../1. MVPA/T1/",run_label, "/xfigures/Figure 3.png"), dpi=300, height=h,width=w)
  #ggsave(p, file = paste0(fpath, "Figure 3.png"), dpi=300, height=h,width=w)
  ggsave(p, file = paste0(fpath, "Figure 3.pdf"), dpi=300, height=h,width=w)
  

  
}

