make_figure6 <- function(fpath,h,w,run_label){
  
  types = unique(d_t1$avg_data_age$type)
  chance_levels = c(.50,.50,.33,.11)
  y_coord_min = c(.4,.4,.28,.1)
  y_coord_max = c(.8,.8,.48,.18)
  titles = c("Face Selectivity", "Orientation Selectivity", "Identity Selectivity", "Image Selectivity")
  p_grid = list();
  group.labs <- c("NT", "Autistic")
  names(group.labs) <- c("td", "asd")
  
  chance_rects = c(.4,.4,.28,.10)
  chance_rects_widths = c(.02,.02,.02/2,.005)
  i=0
  for(i in 1:length(types)){
    rs=chance_rects[i]
    chance_diff = d_t1$avg_data_age %>% filter(type ==types[i]) %>% filter(chance == "yes" & time > 0)
    chance_rect <- data.frame(xmin = c(as.numeric(chance_diff$time)-5), xmax = c(as.numeric(chance_diff$time)), ymin = c(rs), ymax = c(rs + chance_rects_widths[i]), age_cat = chance_diff$age_cat, group = chance_diff$group)
    
    p = ggplot(
      d_t1$avg_data_age %>% filter(type == types[i]) %>% mutate(group=factor(group, levels = c("td","asd")))) +
      geom_hline(yintercept = chance_levels[i], color = "grey") +
      #geom_label(colour = "white", fontface = "bold", label = "Chance Level Accuracy", x = 375, y = chance_levels[i], size = 2, fill="black", label.padding = unit(0.15, "lines")) +
      geom_ribbon(aes(x=time,y=mean_acc,color=age_cat,fill=age_cat,ymin=mean_acc-se, ymax=mean_acc+se), alpha=.4,color=NA) +
      geom_line(inherit.aes=F, show.legend = T,aes(x=time,y=mean_acc,color=age_cat)) +
      geom_rect(inherit.aes = F,data = chance_rect, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=age_cat), position = ggstance::position_dodgev(height=chance_rects_widths[i]), alpha=.7) +
      
      facet_wrap(~group,
                 labeller = labeller(group = group.labs)) +
      
      scale_color_iterm() +
      scale_fill_iterm() +
      scale_x_continuous(breaks=c(seq(-200,500,100))) +
      scale_y_continuous(labels = scales::percent_format())+
      
      coord_cartesian(ylim = c(y_coord_min[i], y_coord_max[i]))+
      
      geom_label(
        data = data.frame(
          group = factor("asd", levels = c("td","asd")),   # keep facet order
          x = 499,
          y = y_coord_max[i] * 1.01,
          label = titles[i]
        ),
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 1,
        size = 2.5,
        color = "black"
      ) +
      
      theme_cowplot() + 
      theme(plot.background = element_rect(fill = "white"), 
            axis.title = element_text(size=12),
            axis.text.x = element_text(size=6, angle=45),
            legend.position = c(.01,.9),
            legend.title = element_blank(),
            plot.margin = margin(1,.5,0,.5, "cm"),
            #plot.title = element_text(size=30),
            legend.text = element_text(size = 8),
            legend.key.width = unit(.5,"cm"),
            legend.key.size=unit(.2,"lines")) +

      
      guides(color = guide_legend(override.aes = list(linewidth=2.5)),
             fill = FALSE) +
      ylab("Decoding Accuracy") +
      xlab("Time")
    
    p_grid[[i]] = p
  }
  
  a1 = "../src/stim_images/production_cropped/array1.png"
  a2 = "../src/stim_images/production_cropped/array2.png"
  a3 = "../src/stim_images/production_cropped/array3.png"
  a4 = "../src/stim_images/production_cropped/array4.png"
  
  p1 = ggdraw() + draw_plot(p_grid[[4]]) +                                                      draw_image(a1,x = 1, y = 1, hjust = 1.1, vjust = .65, halign =  0,valign = .6, width = .7)
  p2 = ggdraw() + draw_plot(p_grid[[1]] + theme(axis.title.y = element_text(color="white"))) +  draw_image(a2,x = 1, y = 1, hjust = 1.1, vjust = .65, halign =  0,valign = .6, width = .7)
  p3 = ggdraw() + draw_plot(p_grid[[3]]) +                                                      draw_image(a3,x = 1, y = 1, hjust = 1.1, vjust = .65, halign =  0,valign = .6, width = .7)
  p4 = ggdraw() + draw_plot(p_grid[[2]] + theme(axis.title.y = element_text(color="white"))) +  draw_image(a4,x = 1, y = 1, hjust = 1.1, vjust = .65, halign =  0,valign = .6, width = .7)
  p = plot_grid(p1,p2,p3,p4,ncol=2, labels = c("A", "B", "C", "D"))
  
  ggsave(p, file = paste0(fpath, "Figure 4.pdf"), dpi=300, height=h,width=w)
  
  }
