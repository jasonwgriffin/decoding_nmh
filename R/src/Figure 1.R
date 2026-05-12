library(png)
library(ggplot2)
library(ggpubr)
library(tibble)
library(scales)
library(viridisLite)
show_col(viridis_pal(option="A")(20))
# Import the image
#setwd("~/Dropbox/asd_brain_map")
img <- png::readPNG("src/src/brain_schem_tp.png")

df <- readxl::read_excel("src/src/brain_atlas_ASD.xlsx")
df$area <- factor(df$area)
df <- df %>% dplyr::distinct(area,.keep_all = T)
df$area <- stringr::str_to_sentence(df$area)
df$area <- stringr::str_to_title(df$area)

## add alpha based on sub_lobe2


#df <- df %>% arrange(desc(face_network)) %>%
# Reorder column "Name" with the custom order followed by alphabetical sorting

custom_order <- c("Fusiform Gyrus", "Superior Temporal Sulcus", "Inferior Occipital Gyrus")
  
df <- df %>%
  mutate(area = factor(area, levels = c(custom_order, setdiff(sort(unique(area)), custom_order)))) %>%
  arrange(area)



pt_size = 8;
shape = 21;
c_stroke = "black";
legend_text_size = 12;
stroke_size = 1.5;
legend_p_size = 7;

# p1 <- 
#   ggplot(df, aes(area,x=x,y=y, fill = sub_lobe, color = area, alpha = division)) +
#   background_image(img)+
#   scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
#   scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
#   #scale_alpha_manual(values=c(.33,.66,.99), guide = "none") +
#   geom_point(stroke = stroke_size, size=pt_size,color = c_stroke, shape=shape) +
#   guides(fill = guide_legend(override.aes = list(size=legend_p_size), ncol=1),
#          alpha = "none") +
#   scale_alpha_manual(values = c(1,.2)) +
#   #scale_shape_manual(values = c(21,22,23,24,25)) +
#   scale_fill_viridis_d(option = "A") +
#   scale_color_viridis_d(option = "A") +
#   #cowplot::theme_cowplot()+
#   theme(
#     legend.text=element_text(size=legend_text_size),
#     legend.margin=margin(l=-15,t=-30),
#     axis.title = element_blank( ),
#     axis.ticks = element_blank(),
#     axis.text = element_blank(),
#     plot.background = element_rect(fill = "white"),
#     legend.title = element_blank(),
#     panel.background = element_rect(fill="white"),
#     legend.key = element_rect(fill = NA),
#     legend.background = element_rect(fill="white")
#     #legend.key.size = unit(.2, 'cm')
#   )
# p1

p1 <-
  ggplot(df, aes(area,x=x,y=y, fill = area, alpha = cortical_depth)) +
  background_image(img)+
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
  scale_alpha_manual(values=c(.33,.66,.99), guide = "none") +
  geom_point(shape=shape,color=c_stroke,stroke = stroke_size, size=pt_size) +
  guides(fill = guide_legend(override.aes = list(size=legend_p_size), ncol=2)) +
  scale_fill_viridis_d(option = "A") +
  #scale_color_viridis_d(option = "A") +
  #cowplot::theme_cowplot()+
  theme(
    legend.text=element_text(size=legend_text_size),
    legend.margin=margin(l=-15,t=-30),
    axis.title = element_blank( ),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    panel.background = element_rect(fill="white"),
    legend.key = element_rect(fill = NA),
    legend.background = element_rect(fill="white")
    #legend.key.size = unit(.2, 'cm')
  )
p1
#sts_cols <- c("#5f187f","#982d80","#feba80")



# p2 <- p2 + annotate("label", x = 76,y=60, label = "Face Perception:\n Core system",size=7,
#                     fill="cornflowerblue",color = "black")
fpath <- "~/Dropbox/UH/_Research/Publications/Decoding/decoding_manuscript/Submission/production/Figures/"
ggsave(p2, file = "brain_fig_sts_stream.png", dpi=300,height=7,width = 12)

ggsave(p1, file = paste0(fpath, "Figure 1.pdf"), dpi=300,height=7,width = 12)
ggsave(p1, file = paste0(fpath, "Figure 1.pdf"), dpi=300,height=7,width = 12)

#ggsave(p1_final, file = "brain_fig_final.tiff", dpi=300,height=7,width = 12)