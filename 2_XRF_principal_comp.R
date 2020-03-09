# # # # # #
# Starting to process XRF data
# Jordi F. Pagès
# 07-10-2019
# University of Barcelona
# # # # # # 

# Libraries
library(tidyverse)
library(tidylog)
library(cowplot)
library(magick)
library(robCompositions)
library(ggsci)
library(RColorBrewer)

# Data input. 
# Go to 1_XRF_processing.R and choose the file_name to be loaded. 
     #
   # # #
 # # ! # # 
# # # # # # 
     #
     #
# Save changes to 1_XRF_processing.R before running source(1_XRF_processing.R)
source("1_XRF_processing.R")

# File_name to be used for printing figures
file_name <- "ErmsTancada"
# file_name <- "Llanada"
# file_name <- "PanissosBIS"
# file_name <- "Secanella"
# file_name <- "Encanyissada"



# # # # # # # # # # # # # # # # # # # # # # # 
# 0) Plots with counts for each element ---- 
# # # # # # # # # # # # # # # # # # # # # # # 

ggplot(XRFdata_tidy_filtered) +
  geom_line(aes(x = -Sample, y = Area, col = Section), lwd = 0.3) +
  xlab("Depth (mm)") +
  ylab("Counts per second") +
  coord_flip() +
  ggtitle(label = file_name) +
  facet_wrap(~Element, scales = "free") +
  theme_bw() +
  theme(legend.position = 'none')
# ggsave2(filename = paste("Figures/", file_name, "_counts_per_element.pdf", sep = ""))




# # # # # # # # # # # # # # # # # # # # # # # 
# 1) Preparing the data set for the PCA ---- 
# # # # # # # # # # # # # # # # # # # # # # # 

# We check if there are any zeros or negative values. Compositional data must be positive,
# according to the book Applied Compositional Data Analysis, by Filzmoser et al. 2018.
XRFdata_wide <- XRFdata_tidy_filtered %>% 
  mutate(Area = ifelse(Area < 100, NA, Area)) %>% 
  select(-Chi, -Std) %>% 
  pivot_wider(names_from = Element, values_from = Area) %>%  # This makes the data set wide again, to compute PCA more easily
  na.omit()

# We will add the minimum value + 1 to all the samples. This doesn't affect anything, because of the 
# scale invariance property of compositional data (Filzmoser et al. 2018, Book), and in this way
# we get rid of zeros and negative values
# XRFdata_wide <- XRFdata_tidy_filtered %>%
  # mutate(Area2 = Area + abs(min(Area)) +1) %>%  # This gets rid of zero and negative values in our data set.
  # select(-Chi, -Std, -Area) %>%
  # pivot_wider(names_from = Element, values_from = Area2) # This makes the data set wide again, to compute PCA more easily

# To check that our results give exactly the same using the _tidy_filtered or _wide version
# of the data set we plot Aluminium, for example.
# XRFdata_tidy_filtered %>% 
#   filter(Element == "Al") %>% 
#   ggplot() +
#   geom_line(aes(x = Sample, y = Area))
# 
# XRFdata_wide %>% 
#   select(Al, Sample) %>% 
#   ggplot() +
#   geom_line(aes(x = Sample, y = Al))
# All good. Playing with the the data set, doesn't change anything (apart from the y axis,
# which has been increased due to the addition of min(Area), to get rid of zeros and negative values)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 2) Running a robust PCA using pcaCoDa() from robCompositions ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# pcaCoDa uses isometric log ratios, which are the best transformation to work with compositional data.
# However, the space of ilr transformation is difficult to interpret, therefore,  pcaCoDa transforms the results from the 
# PCA, back to centre-log-ratio space, which is easier to interpret (equivalent to the original compositions space).
# Moreover, pcaCoDa uses a robust estimator of the covariance matrix, which is not influenced by the presence of outliers.
# Since covariance matrix using the robust estimator uses a random number generator, to make results reproducible, we
# use the function set.seed()

# robCompositions uses data.frames instead of tibbles. Let's use convert the dataset into a data.frame
XRFdata_wide_df <- XRFdata_wide %>% 
  select(-CoreID, -Section, -Sample, -Red, -Green, -Blue, -a, -l, -b)
XRFdata_wide_df <- as.data.frame(XRFdata_wide_df)
rownames(XRFdata_wide_df) <- XRFdata_wide$Sample

# Compositional PCA, non-robust, for comparison.
p_comp <- pcaCoDa(XRFdata_wide_df, method = "classical")

# Optional, if we include colour_data: We use the mult_comp argument because in the XRFdata_wide_df we have columns 1:3 that 
# are from 1 composition (colour data, LAB, if you have L and A you can get B),
# and columns 4:16 is another composition (the "concentrations" of the elements)
# p_comp <- pcaCoDa(XRFdata_wide_df, method = "classical", mult_comp = list(c(1:3, 4:length(names(XRFdata_wide_df)))))

# Compositional PCA, robust. The good one.
set.seed(234)
p_comp_rob <- pcaCoDa(XRFdata_wide_df, method = "robust")
# p_comp_rob <- pcaCoDa(XRFdata_wide_df, method = "robust", mult_comp = list(c(1:3, 4:length(names(XRFdata_wide_df)))))

summary(p_comp_rob)
plot(p_comp_rob, type = "l") # 3 components seem to be needed.

# We extract the % variance explained
variance_explained <- p_comp_rob$eigenvalues/sum(p_comp_rob$eigenvalues)

# We choose the number of components using a cutoff: a PC will be chosen if it explains 
# more than what 2 variables would explain by chance (2/number of columns)
components <- filter(enframe(variance_explained), variance_explained > 2/length(names(XRFdata_wide_df)))
AutVal <- length(components$name)

par(mfrow = c(1,2))
biplot(p_comp, xlabs = rownames(XRFdata_wide_df))
biplot(p_comp_rob, xlabs = rownames(XRFdata_wide_df))



# # # # # # # # # # # # # # # # # # # # # # # # 
# 3) Making beautiful biplots with ggplot ---- 
# # # # # # # # # # # # # # # # # # # # # # # #

# Using code from Carmen Leiva.

# We first extract the loadings for the components of interest from p_comp_rob
dbloadings <- as.data.frame(p_comp_rob$loadings[,1:AutVal])

# Then we extract the scores for the components of interest from p_comp_rob
dbscores <- as.data.frame(p_comp_rob$scores[,1:AutVal])
dbscores$Sample <- XRFdata_wide$Sample
dbscores$Section <- as.factor(XRFdata_wide$Section)

# Biplot for PC1 vs PC2
ggplot(data = dbloadings)+
  geom_point(data = dbscores, aes(x = Comp.1, y = Comp.2, colour = Section)) + 
  # coord_fixed(ratio = 1) +
  annotate(geom = "text", 
           x = (dbloadings$Comp.1*3.4), 
           y = (dbloadings$Comp.2*3.4), 
           label = rownames(dbloadings), 
           col = "black") +
  xlab(paste("PC1 (", round(variance_explained[1]*100), "%)", sep = "")) +
  ylab(paste("PC2 (", round(variance_explained[2]*100), "%)", sep = "")) +
  geom_segment(data = dbloadings, 
               aes(x =0 , y = 0, xend = Comp.1*3, yend = Comp.2*3), 
               arrow = arrow(length = unit(1/2, 'picas')), 
               color = "black", linetype = 'solid', size = 0.5) +
  ggtitle(label = file_name) +
  scale_color_d3(palette = "category20") +
  theme_bw(base_size = 20,
           base_line_size = 0.5) +
  theme(legend.position="none")
# ggsave2(filename = paste("Figures/", file_name, "_biplot_PC1vsPC2.pdf", sep = ""))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 4) Making heatmap plot of the loadings for each element ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

names(dbloadings) <- c("PC1", "PC2")
dbloadings$Elements <- rownames(dbloadings)

dbloadings %>%
  pivot_longer(cols = starts_with("PC"), names_to = "components", values_to = "loadings") %>% 
  ggplot(aes(x = Elements, y = components)) +
  geom_tile(aes(fill = loadings)) +
  scale_fill_gradient2(low = "#2467AD", high = "#B2182B", mid = "white", midpoint = 0) +
  xlab("") +
  ylab("") +
  ggtitle(label = file_name) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 17,
           base_line_size = 0.5)
# ggsave2(filename = paste("Figures/", file_name, "_loadings_heatmap.pdf", sep = ""))



# # # # # # # # # # # # # # # # # # # # # # # # # # 
# 5) Making ggplots of the pca scores downcore ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # #

names(dbscores) <- c("PC1", "PC2", "Sample", "Section")

dbscores_plot <- dbscores %>%
  left_join(XRFdata_wide, by = "Sample") %>% 
  select(starts_with("PC"), "Sample", "Section.x", "CoreID", "l") %>% 
  pivot_longer(cols = c(starts_with("PC"), l), names_to = "Components", values_to = "Scores") %>%
  mutate(Sample = Sample-min(Sample), # to make the plot start at 0 cm, approximately where the 1st XRF measurement was taken
         Sample_cm = Sample/10,
         labs = ifelse(Components == "PC1", 
                       paste("PC1 (", round(variance_explained[1]*100), "%)", sep = ""),
                       ifelse(Components == "PC2", 
                              paste("PC2 (", round(variance_explained[2]*100), "%)", sep = ""),
                              "Core lightness"))) %>% 
  mutate(labs = factor(labs, 
                       levels = c(paste("PC1 (", round(variance_explained[1]*100), "%)", sep = ""),
                                  paste("PC2 (", round(variance_explained[2]*100), "%)", sep = ""),
                                  "Core lightness")))

# We find the points where we change section using diff. We divide by 10 to get it in cm instead of mm.
section_break <- dbscores_plot$Sample[which(diff(as.numeric(dbscores_plot$Section.x)) == 1)]/10

p <- ggplot(dbscores_plot, aes(x = -Sample_cm, y = Scores)) +
  # geom_point(aes(colour = Components)) +
  geom_line(colour = "#757575", lwd = 0.3) +
  geom_smooth(aes(colour = Components, fill = Components), span = 0.08) +
  scale_colour_manual(values = c("#1E6DA8","#FE7223", "#2C9431")) + 
  scale_fill_manual(values = c("#1E6DA8","#FE7223", "#2C9431")) +
  geom_vline(aes(xintercept = -section_break), lty = 2, colour = "darkgrey") +
  # geom_vline(aes(xintercept = -section_break[1]), lty = 2, colour = "darkgrey") + # For PanissosBIS
  # geom_vline(aes(xintercept = -section_break[2]), lty = 2, colour = "darkgrey") + # For PanissosBIS
  # scale_color_d3(palette = "category20") +
  xlab("Core depth (cm)") +
  ylab("") +
  facet_wrap(~labs, scales = "free") +
  ggtitle(label = file_name) +
  coord_flip() +
  theme_bw(base_size = 17,
           base_line_size = 0.5) +
  theme(legend.position = "none")

ggdraw(p) + 
  draw_label("Scores (clr-robust)", x = 0.215, y = 0.03) +
  draw_label("Scores (clr-robust)", x = 0.54, y = 0.03) +
  draw_label("Lightness (%)", x = 0.865, y = 0.03)
# ggsave2(filename = paste("Figures/", file_name, "_scores_downcore.pdf", sep = ""))



# 
# core_img <- image_read(path = "Data/Images_processed/ErmsTancada1i2_joined_cropped.jpg") %>% 
#   image_rotate(degrees = 90) #%>% 
#   # image_resize("1800x5000!")
# 
# ggdraw() + 
#   draw_plot(p) +
#   draw_image(core_img, y = 0.13, height = 0.75, width = 1.765)
# ggsave2(filename = "Figures/ErmsTancada_scores_downcoreTEST.pdf")
