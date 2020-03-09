# # # # # #
# Starting to process XRF data
# Jordi F. Pagès
# 30-09-2019
# University of Barcelona
# # # # # # 

# Libraries
library(tidyverse)
library(tidylog)
library(cowplot)
library(magick)

# Auxiliary functions
source("Functions_for_XRF.R")

# XRF data input
file_name <- "ERMSTANCADA1&2_edited.csv"
# file_name <- "LLANADA1&2_edited.csv"
# file_name <- "PANISSOS1BIS,2BIS,3BISa_edited.csv"
# file_name <- "SECANELLA1&2a_edited.csv"
# file_name <- "ENCANYISSADA1a_edited.csv"

XRFdata <- read_csv(file = paste("Data/XRF_data/", file_name, sep = ""))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 1) Checking data quality - Coefficients of variation using replicates ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# 3 replicates were taken every 30 samples during XRF core-scanning
index <- duplicated(XRFdata$Sample) | duplicated(XRFdata$Sample, fromLast = TRUE)
XRFdata_replicates <- XRFdata[which(index == T),]

# Now we calculate the CV (which is just the mean/sd (see cv function in Functions_for_XRF.R file))
discarded_by_cv_replic <- XRFdata_replicates %>% 
  select(-CoreID, -Section) %>%
  select(Sample, ends_with("Area")) %>% 
  group_by(Sample) %>%
  summarise_all(cv) %>% 
  pivot_longer(cols = -Sample, 
               names_to = "Element") %>% 
  group_by(Element) %>% 
  summarise(mean_cv = mean(value)) %>% 
  filter(mean_cv>15.2 | mean_cv<0) %>% 
  mutate(Element_short = str_extract(Element, pattern = "[^-]*")) %>% 
  select(Element_short)
  
# The following elements had either negative counts, or CV>10%, both of which is problematic:
# ErmsTancada: Ar, As, Au, Cr, Cu, Ga, Mn, Mo, Nb, P, Pb, Rh, Rh, Se, V
# Llanada:     Ar, As, Cl, Cr, Mo, Nb, P, Pb, Rh, Rh, S, Se, V
# PanissosBIS: Ar, As, Cr, Cu, Mo, P, Rh, S,  Se, V
# Secanella:   Ar, As, Cr, Mo, Nb, P, Pb, Rh, S, Se, V
# Encanyissada: NOTHING... this core does not have replicates.



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 2) Checking data quality - Coefficients of variation for each sample ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# We delete coherence/incoherence data to simplify working with the data set. 
# REMEMBER TO CHECK COH-INCOH ratio afterwards!
XRFdata <- XRFdata %>% 
  select(-starts_with("Rh"))

# XRFdata data set is not TIDY! Why? Because we should have 1 column for variable AREA, 
# 1 for CHI and one for STD. And then rows for each element.
# Let's make it TIDY with pivot_longer(). When we do pivot tables, we lose the replicates for XRF data, because
# pivot_ functions don't allow rows with the same values in all columns.
XRFdata_tidy <- XRFdata %>%
  pivot_longer(cols = c(-CoreID, -Sample, -Section),
               names_to = c("Element", "Rest"),
               names_sep = "-") %>% 
  mutate(Rest = recode(Rest,
                       KaChi2 = "Chi",
                       KbChi2 = "Chi",
                       LaChi2 = "Chi",
                       KaArea = "Area",
                       KbArea = "Area",
                       LaArea = "Area",
                       KaAreaStd = "Std",
                       KbAreaStd = "Std",
                       LaAreaStd = "Std")) %>% 
  pivot_wider(names_from = Rest, values_from = value, values_fn = list(value = mean))

# To check that our results give exactly the same using the tidy or untidy version of the data set
# we plot the Aluminium for example.
XRFdata_tidy %>% 
  filter(Element == "S") %>% 
  ggplot() +
  geom_line(aes(x = Sample, y = Area))

XRFdata %>% 
  select(`S-KaArea`, Sample) %>% 
  ggplot() +
  geom_line(aes(x = Sample, y = `S-KaArea`))
# All good, tidying the data set, just removes the replicates every 30 samples. 
# Apart from that, all the rest stands the same.


# Coefficients of variation for each sample
# The coefficient of variation, is just the sd/mean, and this is what we have for each sample (each cm of the core),
# because for each element we have an AREA (=mean) and an STD (=sd).
discarded_by_cv_sample <- XRFdata_tidy %>% 
  mutate(cv_sample = (Std/Area)*100) %>% 
  group_by(Element) %>% 
  summarise(n = n(),
            mean_cv = mean(cv_sample, na.rm = TRUE)) %>% 
  filter(mean_cv<0 | mean_cv>15) %>% 
  select(Element)
# Problematic elements according to this quality check:
# ErmsTancada:  Al (just 11.1), Ar, As, Au, Cr, Cu, Ga, Mn (just 12.2), Mo, Nb, Pb, Se, V
# Llanada:      Ar, As, Au, Cl, Cr, Mo, Nb, P, S, Se, V
# PanissosBIS:  Ar, As, Cr, Cu, Mn, Mo, Nb, P, S, Se, V
# Secanella:    Ar, As, Au, Mo, Nb, P, S, Se, V
# Encanyissada: Ar, As, Cl, Cr, Cu, Ga, Mo, Nb, P, Pb, Se, V

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 3) Checking data quality - Visual inspection of counts vs. sample for each element ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

XRFdata_tidy %>% 
  ggplot() +
  geom_line(aes(x = -Sample, y = Area)) +
  coord_flip() +
  facet_wrap(~Element, scales = "free")

discarded_by_counts <- XRFdata_tidy %>% 
  group_by(Element) %>% 
  summarise(mean_area = mean(Area, na.rm = T),
            max_area = max(Area, na.rm = T),
            median_area = median(Area, na.rm = T)) %>%
  filter(mean_area<1200) %>% 
  select(Element)
# these elements, either have a very noisy plot, or low counts.
# ErmsTancada:  Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Zn
# Llanada:      Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Zn
# PanissosBIS   Ar, As, Au, Cr, Cu, Ga, Mn, Mo, Nb, Ni, P, Pb, Se, V, Zn
# Secanella:    Ar, As, Au, Br, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, S, Se, V, Zn
# Encanyissada: Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Y, Zn


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 4) Final list of elements to be removed due to bad quality ---- 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Lists of problematic elements according to each quality check
discarded_by_cv_replic$Element_short
discarded_by_cv_sample$Element
discarded_by_counts$Element

# Loop to get a list of elements that have been labelled as problematic by ALL quality checkings above
All_elements <- unique(XRFdata_tidy$Element)
discarded_elements <- NULL
for(i in 1:length(All_elements)){
  if(!is.na(match(All_elements[i], discarded_by_cv_replic$Element_short)) | 
     !is.na(match(All_elements[i], discarded_by_cv_sample$Element)) |
     !is.na(match(All_elements[i], discarded_by_counts$Element))){
    discarded_elements <- c(discarded_elements, All_elements[i])
  }
}

# AUTOMATIC List of problematic elements taking all 3 quality checks together 
# (elements present in 1 of the 3 lists will be discarded):
# ErmsTancada:  Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Zn
# Llanada:      Ar, As, Au, Cl, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, S, Se, V, Zn
# PanissosBIS:  Ar, As, Au, Cr, Cu, Ga, Mn, Mo, Nb, Ni, P, Pb, S, Se, V, Zn
# Secanella:    Ar, As, Au, Br, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, S, Se, V, Zn
# Encanyissada: Ar, As, Au, Cl, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Y, Zn

# SUPERVISED list of problematic elements taking all 3 quality checks together:
# ErmsTancada:  Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Zn (same as AUTOMATIC)
# Llanada:      Ar, As, Au, Cl, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Zn (we want to keep S, so use SUPERVISED)
# PanissosBIS:  Ar, As, Au, Cr, Cu, Ga, Mn, Mo, Nb, Ni, P, Pb, S, Se, V, Zn (same as AUTOMATIC)
# Secanella:    Ar, As, Au, Br, Cl, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, S, Se, V, Zn (we want to keep Cl, so use SUPERVISED)
# Encanyissada: Ar, As, Au, Cr, Cu, Ga, Mo, Nb, Ni, P, Pb, Se, V, Y, Zn (same as AUTOMATIC, without Cl, but USE AUTOMATIC)

# We plot all the elements with good quality
if(file_name == "LLANADA1&2_edited.csv"){
   XRFdata_tidy_filtered <- XRFdata_tidy %>% 
   filter(Element %!in% c("Ar", "As", "Au", "Cl", "Cr", "Cu", "Ga", "Mo", "Nb", "Ni", "P", "Pb", "Se", "V", "Zn")) # Llanada
}

if(file_name == "SECANELLA1&2a_edited.csv"){
  XRFdata_tidy_filtered <- XRFdata_tidy %>% 
  filter(Element %!in% c("Ar", "As", "Au", "Br", "Cl", "Cr", "Cu", "Ga", "Mo", "Nb", "Ni", "P", "Pb", "S","Se", "V", "Zn")) # Secanella
}

if(file_name != "SECANELLA1&2a_edited.csv" & file_name != "LLANADA1&2_edited.csv"){
  XRFdata_tidy_filtered <- XRFdata_tidy %>% 
    filter(Element %!in% discarded_elements) # AUTOMATIC list
}



# # # # # # # # # # # # # # # # # # # # #
# 5) Adding colour data to XRF data ---- 
# # # # # # # # # # # # # # # # # # # # #

# Colour data input
colour_data <- read_csv(file = paste("Data/Colour_data/", file_name, sep = ""))

# Summarising colour data for each centimeter, to be able to join colour_data with XRFdata
colour_data_summary <- colour_data %>% 
  mutate(Sample_cm = round(Sample/10)*10) %>% 
  group_by(Sample_cm) %>% 
  summarise(Red = mean(Red),
            Green = mean(Green),
            Blue = mean(Blue),
            l = mean(l), 
            a = mean(a),
            b = mean(b))

# Joining colour_data to XRFdata_tidy_filtered
XRFdata_tidy_filtered <- XRFdata_tidy_filtered %>% 
  left_join(colour_data_summary, by = c("Sample" = "Sample_cm"))


# Clear those data sets that I won't need.
rm(list = c("file_name", "discarded_by_counts", "discarded_by_cv_replic", "discarded_by_cv_sample",
            "XRFdata", "XRFdata_replicates", "All_elements", "i", "index", "XRFdata_tidy",
            "discarded_elements", "colour_data_summary", "colour_data"))


############################
############################
# to add image
# https://stackoverflow.com/questions/9917049/inserting-an-image-to-ggplot2

# p <- XRFdata_tidy %>% 
#   filter(Element == "Al") %>% 
#   ggplot() +
#   geom_line(aes(x = Sample, y = Area))
# 
# logo_file <- system.file("extdata", "logo.png", package = "cowplot")
# 
# ggdraw(p) + 
#   draw_image(logo_file, x = 0.1, y = 0.1, hjust = 0, vjust = 0, width = 0.13, height = 0.2)
# ggsave(filename = "Figures/test2.pdf")
