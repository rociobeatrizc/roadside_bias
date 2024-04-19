# Loading required packages
install.packages("devtools")
install.packages("gdata")
library(sf)
library(ClimDatDownloadR)
# Alternatively, a direct link with GitHub can be created
if(!require(devtools)) install.packages("devtools")
library(devtools)
devtools::install_github("HelgeJentsch/ClimDatDownloadR")
library(raster)
library(viridis)
library(corrplot)
library(virtualspecies)
library(ggplot2)
library(tidyverse)
library(terra)
library(ecospat)
library(ade4)
library(hypervolume)
library(gdata)
install.packages("hypervolume")


devtools::install_github("bblonder/hypervolume")
## Hypervolume: git add hypervolume.R
# E mo??

## Not run: 
# 3000 data point hypervolume 
data(quercus) 
hv_quercus = hypervolume(quercus[,c(2,3)]) 
hypervolume::resam
# the seq argument is equivalent to a length 30 vector {10, 139, ... , 3649, 3779} # 6hr sequential
quercus_bootstrap_seq <- resample('quercus_bootstrap_seq', hv_quercus, method = 'bootstrap seq', points_per_resample = "sample_size", seq = floor(seq(10, 3779, length.out = 30)), cores = 20) 

# Compatible with ggplot syntax when used with as_table = FALSE 
hypervolume_funnel(quercus_bootstrap_seq, title = 'Resampled volumes of Quercus', func = get_volume) + 
  geom_line(aes(y = get_volume(hv_quercus))) + ylab("Volume") 