# Loading required packages
install.packages("devtools")
install.packages("parallel")
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
library(grid)
library(gridExtra)
library(ecospat)
library(ade4)
library(hypervolume)

install.packages("hypervolume")
devtools::install_github("bblonder/hypervolume")

# git add 

## Upload input data
# Set working directory
setwd("C:/chelsa")

# Upload shapefile
aoi <- st_read("parma.shp")
aoi <- aoi$geometry
plot(aoi)

# String containing the names of raster files
rastlist <- list.files(path ="strade_parma/bio/ChelsaV2.1Climatologies/clipped_2024-02-28_09-42-52", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual range air T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# Plot all data
plot(mydata)

# bio1
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

detectCores()

## Random Virtual Species
# Suitability map generation
random.sp <- generateRandomSp(raster.stack = mydata,
                              convert.to.PA = FALSE,
                              # How to combine response functions
                              species.type = "multiplicative",
                              # Random approach between PCA and response function
                              approach = "random",
                              # Response function
                              relations = "gaussian",
                              # Realistic species
                              realistic.sp = TRUE,
                              plot = FALSE)
plot(random.sp$suitab.raster, col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

# Response functions
plotResponse(random.sp)

# Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE,
                       species.prevalence = 0.1)
plot(new.pres$pa.raster)

# Occurences
presence.points <- sampleOccurrences(new.pres,
                                     n = 100,
                                     type = "presence only",
                                     sample.prevalence = 0.9,
                                     error.probability = 0,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = FALSE)
# Plot
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(presence.points$sample.points, col = "black", pch = 19)



### Z transform
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}

plot(mydata)

## Preliminary Steps for Niche Analysis
# Raster with presence/absence points (step 2): used to create the realized niche
# It must be a RasterLayer object
raster01 <- new.pres$pa.raster %>% raster()

# The bioclimatic layers are extracted one by one from the stack
r1 <- mydata$mean.annual.T
r2 <- mydata$annual.range.air.T
r3 <- mydata$annual.precip
r4 <- mydata$amount.of.prec..wettest.month
r5 <- mydata$amount.prec..driest.month

# A stack is created containing the bioclimatic variables and the raster of presence/absence (realized niche)
stack_pa <- brick(r1, r2, r3, r4, r5, raster01)

# The values are extracted from the stack, converted into a dataframe, and only the presence pixels (1) are retained
values <- stack_pa %>% rasterToPoints() %>% as.data.frame()
filtered_pa <- values %>% filter(., lyr.1 == 1) %>% as.data.frame()
filtered_pa


# Conversion from a dataframe to a raster
raster_pa <- filtered_pa %>% .[,-8] %>% rasterFromXYZ()
plot(raster_pa)

filtered_pa <- filtered_pa[,-8]
# The raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]

# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r2, r3, r4, r5)
values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
filtered_occ <- filtered_occ[,-(8:9)]
filtered_occ

# Nicchia realizzata
hv_pa = hypervolume(filtered_pa)

# Occorrenze: null model
hv_occ = hypervolume(filtered_occ)




# Warning messages:
# 1: In hypervolume(filtered_pa) : 
# Consider removing some axes.
# 2: In hypervolume(filtered_pa) :
# Log number of observations (6.96) is less than or equal to the number of dimensions (7).
# You may not have enough data to accurately estimate a hypervolume with this dimensionality.
# Consider reducing the dimensionality of the analysis.

pa_seq = hypervolume_resample("pa_seq", hv_pa, "bootstrap seq", n = 3, seq = seq(10, 100, 10), cores = 1)
occ_seq = hypervolume_resample("occ_seq", hv_occ, "bootstrap seq", n = 3, seq = seq(10, 100, 10), cores = 1)

?hypervolume_resample
# Funnel Plots
occ_plot = hypervolume_funnel(occ_seq) + 
  geom_point(aes(y = upperq)) + 
  geom_point(aes(y = lowerq)) + 
  geom_point(aes(y = sample_mean), col = "blue") + 
  theme_bw() + 
  labs(title = "a)", subtitle = NULL) + 
  ylab("Volume") + 
  ylim(0, 0.3) + 
  xlim(0, 100)
pa_plot = hypervolume_funnel(pa_seq) + 
  geom_point(aes(y = upperq)) + 
  geom_point(aes(y = lowerq)) + 
  geom_point(aes(y = sample_mean), col = "blue") + 
  theme_bw() + 
  labs(title = "b)", subtitle = NULL) + 
  ylab("Volume") + 
  ylim(0, 0.3) + 
  xlim(0, 100)
grid.arrange(pa_plot, occ_plot, nrow = 1)
