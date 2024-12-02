# area of applicability
# output: raster stack with AOA of unbiased, biased and all original points + DI of unbiased, biased and all original points

library(CAST)
library(caret)
library(sf)
library(devtools)
library(raster)
library(viridis)
library(ggplot2)
library(tidyverse)
library(terra)
library(geodata)

setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo")


# let's start from the csv we created before (one of them)
# this is the species n. 2, with species prevalence 0.3, sample prevalence 0.9 and 200 random sampled occurrences
sp2_sp_prev0.3_sample_prev0.9_nocc100 <- read.csv("species_2_sp_prevalence_0.3_sample_prev_0.9_n_occ_100.csv",
                                                  row.names = NULL)

# upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# environmental variables 
# bioclimatic variables from worldclim
tmin <- worldclim_country("Ita", "tmin", path=tempdir(), res = 0.5, version = "2.1")
tmax <- worldclim_country("Ita", "tmax", path=tempdir(), res = 0.5, version = "2.1")
prec <- worldclim_country("Ita", "prec", path=tempdir(), res = 0.5, version = "2.1")

# first month only
tmin <- tmin$ITA_wc2.1_30s_tmin_1
tmax <- tmax$ITA_wc2.1_30s_tmax_1
prec <- prec$ITA_wc2.1_30s_prec_1

# stack
r_list <- c(tmin, tmax, prec)
mydata <- raster::stack(r_list)

# crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi_abruzzo)
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)

# csv check 
head(sp2_sp_prev0.3_sample_prev0.9_nocc100)
nrow(sp2_sp_prev0.3_sample_prev0.9_nocc100)

# split the dataset in unbiased and biased
unbiased_random <- sp2_sp_prev0.3_sample_prev0.9_nocc100[sp2_sp_prev0.3_sample_prev0.9_nocc100$UNBIASED == TRUE & sp2_sp_prev0.3_sample_prev0.9_nocc100$BIASED == FALSE, ]
biased <- sp2_sp_prev0.3_sample_prev0.9_nocc100[sp2_sp_prev0.3_sample_prev0.9_nocc100$BIASED == TRUE, ]

head(unbiased_random)
head(biased)

# now we will simulate the same sampling effort, i.e. the two datasets should have more or less the
# same number of points 
# let's say that, if the number of biased points is n, the number of unbiased points will be n + 0.2 * n 
# please note that the unbiased sample must contain the biased one

# n + 0.2 * n 
n_random_to_add <- ceiling(0.2 * nrow(biased))

# random extraction of points BIASED = FALSE (UNBIASED = TRUE), without overlaps with BIASED = TRUE
set.seed(123) 
random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_random_to_add), ]

# combine datasets
unbiased_20 <- rbind(biased, random_points)

# for each of these datasets (unbiased, biased, total occurrences)
# just the environmental variables are needed
head(biased)
head(unbiased_20)
head(sp2_sp_prev0.3_sample_prev0.9_nocc100)

# drop useless variables
drops <- c("distance", "ID", "probability", "UNBIASED", "BIASED", "ID.1")

# subset with environmental data
biased_env <- biased[ , !(names(biased) %in% drops)]
unbiased_20_env <- unbiased_20[ , !(names(unbiased_20) %in% drops)]
sp2_env <- sp2_sp_prev0.3_sample_prev0.9_nocc100[ , !(names(sp2_sp_prev0.3_sample_prev0.9_nocc100) %in% drops)]


# area of applicability
# https://hannameyer.github.io/CAST/articles/cast02-AOA-tutorial.html


# train data: null, biased, all
unbiased_20_env <- unbiased_20_env %>% st_as_sf(., coords = c("X", "Y"), crs = 4326) %>%  as.data.frame()
biased_env <- biased_env %>% st_as_sf(., coords = c("X", "Y"), crs = 4326) %>% as.data.frame()
sp2_env <- sp2_env %>% st_as_sf(., coords = c("X", "Y"), crs = 4326) %>% as.data.frame()

# rasterize environmental data
mydata <- rast(mydata)

# name of the variables to take into account
vars <- names(unbiased_20_env[, !names(unbiased_20_env) %in% c("suitability", "geometry")])

# area of applicability: null, biased, all
aoa_null <- aoa(mydata, train = unbiased_20_env, variables = vars)
aoa_biased <- aoa(mydata, train = biased_env, variables = vars)
aoa_all <- aoa(mydata, train = sp2_env, variables = vars)

# name
base_name <- "sp2_sp_prev0.3_sample_prev0.9_nocc100"

# stack with all the raster
all_rasters <- c(aoa_null$AOA, aoa_null$DI, aoa_biased$AOA, aoa_biased$DI, aoa_all$AOA, aoa_all$DI)

# rename
names(all_rasters) <- c(
  "AOA_null", "DI_null",
  "AOA_biased", "DI_biased",
  "AOA_all", "DI_all"
)

# output as .tif
output_file <- paste0(base_name, "_aoa.tif")
writeRaster(all_rasters, filename = output_file, overwrite = TRUE)
