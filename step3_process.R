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

# rasterize environmental data
mydata <- rast(mydata)


setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo/virtualspecies_output")


# function: from csv to raster
process_csv_to_raster <- function(input_folder) {
  
  # all csv files with the same structure
  csv_files <- list.files(input_folder, pattern = "^species_\\d+_sp_prevalence_\\d+\\.\\d+_sample_prev_\\d+\\.\\d+_n_occ_\\d+\\.csv$", full.names = TRUE)
  
  # loop
  for (csv_file in csv_files) {
    # name
    base_name <- tools::file_path_sans_ext(basename(csv_file))
    
    # read csv
    data <- read.csv(csv_file, row.names = NULL)
    
    # unbiased and biased
    unbiased_random <- data[data$UNBIASED == TRUE & data$BIASED == FALSE, ]
    biased <- data[data$BIASED == TRUE, ]
    
    # same sampling effort
    n_random_to_add <- ceiling(0.2 * nrow(biased))
    set.seed(123)
    random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_random_to_add), ]
    
    # combine
    unbiased_20 <- rbind(biased, random_points)
    
    # remove useless columns
    drops <- c("distance", "ID", "probability", "UNBIASED", "BIASED", "ID.1")
    biased_env <- biased[ , !(names(biased) %in% drops)]
    unbiased_20_env <- unbiased_20[ , !(names(unbiased_20) %in% drops)]
    sp2_env <- data[ , !(names(data) %in% drops)]
    
    # into sf 
    unbiased_20_env <- unbiased_20_env %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% as.data.frame()
    biased_env <- biased_env %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% as.data.frame()
    sp2_env <- sp2_env %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% as.data.frame()
    
    # drop variables
    vars <- names(unbiased_20_env[, !names(unbiased_20_env) %in% c("suitability", "geometry")])
    
    # aoa
    aoa_null <- aoa(mydata, train = unbiased_20_env, variables = vars)
    aoa_biased <- aoa(mydata, train = biased_env, variables = vars)
    aoa_all <- aoa(mydata, train = sp2_env, variables = vars)
    
    # raster stack
    all_rasters <- c(aoa_null$AOA, aoa_null$DI, aoa_biased$AOA, aoa_biased$DI, aoa_all$AOA, aoa_all$DI)
    names(all_rasters) <- c(
      "AOA_null", "DI_null",
      "AOA_biased", "DI_biased",
      "AOA_all", "DI_all"
    )
    
    # save raster files
    output_file <- file.path(input_folder, paste0(base_name, "_aoa.tif"))
    writeRaster(all_rasters, filename = output_file, overwrite = TRUE)
    
    cat("File salvato:", output_file, "\n")
  }
}


# let's go
process_csv_to_raster(".")
