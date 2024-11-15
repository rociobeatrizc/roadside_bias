# first step: random species generation
# each species has a number of occurrences
# the output is a csv of points in which we can find: coordinates of the points, distance from the roads, probability to be sampled,
# bioclimatic variables, type of point according to the lazy sampler simulation (unbiased/biased)
# parameters to be setted before: number of starting occurrences, species prevalence, sample prevalence
# those infos will be inside the name of the csv

library(sf)
library(raster)
library(virtualspecies)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(terra)
library(geodata)
library(osmdata)
library(osmextract)


# set wd
setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo")

# upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# plot region
plot(aoi_abruzzo)

# bounding box 
abruzzo_bb <- st_bbox(aoi_abruzzo)

# from OSM select type of roads: primary, secondary, tertiary (paths)
ht_secondary <- "secondary"

# download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_secondary, ]

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

# original data: will be useful later (AOA)
mydata_backup <- mydata


# distances
roads_vect <- terra::vect(osm_abruzzo_roads$geometry)
raster_roads <- as(mydata_backup[[1]], "SpatRaster")

r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 

# raster with distances
d_raster <- d %>% raster()

# probability to be sampled
distances <- d_raster %>% as.data.frame()
c <- 1
sampling_prob <- 1 - (((log(c * distances)) / (log(max(c * distances))))) %>% as.data.frame()

sampling_prob[sampling_prob == Inf] <- 1
sampling_prob[sampling_prob > 1] <- 1

# probability as a raster
prob_raster <- classify(d, cbind(values(d), sampling_prob))

# parameters to be saved
# species prevalence
sp_prevalence <- 0.1

# sample prevalence
sample_prev <- 0.9

# n. occurrences
n_occ <- 300

# number of species
n_species <- 4

# cicle for virtual species generation 
for (species_num in 1:n_species) {
  
  # virtual species: suitability
  random.sp <- generateRandomSp(raster.stack = mydata,
                                convert.to.PA = FALSE,
                                species.type = "multiplicative",
                                approach = "response",
                                relations = "gaussian",
                                realistic.sp = TRUE,
                                plot = FALSE)
  
  # presence/absence
  new.pres <- convertToPA(random.sp,
                          beta = "random",
                          alpha = -0.05, plot = FALSE,
                          species.prevalence = sp_prevalence)
  
  # random occurrences
  presence.points <- sampleOccurrences(new.pres,
                                       n = n_occ,
                                       type = "presence only",
                                       sample.prevalence = sample_prev,
                                       error.probability = 0,
                                       detection.probability = 1,
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  
  # Z-transform bioclimatic variables
  for (i in 1:nlayers(mydata)) {
    mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd')
  }
  
  # occurrences as points
  initial_occ <- presence.points$sample.points %>%
    as.data.frame() %>%
    .[.$Real == 1 & .$Observed == 1, ] %>%
    st_as_sf(., coords = c("x","y"))
  
  # add bioclimatic values
  drops <- c("Real", "Observed")
  initial_occ_bio_var <- terra::extract(mydata, initial_occ) %>%
    cbind(., initial_occ) %>%
    st_set_crs(4326) %>%
    .[, !(names(.) %in% drops)]
  
  # add distances
  initial_occ_bio_var <- d_raster %>%
    terra::extract(., initial_occ_bio_var) %>%
    cbind(initial_occ_bio_var, .)
  
  # rename
  names(initial_occ_bio_var)[names(initial_occ_bio_var) == "."] <- "distance"
  
  # add probability to be sampled for each point
  initial_occ_bio_var <- prob_raster %>%
    terra::extract(., initial_occ_bio_var) %>%
    cbind(initial_occ_bio_var, .)
  
  # rename
  names(initial_occ_bio_var)[names(initial_occ_bio_var) == "layer"] <- "probability"
  
  # unbiased and biased columns
  # every point is initially unbiased (random sampling)
  # if the probability to be sampled is == 1, the point could fall into the 'biased' category
  initial_occ_bio_var <- initial_occ_bio_var %>%
    mutate(
      UNBIASED = TRUE, 
      BIASED = probability == 1 
    )
  
  
  # dataframe
  df <- as.data.frame(st_coordinates(initial_occ_bio_var))
  
  # add variables
  df <- cbind(df, initial_occ_bio_var[, !(names(initial_occ_bio_var) %in% c("geometry"))])
  
  # save as csv
  # the name will include: the number of the species, species prevalence, sample prevalence and number of initial occurrences
  file_name <- paste0("species_", species_num, "_sp_prevalence_", sp_prevalence, 
                      "_sample_prev_", sample_prev, "_n_occ_", n_occ, ".csv")
  write.csv(df, file_name, row.names = FALSE)

}

