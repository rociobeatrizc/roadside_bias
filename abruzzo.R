install.packages("osmextract")
install.packages("ggmap")
install.èackages("tidylog")

library(sf)
library(ClimDatDownloadR)
# Alternatively, a   direct link with GitHub can be created
if(!require(devtools)) install.packages("devtools")
library(devtools)
# devtools::install_github("HelgeJentsch/ClimDatDownloadR")
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
library(ggmap)
library(osmdata)
library(osmextract)
# install.packages("remotes")
remotes::install_github("ropensci/osmextract")
library(osmextract)

setwd("/media/r_projects/phd_rocio/hackathon_project_8/centro_ita")

# Upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry %>% st_union()
plot(aoi_abruzzo)
dev.off()
centro_italia_bb <- st_bbox(aoi_abruzzo)
ht_primary <- "tertiary"
# Download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_primary, ]

plot(osm_abruzzo_roads$geometry)

# Download bioclimatic variables from CHELSA
Chelsa.Clim.download(
  # Starting from the workind directory, specify the path
#  save.location = "strade_parma",
  # 'bio' contains all the bioclimatic variables
  parameter = "bio",
  # Some variables are chosen from the 19 available
  bio.var = c(1, 7, 13, 14),
  # Version
  version.var = "2.1",
  # Cropping along the area of interest
  clipping = TRUE,
  clip.shapefile = aoi_abruzzo,
  # Insert the coordinates of the area of interest (bounding box)
  clip.extent = centro_italia_bb,
  # Buffer, if needed
  # buffer = 3,
  # Other commands
  convert.files.to.asc = FALSE,
  stacking.data = TRUE,
  combine.raw.zip = FALSE,
  delete.raw.data = FALSE,
  save.bib.file = TRUE
)


##  Upload
# String containing the names of raster files
rastlist <- list.files(path ="/media/r_projects/phd_rocio/hackathon_project_8/centro_ita/bio/ChelsaV2.1Climatologies/clipped_2024-05-08_15-30-29/", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# Plot all data
plot(mydata)

# bio1: temperature
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

# Crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi_abruzzo)
aoi_sp
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)

plot(mydata)

# Original data: will be useful later
mydata_backup <- mydata

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
# plotResponse(random.sp)

# Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE,
                       species.prevalence = 0.1)
plot(new.pres$pa.raster)

# Occurences
presence.points <- sampleOccurrences(new.pres,
                                     n = 200,
                                     type = "presence only",
                                     sample.prevalence = 0.9,
                                     error.probability = 0,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = FALSE)

# Plot
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(presence.points$sample.points, col = "black", pch = 19)

dev.off()

### Z transform: preliminary step for hypervolume building
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}


plot(mydata_backup[[4]])
plot(mydata[[4]])

# Base pipe |>
# %<>% assignment pipe
# %$% exposition pipe
# Tee pipe %T>%
# tidylog package: comment the code

## Preliminary Steps for Niche Analysis
# Raster with presence/absence points (step 2): used to create the realized niche
# It must be a RasterLayer object
raster01 <- new.pres$pa.raster %>% raster()
raster01

# The bioclimatic layers are extracted one by one from the stack
r1 <- mydata$mean.annual.T
# r2 <- mydata$annual.range.air.T
r3 <- mydata$annual.precip
r4 <- mydata$amount.of.prec..wettest.month
r5 <- mydata$amount.prec..driest.month

# A stack is created containing the bioclimatic variables and the raster of presence/absence (realized niche)
stack_pa <- brick(r1, r3, r4, r5, raster01)

# The raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]
raster_occurences

# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r3, r4, r5)
values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
values_occ

filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))

occurrences_values <- filtered_occ[,-c(1:2, 7:8)]
occurrences_values
############################# Function for hypervolume ############
# Hypervolume
calcola_ipervolume <- function(data) {
  hv_occ <- hypervolume_gaussian(data)
  return(hv_occ@Volume)
}

# Random increment
pippo <- function(x, no){
  # Una riga a caso dal dataset occorrenze: inizializza
  fx <- x %>% 
    sample_n(size = 1) 
  
  ipervolumi <- 0
  num_occurrences <- 0
  
  for (i in 1:1000) {
    
    # Al valore di inizio (una riga)
    # Scelgo a caso no valori
    # Li lego ad fx
    # Mantengo valori univoci
    fx <- x %>% 
      sample_n(size = no) %>% 
      bind_rows(fx) %>% 
      distinct()
    
    # Ipervolume per subset
    hv <- calcola_ipervolume(fx)
    
    # Salva l'ipervolume e il numero di occorrenze
    ipervolumi <- c(ipervolumi, hv)
    num_occurrences <- c(num_occurrences, nrow(fx))
    
    # Condizione
    # Interrompi quando il subset ha lo stesso numero di occorrenze del set originale
    if(nrow(fx) == nrow(x)) break   
  }
  bind_cols(iperv = ipervolumi, n_occ = num_occurrences)
  
}


#################### Hypervolume of occurrences (random sampled: null model) ###############
nrow(occurrences_values)

# List of the occurrences we want to test
valori_n_occ <- c(20, 50, 80, 100, 150, 200)

lista_output_occ <- list()

# Loop for the function 
for (i in seq_along(valori_n_occ)) {
  pluto <- pippo(occurrences_values, valori_n_occ[i])
  lista_output_occ[[i]] <- pluto
}

##################### Output ##################
lista_output_occ
length(lista_output_occ)

# All together
df_total_occ <- do.call(rbind, lista_output_occ)

# Plot
ggplot(df_total_occ, aes(x = n_occ, y = iperv)) +
  geom_point() + # Punti
  geom_smooth(method = "loess", se = TRUE) +  # Linea di interpolazione
  labs(x = "n_occ", y = "iperv", title = "graph") +
  theme_minimal()

###### Roadside Bias #######
plot(mydata_backup[[1]])
plot(osm_abruzzo_roads$geometry)

dev.off()

############# Create raster with distances from roads #######################
roads_vect <- terra::vect(osm_abruzzo_roads$geometry)

raster_roads <- as(mydata_backup[[1]], "SpatRaster")

r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 

# %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
plot(d)

# Extract distances
d_raster <- d %>% raster()
distances <- d_raster %>%  as.data.frame()
distances

# Sampling probability
c <- 1
sampling_prob <- 1-(((log(c*distances))/(log(max(c*distances)))))
sampling_prob <- as.data.frame(sampling_prob)
sampling_prob

# Some values are: Inf. Replace those values with 1
sampling_prob[sampling_prob == Inf] <- 1
sampling_prob[sampling_prob > 1] <- 1

# New raster with probability to be sampled instead of distances
prob_raster <- classify(d, cbind(values(d), sampling_prob))
plot(prob_raster)

########## Occurrences as points
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
coord_occ
points(coord_occ)

############## Probability of each point to be sampled
probabilities_occ <- terra::extract(prob_raster, coord_occ, ID = TRUE)
probabilities_occ
# Add the probability value to the points
occ_with_prob <- cbind(coord_occ, probabilities_occ)
print(occ_with_prob)

#sample_frac# Points with 100% of probability to be sampled (the one in the roads and within 1 km)
points_biased <- occ_with_prob[occ_with_prob$layer == 1, ]
# points_biased <- sample(points_biased, 60)
points_biased

plot(mydata[[1]])
points(coord_occ)
points(points_biased, col = "red")

dev.off()

points_biased

# Dataframe to build hypervolume (biased)
biased_df <- points_biased %>% as.data.frame()
biased_df <- biased_df[,-c(5:8)]
biased_df

nrow(biased_df)

# List of the occurrences we want to test
valori_n_occ_biased <- c(10, 20, 30, 50)

lista_output_biased <- list()

# Loop for the function 
for (i in seq_along(valori_n_occ_biased)) {
  pluto <- pippo(biased_df, valori_n_occ_biased[i])
  lista_output_biased[[i]] <- pluto
}

lista_output_occ

length(lista_output_biased)

lista_output_biased

# Plot
ggplot(lista_output_biased[[1]], aes(x = n_occ, y = iperv)) +
  geom_point() + # Punti
  geom_smooth(method = "loess", se = TRUE) +  # Linea di interpolazione
  labs(x = "n_occ", y = "iperv", title = "graph") +
  theme_minimal()

# Unico dataframe
dataframe_outputs_occ <- do.call(rbind, lista_output_occ)
dataframe_outputs_biased <- do.call(rbind, lista_output_biased)

# Combinati
plot_combined <- ggplot() +
  geom_point(data = dataframe_outputs_occ, aes(x = n_occ, y = iperv, color = "unbiased")) +
  geom_point(data = dataframe_outputs_biased, aes(x = n_occ, y = iperv, color = "biased")) +
  geom_smooth(data = dataframe_outputs_occ, aes(x = n_occ, y = iperv, color = "blue"), method = "loess", se = TRUE) +
  geom_smooth(data = dataframe_outputs_biased, aes(x = n_occ, y = iperv, color = "red"), method = "loess", se = TRUE) +
  labs(x = "n_occ", y = "iperv", color = "Dataset") +  
  ggtitle("ipervolume") +  
  theme_minimal() + 
  theme_bw() +
  scale_color_manual(values = c("unbiased" = "blue", "biased" = "red"),
                     labels = c("unbiased" = "Unbiased", "biased" = "Biased"))

print(plot_combined)
