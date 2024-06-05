library(CAST)
library(caret)
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

remotes::install_github("ropensci/osmextract")
library(osmextract)

setwd("/media/r_projects/phd_rocio/hypervolume")

###### Run manually ######

# Upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# %>% st_union() if there are more regions to merge
plot(aoi_abruzzo)
dev.off()

# Bounding Box 
centro_italia_bb <- st_bbox(aoi_abruzzo)

# From OMSS select type of roads: primary, secondary, tertiary (paths)
ht_primary <- "primary"

# Download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_primary, ]

plot(osm_abruzzo_roads$geometry)

########## Download bioclimatic variables from CHELSA #############################
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


#################  Upload Bioclimatic ################
# String containing the names of raster files
rastlist <- list.files(path ="/media/r_projects/phd_rocio/hypervolume", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# Plot all data
plot(mydata)

# bio1: temperature
plot(mydata[[4]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE, bty="n", box=FALSE)

# Crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi_abruzzo)
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)



########### For plot purposes ####### #########
# Labels 
titles <- c("Mean Annual Temperature", "Annual Precipitation", 
            "Amount of Precipitation in Wettest Month", "Amount of Precipitation in Driest Month")

# Plot all together
par(mfrow=c(2,2), mar=c(2,2,2,0.5))
for (i in 1:nlayers(mydata)) {
  plot(mydata[[i]], main=titles[i], col=magma(500, alpha = 1, begin = 0, end = 1, direction = 1), legend.width=1.5, legend.shrink=0.75)
}



# Original data: will be useful later
mydata_backup <- mydata



########### Random Virtual Species: run every time you want to create a virtual species, from the beginning.
# Suitability map generation
random.sp <- generateRandomSp(raster.stack = mydata,
                              convert.to.PA = FALSE,
                              # How to combine response functions
                              species.type = "multiplicative",
                              # Random approach between PCA and response function
                              approach = "response",
                              # Response function
                              relations = "gaussian",
                              # Realistic species
                              realistic.sp = TRUE,
                              plot = FALSE)


plot(random.sp$suitab.raster, col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability Map", outer=TRUE, line=-1)

random.sp$suitab.raster

# Response functions
# plotResponse(random.sp)

# Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE,
                       species.prevalence = 0.1)



# Plot purposes 
plot(random.sp$suitab.raster)
plot(new.pres$pa.raster, col = c("lightgreen", "orange"))
title("Presence-Absence Map", outer=TRUE, line=-1)

# Occurences
presence.points <- sampleOccurrences(new.pres,
                                     n = 200,
                                     type = "presence only",
                                     sample.prevalence = 0.9,
                                     error.probability = 0,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = FALSE)
dev.off()

# Plot
par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
points(presence.points$sample.points, col = "black", pch = 19, cex=0.3)
title("Occurrences Map", outer=TRUE, line=-1)
dev.off()



######## Preliminary Steps for Niche Analysis #####################
#### Z transform for hypervolume building
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}

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


########### Functions for Hypervolume ######
# Hypervolume
calcola_ipervolume <- function(data) {
  hv_occ <- hypervolume_gaussian(data)
  return(hv_occ@Volume)
}

# Random increment
pippo <- function(x, no){
  # Starts with a random row
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
# Num. simulations each species
num_simulazioni <- 10

# List with the occurrences we want to test
valori_n_occ <- c(50, 80, 100, 150, 200)

# Empty list 
tutte_simulazioni <- list()


# For cycle for simulations
for (sim in 1:num_simulazioni) {
    
     lista_output_occ <- list()
  
      for (i in seq_along(valori_n_occ)) {
        pluto <- pippo(occurrences_values, valori_n_occ[i])
        lista_output_occ[[i]] <- pluto
      }
  
    # Aggiungi la lista delle occorrenze a tutte_simulazioni
    tutte_simulazioni[[sim]] <- lista_output_occ

   }



# All simulations in one df
combined_df <- do.call(rbind, lapply(seq_along(tutte_simulazioni), function(sim) {
   do.call(rbind, lapply(tutte_simulazioni[[sim]], function(df) {
      df$sim <- sim
      df
      }))
    }))


# Mean predictions (LOESS): x sequence 
x_seq <- seq(min(combined_df$n_occ), max(combined_df$n_occ), length.out = 100)

# Mean predictions: LOESS method
loess_predictions <- lapply(unique(combined_df$n_occ), function(n) {
   preds <- sapply(tutte_simulazioni, function(lista) {
      loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
      predict(loess_fit, newdata = data.frame(n_occ = n))
    })
   
   data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))

})


# Mean df
predizioni_media <- do.call(rbind, loess_predictions)

# Plot
ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = predizioni_media, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  labs(title = "Mean (unbiased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()



################# Roadside Bias #########################

# Create raster with distances from roads
roads_vect <- terra::vect(osm_abruzzo_roads$geometry)

raster_roads <- as(mydata_backup[[1]], "SpatRaster")

r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 


# Plot purposes 
d_rast <- d %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(d_rast, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = -1))
title("Distance from Roads (km)", outer=TRUE, line=-1)

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


# Plot purposes
prob_r <- prob_raster %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(prob_r, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = -1))
title("Probability to be sampled", outer=TRUE, line=-1)


# Occurrences as points
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
coord_occ
points(coord_occ)

# Probability of each point to be sampled
probabilities_occ <- terra::extract(prob_raster, coord_occ, ID = TRUE)
probabilities_occ

# Add the probability value to the points
occ_with_prob <- cbind(coord_occ, probabilities_occ)
print(occ_with_prob)

#sample_frac # Points with 100% of probability to be sampled (the one in the roads and within 1 km)
points_biased <- occ_with_prob[occ_with_prob$layer == 1, ]

# If we need a subset: points_biased <- sample(points_biased, 60)


# Plot purposes
par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(prob_r, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = -1))
# title("Probability to be sampled", outer=TRUE, line=-1)
points(coord_occ, cex = 0.5)
points(points_biased, col = "red", cex = 0.5)
# Aggiungi la legenda
legend("topright", legend = c("Unbiased", "Biased"), col = c("black", "red"), pch = 19, cex = 0.8,
       xpd = TRUE, y.intersp = 0.8)

dev.off()

#################### Hypervolume of biased occurrences (road driven: biased sampling) ###############

biased_df <- points_biased %>% as.data.frame()
biased_df <- biased_df[,-c(5:8)]
biased_df

nrow(biased_df)

# Occurrences list
valori_n_occ_biased <- c(20, 34, 57, 68)

# Empty list
tutte_simulazioni_biased <- list()

# Simulations biased df
for (sim in 1:num_simulazioni) {
  
     lista_output_biased <- list()
  
     for (i in seq_along(valori_n_occ_biased)) {
       pluto <- pippo(biased_df, valori_n_occ_biased[i])
       lista_output_biased[[i]] <- pluto
     }
  
  # Aggiungi la lista delle occorrenze biased a tutte_simulazioni_biased
  tutte_simulazioni_biased[[sim]] <- lista_output_biased
  
  }


# Combined df
combined_df_biased <- do.call(rbind, lapply(seq_along(tutte_simulazioni_biased), function(sim) {
    
    do.call(rbind, lapply(tutte_simulazioni_biased[[sim]], function(df) {
    df$sim <- sim
    df
     
    }))
  
  }))


# Mean LOESS
loess_predictions_biased <- lapply(unique(combined_df_biased$n_occ), function(n) {
  
    preds <- sapply(tutte_simulazioni_biased, function(lista) {
      loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
      predict(loess_fit, newdata = data.frame(n_occ = n))
    })
    
  data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  
 })


# Mean in one df
predizioni_media_biased <- do.call(rbind, loess_predictions_biased)

# Plot dei risultati biased
# ggplot(combined_df_biased, aes(x = n_occ, y = iperv)) +
#  geom_point(color = "darkgreen") +
#  geom_smooth(method = "loess", se = TRUE, color = "olivedrab4") +
#  labs(x = "Occurrences", y = "Hypervolume") +
#  theme_minimal()


# Plot
ggplot() +
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = predizioni_media_biased, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
  labs(title = "Mean (biased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()


############## Unbiased & Biased #############
# Same plot
combined_df$total <- "unbiased"
combined_df_biased$total <- "biased"

combined_data <- rbind(combined_df, combined_df_biased)

# Filter NA
predizioni_media <- predizioni_media %>% filter(!is.na(n_occ) & !is.na(iperv_mean))
predizioni_media_biased <- predizioni_media_biased %>% filter(!is.na(n_occ) & !is.na(iperv_mean))

# Function that finds intersection points bewteen curves
find_intersection <- function(df1, df2) {
  intersection_points <- data.frame()
  for (i in 2:nrow(df1)) {
    x1 <- df1$n_occ[i-1]
    x2 <- df1$n_occ[i]
    y1 <- df1$iperv_mean[i-1]
    y2 <- df1$iperv_mean[i]
    
    x3 <- df2$n_occ[i-1]
    x4 <- df2$n_occ[i]
    y3 <- df2$iperv_mean[i-1]
    y4 <- df2$iperv_mean[i]
    
    denominator <- (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if (!is.na(denominator) && denominator != 0) {
      intersect_x <- ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / denominator
      intersect_y <- ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / denominator
      
      if (min(x1, x2) <= intersect_x && intersect_x <= max(x1, x2) &&
          min(x3, x4) <= intersect_x && intersect_x <= max(x3, x4)) {
        intersection_points <- rbind(intersection_points, data.frame(n_occ = intersect_x, iperv = intersect_y))
      }
    }
  }
  return(intersection_points)
}

# Calc. intersection points
intersection_points <- find_intersection(predizioni_media, predizioni_media_biased)

# Combined plot
plot_combined <- ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = predizioni_media, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  geom_line(data = predizioni_media_biased, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
  geom_point(data = intersection_points, aes(x = n_occ, y = iperv), color = "red", size = 2) +
  labs(title = "Hypervolume (Unbiased vs Biased)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()

print(plot_combined)
intersection_points

######### Useful values 
# Convergence values
final_unbiased <- tail(predizioni_media, 1)
final_biased <- tail(predizioni_media_biased, 1)
final_results <- data.frame(
  type = c("unbiased", "biased"),
  n_occ = c(final_unbiased$n_occ, final_biased$n_occ),
  iperv = c(final_unbiased$iperv_mean, final_biased$iperv_mean)
)

# Intersection values
final_results <- rbind(
  final_results,
  data.frame(type = rep("intersection", nrow(intersection_points)), intersection_points)
)

final_results

# CSV
write.csv(final_results, file = "specie_1.csv", row.names = FALSE)



############ Estimating the Area Of Applicability of spatial prediction models ###############
# https://hannameyer.github.io/CAST/articles/cast02-AOA-tutorial.html
################### AOA for Spatially Clustered Data: Null Model vs Biased ###################
##############################################################################################

#################################  Null model ################################################
points(presence.points$sample.points, col = "black", pch = 19, cex=0.3)


## Model training
# A machine learning algorithm will be applied to learn the relationships between predictors and response

## Train data: must be converted in the format required by terra::extract
pa_points <- presence.points$sample.points[,-(3:4)] %>% as.data.frame() %>% st_as_sf(., coords = c("x","y"), crs = 4326)
mydata_aoa <- rast(mydata_backup)

# From raster, extract corresponding values 
trainDat_null <- terra::extract(mydata_aoa, pa_points, na.rm = FALSE)

# From raster, extract suitability values, NA omit, assign spatial reference
trainDat_null$response <- terra::extract(random.sp$suitab.raster, pa_points, na.rm=FALSE, ID=FALSE)
trainDat_null <- data.frame(trainDat_null, pa_points)
trainDat_null <- na.omit(trainDat_null)

## Train model for Spatially Clustered Data
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
trainDat_null
folds_null <- CreateSpacetimeFolds(trainDat_null, spacevar = "geometry", k = 10)


set.seed(15)
model_null <- train(trainDat_null[,names(mydata_aoa)],
                       trainDat_null$response$`VSP suitability`,
                       method = "rf",
                       importance = TRUE,
                       tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                       trControl = trainControl(method ="cv", index = folds$index))

print(model_null)

# Variable Importance of each predictor
plot(varImp(model_null, scale = F), col="black")
plotResponse(random.sp)

## Predict and calculate error 
# The trained model is then used to make predictions for the entire area of interest
prediction_null <- predict(mydata_aoa, model_null, na.rm=T)

# Difference bewteen prediction and reference: true prediction error 
truediff_null <- abs(prediction_null - random.sp$suitab.raster)

# Plot Prediction, Reference and Difference
plot(rast(list(prediction_null, random.sp$suitab.raster, truediff_null)), main = c("Prediction", "Reference", "Difference"), col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

## The AOA calculation takes the model as input to extract the importance of the predictors 
# used as weights in multidimensional distance calculation.
AOA_null <- aoa(mydata_aoa, model_null, LPD = TRUE, verbose = FALSE)

# Features: DI, LPD, AOA
print(AOA_null)

# Plotting the aoa object 
# Shows the distribution of DI values within the training data and the DI of the new data.
plot(AOA_null)

dev.off()

plot(truediff_null, col = viridis(100), main = "True Prediction Error")

# DI: normalized and weighted minimum distance to a nearest training data point 
# divided by the average distance within the training data
plot(AOA_null$DI, col = viridis(100), main = "DI")

# LPD: absolute count of training data points
plot(AOA_null$LPD, col = viridis(100), main = "LPD")

# AOA: derived from the DI by using a threshold.
plot(prediction_null, col=viridis(100), main = "Prediction for AOA")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))

dev.off()


##################################### Biased points ###########################################
###############################################################################################

points(points_biased, col = "red", cex = 0.5)
biased_sp_points <- points_biased %>% st_as_sf(., crs = 4326)
biased_sp_points <- biased_sp_points[,-(1:8)]
biased_sp_points

# From raster, extract corresponding values 
trainDat_biased <- terra::extract(mydata_aoa, biased_sp_points, na.rm=FALSE)

# From raster, extract suitability values 
trainDat_biased$response <- terra::extract(random.sp$suitab.raster, biased_sp_points, na.rm = FALSE, ID=FALSE)
trainDat_biased <- data.frame(trainDat_biased, biased_sp_points)

# Omit NA
trainDat_biased <- na.omit(trainDat_biased)
trainDat_biased
## Train model
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
folds_biased <- CreateSpacetimeFolds(trainDat_biased, spacevar = "geometry", k = 10)
set.seed(15)
model_biased <- train(trainDat_biased[,names(mydata_aoa)],
                    trainDat_biased$response$`VSP suitability`,
                    method = "rf",
                    importance = TRUE,
                    tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                    trControl = trainControl(method ="cv", index = folds_biased$index))

print(model_null)

# Variable Importance of each predictor
plot(varImp(model_biased, scale = F), col="black")
plotResponse(random.sp)

## Predict and calculate error 
# The trained model is then used to make predictions for the entire area of interest
prediction_biased <- predict(mydata_aoa, model_biased, na.rm=T)

# Difference bewteen prediction and reference: true prediction error 
truediff_biased <- abs(prediction_biased - random.sp$suitab.raster)

# Plot Prediction, Reference and Difference
plot(rast(list(prediction_biased, random.sp$suitab.raster, truediff_biased)), main = c("Prediction", "Reference", "Difference"), col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))
dev.off()

plot(rast(list(prediction_random, prediction_biased)), main = c("Random", "Biased"), col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

## The AOA calculation takes the model as input to extract the importance of the predictors 
# used as weights in multidimensional distance calculation.
AOA_biased <- aoa(mydata_aoa, model_biased, LPD = TRUE, verbose = FALSE)

# Features: DI, LPD, AOA
class(AOA_biased)
print(AOA_biased)

# Plotting the aoa object 
# Shows the distribution of DI values within the training data and the DI of the new data.
plot(AOA_biased)

dev.off()

plot(truediff_biased, col = viridis(100), main = "True Prediction Error")

# DI: normalized and weighted minimum distance to a nearest training data point 
# divided by the average distance within the training data
plot(AOA_biased$DI, col = viridis(100), main = "DI")

# LPD: absolute count of training data points
plot(AOA_biased$LPD, col = viridis(100), main = "LPD")


# AOA: derived from the DI by using a threshold.
plot(prediction_biased, col=viridis(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))


dev.off()
############################ Comparison ###############################################

par(mfrow=c(1,2))
plot(prediction_null, col=viridis(100), main = "Prediction for AOA (null)")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

plot(prediction_biased, col=viridis(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))