# Packages
library(CAST)
library(caret)
library(sf)
library(ClimDatDownloadR)
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

# Set working directory
setwd("/media/r_projects/phd_rocio/hypervolume")

# Upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# %>% st_union() if there are more regions to merge

# Plot region
plot(aoi_abruzzo)

# Bounding Box 
abruzzo_bb <- st_bbox(aoi_abruzzo)

# From OSM select type of roads: primary, secondary, tertiary (paths)
ht_primary <- "primary"

# Download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_primary, ]
plot(osm_abruzzo_roads$geometry)

## Download bioclimatic variables from CHELSA
Chelsa.Clim.download(
  # Starting from the workind directory, specify the path
  #  save.location = "my/path",
  
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
  clip.extent = abruzzo_bb,
  
  # Buffer, if needed
  # buffer = 3,
  
  # Other commands
  convert.files.to.asc = FALSE,
  stacking.data = TRUE,
  combine.raw.zip = FALSE,
  delete.raw.data = FALSE,
  save.bib.file = TRUE
)


##  Upload bioclimatic variables
# String containing the names of raster files
rastlist <- list.files(path ="/media/r_projects/phd_rocio/hypervolume", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# Plot all data
plot(mydata)

# bio1: temperature
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE, bty="n", box=FALSE)

# Crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi_abruzzo)
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)

## Plot: roads on raster
raster_df <- as.data.frame(rasterToPoints(mydata_backup[[1]]), xy = TRUE)
value_column <- names(raster_df)[3]

ggplot() +
  # Add raster
  geom_raster(data = raster_df, aes_string(x = "x", y = "y", fill = value_column)) +
  scale_fill_viridis_c() + 
  # Roads
  geom_sf(data = osm_abruzzo_roads$geometry, color = "black", size = 0.5) +
  theme_minimal() +
  labs(title = "Roads",
       fill = "Values") +
  coord_sf()

## Plot: bioclimatic variables
# Labels 
titles <- c("Mean Annual Temperature", "Annual Precipitation", 
            "Amount of Precipitation in Wettest Month", "Amount of Precipitation in Driest Month")
# Plot all together
par(mfrow=c(2,2), mar=c(2,2,2,0.5))
for (i in 1:nlayers(mydata_backup)) {
  plot(mydata_backup[[i]], main=titles[i], col=magma(500, alpha = 1, begin = 0, end = 1, direction = 1), legend.width=1.5, legend.shrink=0.75, axes=FALSE, box=FALSE)
}
dev.off()

## Original data: will be useful later
mydata_backup <- mydata

## Random Virtual Species: run every time you want to create a virtual species, from the beginning.

## Suitability map generation
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


# Suitability plot
plot(random.sp$suitab.raster, col = plasma(500, alpha = 1, begin = 0, end = 1, direction = 1))
title("Suitability Map", outer=TRUE, line=-1)
dev.off()

# Response functions
# plotResponse(random.sp)

## Presence/Absence: requires defining the parameters alpha, beta, and species prevalence
new.pres <-convertToPA(random.sp,
                       beta = "random",
                       alpha = -0.05, plot = FALSE,
                       species.prevalence = 0.1)


# Presence/Absence plot
plot(random.sp$suitab.raster)
plot(new.pres$pa.raster, col = c("yellowgreen", "deeppink"), box = FALSE, axes = FALSE)
title("Presence-Absence Map", outer=TRUE, line=-1)

## Occurences
presence.points <- sampleOccurrences(new.pres,
                                     n = 200,
                                     type = "presence only",
                                     sample.prevalence = 0.9,
                                     error.probability = 0,
                                     detection.probability = 1,
                                     correct.by.suitability = TRUE,
                                     plot = FALSE)

# Plot: occurrences
par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(random.sp$suitab.raster, col = plasma(500, alpha = 1, begin = 0, end = 1, direction = 1), axes = FALSE, box = FALSE)
points(presence.points$sample.points, col = "black", pch = 19, cex=0.3)
title("Occurrences", outer=TRUE, line=-1)
dev.off()


## Preliminary Steps for Niche Analysis

# Z transform for hypervolume building
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}

# Raster with presence/absence points (step 2): used to create the realized niche
# It must be a RasterLayer object
raster01 <- new.pres$pa.raster %>% raster()
raster01

######### RIVEDERE!!! ############
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

##### RIVEDERE!!!: con funzione extract
# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r3, r4, r5)
values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
values_occ

filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))

occurrences_values <- filtered_occ[,-c(1:2, 7:8)]


## Functions for Hypervolume

# Hypervolume: just the hypervolume value from hypervolume_gaussian function
hyp_calc <- function(data) {
  hv_occ <- hypervolume_gaussian(data)
  return(hv_occ@Volume)
}

# Function to build the accumulation curve with random increment in occurrences
acc_curve <- function(x, no) {
  # Starts with a random row
  fx <- x %>% 
    sample_n(size = 1) 
  
  ipervolumi <- 0
  num_occurrences <- 0
  
  for (i in 1:1000) {
    
   # To the initial value (a row)
   # Random values are selected
   # They are bound to fx
   # Unique values are kept
    fx <- x %>% 
      sample_n(size = no) %>% 
      bind_rows(fx) %>% 
      distinct()
    
    # Hypervolume per subset
    hv <- hyp_calc(fx)
    
    # Save hypervolume & number of occurrences
    ipervolumi <- c(ipervolumi, hv)
    num_occurrences <- c(num_occurrences, nrow(fx))

   # Condition
   # Stop when the subset has the same number of occurrences as the original set
    if(nrow(fx) == nrow(x)) {
      break
    }
  }
  
  result <- bind_cols(iperv = ipervolumi, n_occ = num_occurrences)
  return(list(result))
}



## Roadside bias
# Create raster with distances from roads
roads_vect <- terra::vect(osm_abruzzo_roads$geometry)

# Turn into SpatRaster object
raster_roads <- as(mydata_backup[[1]], "SpatRaster")

# Rasterize distances
r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 

## Plot: distance from roads
d_rast <- d %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
raster_df_dist <- as.data.frame(rasterToPoints(d_rast), xy = TRUE)
value_column <- names(raster_df_dist)[3]

ggplot() +
  # Add raster
  geom_raster(data = raster_df_dist, aes_string(x = "x", y = "y", fill = value_column)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1) +  
  # Add roads
  geom_sf(data = osm_abruzzo_roads$geometry, color = "black", size = 0.5) +
  theme_bw() +
  theme_minimal() +
  labs(title = "Distance from Roads",
       fill = "Distance (km)") +
  coord_sf() +
  theme(
     axis.title.x = element_blank(),
     axis.title.y = element_blank()
   )


## Extract distances
d_raster <- d %>% raster()
distances <- d_raster %>%  as.data.frame()

## Sampling probability: simulation of the lazy sampler
c <- 1
sampling_prob <- 1-(((log(c*distances))/(log(max(c*distances)))))
sampling_prob <- as.data.frame(sampling_prob)

# Some values are: Inf. Replace those values with 1
sampling_prob[sampling_prob == Inf] <- 1
sampling_prob[sampling_prob > 1] <- 1

# New raster with probability to be sampled instead of distances
prob_raster <- classify(d, cbind(values(d), sampling_prob))

## Plot: sampling probability
prob_r <- prob_raster %>% raster() %>% crop(., aoi_sp) %>% mask(., aoi_sp)
raster_df_prob <- as.data.frame(rasterToPoints(prob_r), xy = TRUE)
value_column <- names(raster_df_prob)[3]

ggplot() +
  # Add raster
  geom_raster(data = raster_df_prob, aes_string(x = "x", y = "y", fill = value_column)) +
  scale_fill_viridis_c(alpha = 1, begin = 0, end = 1) +  
  # Add roads
  geom_sf(data = osm_abruzzo_roads$geometry, color = "black", size = 0.5) +
  theme_bw() +
  theme_minimal() +
  labs(title = "Sampling Probability",
       fill = "Probability") +
  coord_sf() +
  theme(
     axis.title.x = element_blank(),
     axis.title.y = element_blank()
   )


## Occurrences as points
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")

# Probability of each point to be sampled
probabilities_occ <- terra::extract(prob_raster, coord_occ, ID = TRUE)

# Add the probability value to the points
occ_with_prob <- cbind(coord_occ, probabilities_occ)

# Points with 100% of probability to be sampled (those on the roads and/or within 1 km)
points_biased <- occ_with_prob[occ_with_prob$layer == 1, ]


## Hypervolume of occurrences (random sampled: null model) 
# Num. simulations each species
num_sim <- 10

# Set stop point according to the number of biased occurrences: same sampling effort
nrow(points_biased)
nrow(occurrences_values)
stop <-  ceiling(nrow(points_biased) + 0.2 * (nrow(points_biased)))

# Random subsample of occurrences from null model 
occurrences_values <- occurrences_values[sample(nrow(occurrences_values), stop), ]

############################# Plot purposes: map with unbiased-biased points #####################
# Index
indices <- rownames(occurrences_values)
indices <- as.numeric(indices)
filtered_coord_occ <- coord_occ[indices, ]

par(mfrow=c(1,1), mar=c(2,2,2,0.5)) 
plot(prob_r, col = viridis(500, alpha = 1, begin = 0, end = 1, direction = 1))
# title("Probability to be sampled", outer=TRUE, line=-1)
points(filtered_coord_occ, cex = 0.6)
points(points_biased, col = "red", cex = 0.6)
# Aggiungi la legenda
legend("topright", legend = c("Unbiased", "Biased"), col = c("black", "red"), pch = 19, cex = 0.8,
       xpd = TRUE, y.intersp = 0.8)

dev.off()

coord_occ
occurrences_values





# List with the occurrences we want to test
valori_n_occ <- c(seq(from = 20, to = stop, by = 20), stop)

# Empty list 
tutte_simulazioni <- list()

convergenza_info <- vector("logical", num_simulazioni)

# For cycle for simulations
for (sim in 1:num_sim) {
  
  lista_output_occ <- list()
  interrotta_per_convergenza_locale <- FALSE
  
  for (i in seq_along(valori_n_occ)) {
    pluto <- pippo(occurrences_values, valori_n_occ[i])
    lista_output_occ[[i]] <- pluto[[1]]
    
    # Aggiorna la variabile di stato locale in base alla presenza di convergenza in questa simulazione
    if (pluto$interrotta_per_convergenza) {
      interrotta_per_convergenza_locale <- TRUE
    }
  }
  
  # Aggiorna la lista delle simulazioni interrotte per convergenza
  convergenza_info[sim] <- interrotta_per_convergenza_locale
  
  # Aggiungi la lista delle occorrenze a tutte_simulazioni
  tutte_simulazioni[[sim]] <- lista_output_occ
  
}

# Output delle informazioni sulla convergenza
print(convergenza_info)

tutte_simulazioni
# All simulations in one df
combined_df <- do.call(rbind, lapply(seq_along(tutte_simulazioni), function(sim) {
   do.call(rbind, lapply(tutte_simulazioni[[sim]], function(df) {
      df$sim <- sim
      df
      }))
    }))

combined_df
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

############################ Plot: unbiased hypervolume ###############################
ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = predizioni_media, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  labs(title = "Mean (unbiased dataset)",
       x = "Occurrences",
       y = "Hypervolume") +
  theme_minimal()

#################### Hypervolume of biased occurrences (road driven: biased sampling) ###############
#####################################################################################################
biased_df <- points_biased %>% as.data.frame()
biased_df <- biased_df[,-c(5:8)]
biased_df

# Stop
stop_biased <- nrow(biased_df)
valori_n_occ_biased <- c(seq(from = 20, to = stop_biased, by = 20), stop_biased)

# Empty list
tutte_simulazioni_biased <- list()
convergenza_info_biased <- vector("logical", num_simulazioni)

# For cycle for simulations
for (sim in 1:num_simulazioni) {
  
  lista_output_biased <- list()
  interrotta_per_convergenza_locale_biased <- FALSE
  
  for (i in seq_along(valori_n_occ_biased)) {
    pluto <- pippo(biased_df, valori_n_occ_biased[i])
    lista_output_biased[[i]] <- pluto[[1]]
    
    # Aggiorna la variabile di stato locale in base alla presenza di convergenza in questa simulazione
    if (pluto$interrotta_per_convergenza) {
      interrotta_per_convergenza_locale_biased <- TRUE
    }
  }
  
  # Aggiorna la lista delle simulazioni interrotte per convergenza
  convergenza_info_biased[sim] <- interrotta_per_convergenza_locale_biased
  
  # Aggiungi la lista delle occorrenze a tutte_simulazioni
  tutte_simulazioni_biased[[sim]] <- lista_output_biased
  
}

# Output delle informazioni sulla convergenza
print(convergenza_info_biased)

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


########################## Plot: biased hypervolume ###############################
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
intersection_points

############################### Combined plot: biased unbiased ###################################
plot_combined <- ggplot() +
  geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_smooth(data = combined_df_biased, aes(x = n_occ, y = iperv, group = sim), 
              method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
  geom_line(data = predizioni_media, aes(x = n_occ, y = iperv_mean), 
            color = "sienna1", size = 1.2) +
  geom_line(data = predizioni_media_biased, aes(x = n_occ, y = iperv_mean), 
            color = "darkgreen", size = 1.2) +
#  geom_point(data = intersection_points, aes(x = n_occ, y = iperv), color = "red", size = 2) +
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


## Model training
# A machine learning algorithm will be applied to learn the relationships between predictors and response

## Train data: must be converted in the format required by terra::extract
pa_points <- presence.points$sample.points[,-(3:4)] %>% as.data.frame() %>% st_as_sf(., coords = c("x","y"), crs = 4326)
mydata_aoa <- rast(mydata_backup)

pa_points
# Subset of the original 200 points
rownames(occurrences_values)
pa_points <- pa_points[rownames(occurrences_values), ]
pa_points


# From raster, extract corresponding values 
trainDat_null <- terra::extract(mydata_aoa, pa_points, na.rm = FALSE)

# From raster, extract suitability values, NA omit, assign spatial reference
trainDat_null$response <- terra::extract(random.sp$suitab.raster, pa_points, na.rm=FALSE, ID=FALSE)
trainDat_null <- data.frame(trainDat_null, pa_points)
trainDat_null <- na.omit(trainDat_null)

## Train model for Spatially Clustered Data
# train from CARET package: data train, data output, method (Random Forest) and Cross Validation 
trainDat_null
folds_null <- CreateSpacetimeFolds(trainDat_null, spacevar = "geometry", k = 4)


set.seed(15)
model_null <- train(trainDat_null[,names(mydata_aoa)],
                       trainDat_null$response$`VSP suitability`,
                       method = "rf",
                       importance = TRUE,
                       tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                       trControl = trainControl(method ="cv", index = folds_null$index))

print(model_null)

# Variable Importance of each predictor
plot(varImp(model_null, scale = F), col="black", main = "Importance of each predictor", axes =FALSE)
plotResponse(random.sp)
dev.off()


## Predict and calculate error 
# The trained model is then used to make predictions for the entire area of interest
prediction_null <- predict(mydata_aoa, model_null, na.rm=T)

# Difference bewteen prediction and reference: true prediction error 
truediff_null <- abs(prediction_null - random.sp$suitab.raster)


# Plot Prediction, Reference and Difference
par(mfrow = c(1, 2)) 
plot(prediction_null, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

dev.off()


plot(truediff_null, main = "Difference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

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
plot(prediction_null, col=inferno(100), main = "Prediction for Area of Applicability")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))

dev.off()


##################################### Biased points ###########################################
###############################################################################################
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
par(mfrow = c(1, 2)) 
plot(prediction_biased, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

dev.off()


# Plot Prediction, Reference and Difference
par(mfrow = c(1, 2)) 
plot(prediction_biased, main = "Prediction with RF", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(random.sp$suitab.raster, main = "Reference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))
plot(truediff_biased, main = "Difference", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))

par(mfrow = c(1, 2)) 
plot(prediction_null, main = "RF Null Model", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1), legend = FALSE)
plot(prediction_biased, main = "RF Biased data", col = inferno(500, alpha = 1, begin = 0, end = 1, direction = 1))
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
plot(prediction_biased, col=inferno(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topleft", box.col = "black", bty = "o", title = "AOA"))


dev.off()
############################ Comparison ###############################################
## Set same scale
par(mfrow=c(1,2))
plot(prediction_null, col=viridis(100), main = "Prediction for AOA (Null)")
plot(AOA_null$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

plot(prediction_biased, col=viridis(100), main = "Prediction for AOA (Biased)")
plot(AOA_biased$AOA, col = c("grey","transparent"), add = T, plg = list(x = "topright", box.col = "black", bty = "o", title = "AOA"))

#######################################################################################
model_null$results
model_biased$results

## Difference? Show in the map (calc. pixel)
plot(prediction_biased)
plot(AOA_biased$AOA)

# Uniased Masked 
masked_raster_null <- mask(prediction_null, AOA_null$AOA, maskvalues=0, updatevalue=NA)

# Biased Masked
masked_raster_biased <- mask(prediction_biased, AOA_biased$AOA, maskvalues=0, updatevalue=NA)

par(mfrow=c(1,3))
plot(masked_raster_null)
plot(masked_raster_biased)

dev.off()

# Pixels in Null Model only
diff_null_only <- ifel(!is.na(masked_raster_null) & is.na(masked_raster_biased), 1, NA)

# Pixels in Biased Model only
diff_biased_only <- ifel(is.na(masked_raster_null) & !is.na(masked_raster_biased), -1, NA)

# Merge
diff_raster <- merge(diff_null_only, diff_biased_only)

# Palette
col_palette <- c("deeppink", "darkgreen")

# Plot
par(mfrow = c(1, 3), mar = c(5, 4, 4, 4) + 0.1)
plot(masked_raster_null, main = "Null", col=viridis(100), legend =FALSE)
plot(masked_raster_biased, main = "Biased", col=viridis(100), legend = FALSE)
plot(diff_raster, col = col_palette, main = "Difference", legend = FALSE)
par(mar = c(5, 4, 4, 4) + 0.1, xpd = TRUE)
legend("topleft", legend = c("Bias - Null", "Null - Bias"), fill = col_palette, cex = 0.8, bty = "n")
par(mfrow = c(1, 1))
dev.off()


#### Spatial difference ####
pixel_values <- values(diff_raster)

# Num. red and blue pixels
num_red_pixels <- sum(pixel_values == -1, na.rm = TRUE)
num_blue_pixels <- sum(pixel_values == 1, na.rm = TRUE)

# Area of each kind of pixel: useless
# area_red_km2 <- num_red_pixels * 1  # Area in km^2
# area_blue_km2 <- num_blue_pixels * 1  # Area in km^2

# Print
cat("N. red pixels:", num_red_pixels)
cat("Area red pixels (km^2):", area_red_km2)
cat("N. blu pixels", num_blue_pixels)
cat("Area blu pixels (km^2):", area_blue_km2)


################## Stabilire un epsilon: condizione di convergenza sulla base di cosa? #########################################################################
################## hypervolume_overlap_statistics(): indici di Sorensen o Jaccard (Similarity) fra due ipervolumi consecutivi nell'iterazione ##################
