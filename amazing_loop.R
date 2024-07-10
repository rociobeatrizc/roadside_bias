########### The amazing loop! #####################
########### No plots ##############################

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

# Upload shapefile
aoi_abruzzo <- st_read("abruzzo.shp") %>% .$geometry 

# Bounding Box 
centro_italia_bb <- st_bbox(aoi_abruzzo)

# From OMSS select type of roads: primary, secondary, tertiary (paths)
ht_primary <- "primary"

# Download roads from OSM 
osm_abruzzo <- oe_get("Abruzzo", stringsAsFactors = FALSE, quiet = TRUE)
osm_abruzzo_roads <- osm_abruzzo[osm_abruzzo$highway %in% ht_primary, ]


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

# Crop and mask by region borders
aoi_sp <- sf::as_Spatial(aoi_abruzzo)
mydata <- mydata %>% crop(., aoi_sp) %>% mask(., aoi_sp)


# Original data: will be useful later
mydata_backup <- mydata


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


######################## Functions for Hypervolume ####################
# Hypervolume
calcola_ipervolume <- function(data) {
  hv_occ <- hypervolume_gaussian(data)
  return(hv_occ@Volume)
}

# Random increment
pippo <- function(x, no
                  # , epsilon = 0.01
                  ) {
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
    if(nrow(fx) == nrow(x)) {
      break
    }
  }
  
  result <- bind_cols(iperv = ipervolumi, n_occ = num_occurrences)
  return(list(result))
         #, interrotta_per_convergenza = interrotta_per_convergenza))
}

############################# Roadside Bias ################################
############################################################################
# Create raster with distances from roads
roads_vect <- terra::vect(osm_abruzzo_roads$geometry)

raster_roads <- as(mydata_backup[[1]], "SpatRaster")

r <- terra::rasterize(roads_vect, raster_roads)
d <- distance(r, unit = "km") 


############################################################################

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

##############################################################################################
######################## Amazing Loop ########################################################
# Define the function

simulate_species <- function(nocc, nsim, n_species) {
  
  final_results_list <- list()
  
  for (species in 1:n_species) {
    cat("\nSpecie", species, "\n")
    
    # Suitability map generation
    random.sp <- generateRandomSp(raster.stack = mydata,
                                  convert.to.PA = FALSE,
                                  species.type = "multiplicative",
                                  approach = "response",
                                  relations = "gaussian",
                                  realistic.sp = TRUE,
                                  plot = FALSE)
    
    # Presence/Absence
    new.pres <- convertToPA(random.sp,
                            beta = "random",
                            alpha = -0.05, plot = FALSE,
                            species.prevalence = 0.1)
    
    # Occurrences
    presence.points <- sampleOccurrences(new.pres,
                                         n = nocc,
                                         type = "presence only",
                                         sample.prevalence = 0.9,
                                         error.probability = 0,
                                         detection.probability = 1,
                                         correct.by.suitability = TRUE,
                                         plot = FALSE)
    
    # Filter occurrences
    raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]
    stack_occ <- brick(r1, r3, r4, r5)
    values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
    filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
    occurrences_values <- filtered_occ[,-c(1:2, 7:8)]
    
    # Occurrences as points
    coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
    
    # Probability of each point to be sampled
    probabilities_occ <- terra::extract(prob_raster, coord_occ, ID = TRUE)
    occ_with_prob <- cbind(coord_occ, probabilities_occ)
    
    # Biased points
    points_biased <- occ_with_prob[occ_with_prob$layer == 1, ]
    
    # Hypervolume of occurrences (random sampled: null model)
    num_simulazioni <- nsim
    nrow(points_biased)
    nrow(occurrences_values)
    stop <- ceiling(nrow(points_biased) + 0.2 * (nrow(points_biased)))
    occurrences_values <- occurrences_values[sample(nrow(occurrences_values), stop), ]
    nrow(occurrences_values)
    
    # List with occurrences to test
    valori_n_occ <- c(seq(from = 40, to = stop, by = 10), stop)
    tutte_simulazioni <- list()
    convergenza_info <- vector("logical", num_simulazioni)
    
    # Simulations loop
    for (sim in 1:num_simulazioni) {
      lista_output_occ <- list()
      
      for (i in seq_along(valori_n_occ)) {
        pluto <- pippo(occurrences_values, valori_n_occ[i])
        lista_output_occ[[i]] <- pluto[[1]]
      }
      convergenza_info[sim] <- interrotta_per_convergenza_locale
      tutte_simulazioni[[sim]] <- lista_output_occ
    }

    
    combined_df <- do.call(rbind, lapply(seq_along(tutte_simulazioni), function(sim) {
      do.call(rbind, lapply(tutte_simulazioni[[sim]], function(df) {
        df$sim <- sim
        df
      }))
    }))
    
    x_seq <- seq(min(combined_df$n_occ), max(combined_df$n_occ), length.out = 100)
    loess_predictions <- lapply(unique(combined_df$n_occ), function(n) {
      preds <- sapply(tutte_simulazioni, function(lista) {
        loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
        predict(loess_fit, newdata = data.frame(n_occ = n))
      })
      data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
    })
    predizioni_media <- do.call(rbind, loess_predictions)
    
    # Hypervolume of biased occurrences (road driven: biased sampling)
    biased_df <- points_biased %>% as.data.frame()
    biased_df <- biased_df[,-c(5:8)]
    stop_biased <- nrow(biased_df)
    valori_n_occ_biased <- c(seq(from = 30, to = stop_biased, by = 20), stop_biased)
    tutte_simulazioni_biased <- list()
    convergenza_info_biased <- vector("logical", num_simulazioni)
    
    for (sim in 1:num_simulazioni) {
      lista_output_biased <- list()

      
      for (i in seq_along(valori_n_occ_biased)) {
        pluto <- pippo(biased_df, valori_n_occ_biased[i])
        lista_output_biased[[i]] <- pluto[[1]]
        
      }

      tutte_simulazioni_biased[[sim]] <- lista_output_biased
    }
    

    
    combined_df_biased <- do.call(rbind, lapply(seq_along(tutte_simulazioni_biased), function(sim) {
      do.call(rbind, lapply(tutte_simulazioni_biased[[sim]], function(df) {
        df$sim <- sim
        df
      }))
    }))
    
    loess_predictions_biased <- lapply(unique(combined_df_biased$n_occ), function(n) {
      preds <- sapply(tutte_simulazioni_biased, function(lista) {
        loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
        predict(loess_fit, newdata = data.frame(n_occ = n))
      })
      data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
    })
    predizioni_media_biased <- do.call(rbind, loess_predictions_biased)
    
    # Combine unbiased and biased data
    combined_df$total <- "unbiased"
    combined_df_biased$total <- "biased"
    combined_data <- rbind(combined_df, combined_df_biased)
    
    predizioni_media <- predizioni_media %>% filter(!is.na(n_occ) & !is.na(iperv_mean))
    predizioni_media_biased <- predizioni_media_biased %>% filter(!is.na(n_occ) & !is.na(iperv_mean))
    
    final_unbiased <- tail(predizioni_media, 1)
    final_biased <- tail(predizioni_media_biased, 1)
    final_results <- data.frame(
      species = species,
      type = c("unbiased", "biased"),
      n_occ = c(final_unbiased$n_occ, final_biased$n_occ),
      iperv = c(final_unbiased$iperv_mean, final_biased$iperv_mean),
      points_biased = nrow(points_biased),
      points_nb_20_percent = nrow(occurrences_values),
      pixels_non_biased_raster = NA,  # placeholder
      pixels_biased_raster = NA,  # placeholder
      biased_null = NA,  # placeholder
      null_biased = NA  # placeholder
    )
    
    
    
    
    # AOA estimation
    # Unbiased
    #################################  Null model ################################################
    
    
    ## Model training
    # A machine learning algorithm will be applied to learn the relationships between predictors and response
    
    ## Train data: must be converted in the format required by terra::extract
    pa_points <- presence.points$sample.points[,-(3:4)] %>% as.data.frame() %>% st_as_sf(., coords = c("x","y"), crs = 4326)
    mydata_aoa <- rast(mydata_backup)
    
    # Subset of the original 200 points
    pa_points <- pa_points[rownames(occurrences_values), ]
    
    
    # From raster, extract corresponding values 
    trainDat_null <- terra::extract(mydata_aoa, pa_points, na.rm = FALSE)
    
    # From raster, extract suitability values, NA omit, assign spatial reference
    trainDat_null$response <- terra::extract(random.sp$suitab.raster, pa_points, na.rm=FALSE, ID=FALSE)
    trainDat_null <- data.frame(trainDat_null, pa_points)
    trainDat_null <- na.omit(trainDat_null)
    
    ## Train model for Spatially Clustered Data
    folds_null <- CreateSpacetimeFolds(trainDat_null, spacevar = "geometry", k = 4)
    
    
    set.seed(15)
    model_null <- train(trainDat_null[,names(mydata_aoa)],
                        trainDat_null$response$`VSP suitability`,
                        method = "rf",
                        importance = TRUE,
                        tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                        trControl = trainControl(method ="cv", index = folds_null$index))
    
    
    ## Predict and calculate error 
    # The trained model is then used to make predictions for the entire area of interest
    prediction_null <- predict(mydata_aoa, model_null, na.rm=T)
    
    
    # Biased
    biased_sp_points <- points_biased %>% st_as_sf(., crs = 4326)
    biased_sp_points <- biased_sp_points[,-(1:8)]
    trainDat_biased <- terra::extract(mydata_aoa, biased_sp_points, na.rm=FALSE)
    trainDat_biased$response <- terra::extract(random.sp$suitab.raster, biased_sp_points, na.rm = FALSE, ID=FALSE)
    trainDat_biased <- data.frame(trainDat_biased, biased_sp_points)
    trainDat_biased <- na.omit(trainDat_biased)
    
    folds_biased <- CreateSpacetimeFolds(trainDat_biased, spacevar = "geometry", k = 10)
    set.seed(15)
    model_biased <- train(trainDat_biased[,names(mydata_aoa)],
                          trainDat_biased$response$`VSP suitability`,
                          method = "rf",
                          importance = TRUE,
                          tuneGrid = expand.grid(mtry = c(2:length(names(mydata_aoa)))),
                          trControl = trainControl(method = "repeatedcv",
                                                   repeats = 10,
                                                   index = folds_biased$index,
                                                   savePredictions = "all"),
                          metric = "RMSE")
    
    
    
    # The trained model is then used to make predictions for the entire area of interest
    prediction_biased <- predict(mydata_aoa, model_biased, na.rm=T)
    AOA_null <- aoa(mydata_aoa, model_null, LPD = TRUE, verbose = FALSE)
    AOA_biased <- aoa(mydata_aoa, model_biased, LPD = TRUE, verbose = FALSE)
    
    masked_raster_null <- mask(prediction_null, AOA_null$AOA, maskvalues=0, updatevalue=NA)
    masked_raster_biased <- mask(prediction_biased, AOA_biased$AOA, maskvalues=0, updatevalue=NA)
    
    count_non_na_null <- sum(!is.na(masked_raster_null[]))
    count_non_na_biased <- sum(!is.na(masked_raster_biased[]))
    
    diff_null_only <- ifel(!is.na(masked_raster_null) & is.na(masked_raster_biased), 1, NA)
    diff_biased_only <- ifel(is.na(masked_raster_null) & !is.na(masked_raster_biased), -1, NA)
    diff_raster <- merge(diff_null_only, diff_biased_only)
    pixel_values <- values(diff_raster)
    num_red_pixels <- sum(pixel_values == -1, na.rm = TRUE)
    num_blue_pixels <- sum(pixel_values == 1, na.rm = TRUE)
    
    final_results$pixels_non_biased_raster <- count_non_na_null
    final_results$pixels_biased_raster <- count_non_na_biased
    final_results$biased_null <- num_red_pixels
    final_results$null_biased <- num_blue_pixels
    
    cat("Specie", species, "\n")
    cat("points biased", nrow(points_biased), "\n")
    cat("points nb 20%", nrow(occurrences_values), "\n")
    print(final_results)
    cat("pixels in non-biased raster:", count_non_na_null, "\n")
    cat("pixels in biased raster:", count_non_na_biased, "\n")
    cat("biased-null", num_red_pixels, "\n")
    cat("null-biased", num_blue_pixels, "\n")
    
    final_results_list[[species]] <- final_results
  }
  
  # Combine all results and save to CSV
  combined_final_results <- do.call(rbind, final_results_list)
  write.csv(combined_final_results, "final_results_300occ_10sim.csv", row.names = FALSE)
  
  return(combined_final_results)
}

# Example call to the function
results <- simulate_species(nocc = 300, nsim = 10, n_species = 3)
