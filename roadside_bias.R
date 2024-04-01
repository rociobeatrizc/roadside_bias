# Loading required packages
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


## Upload input data
# Set working directory
setwd("C:/chelsa")

# Upload shapefile
aoi <- st_read("parma.shp")
aoi <- aoi$geometry
plot(aoi)


# Download bioclimatic variables from CHELSA
Chelsa.Clim.download(
  # Starting from the workind directory, specify the path
  save.location = "strade_parma",
  # 'bio' contains all the bioclimatic variables
  parameter = "bio",
  # Some variables are chosen from the 19 available
  bio.var = c(1, 7, 13, 12, 14),
  # Version
  version.var = "2.1",
  # Cropping along the area of interest
  clipping = TRUE,
  clip.shapefile = aoi,
  # Insert the coordinates of the area of interest (bounding box)
  clip.extent = c(9.439404, 44.34708, 10.50532, 45.04535),
  # Buffer, if needed
  # buffer = 3,
  # Other commands
  convert.files.to.asc = FALSE,
  stacking.data = TRUE,
  combine.raw.zip = FALSE,
  delete.raw.data = FALSE,
  save.bib.file = TRUE
)

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



## Correlation Matrix
# Subsample 10% of pixels and calculate pairwise correlations
r1 <- mydata$mean.annual.T
cor <- cor(sampleRandom(mydata, size= ncell(r1) * 0.30 ), method = "pearson")

# Plot correlation matrix
df <- corrplot(cor, method = "number", col =  magma(30, alpha = 1, begin = 0, end = 1, direction = 1), type = "lower", tl.pos = 'ld')


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

# Conversion from a dataframe to a raster
raster_pa <- filtered_pa %>% .[,-8] %>% rasterFromXYZ()
plot(raster_pa)

# The raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]

# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r2, r3, r4, r5)
values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))


## Introducing roadside bias 
# The occurrences (data.frame) are transformed into a SpatVector object
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")

# For calculating distances in meters, it is necessary to change the CRS from 4326 to 3857 
# Keeping WGS84 will result in distances in degrees
coord_occ_dist <- terra::project(coord_occ, "EPSG:3857")


# Upload the shapefile containing the road network and convert it into a SpatVector object
setwd("C:/chelsa")
roads <- "strade_parma.shp" %>% st_read() %>% .$geometry %>% terra::vect()

# Here, the CRS also needs to be changed
roads <- terra::project(roads, "EPSG:3857")

# Visualization
plot(roads)
points(coord_occ_dist, col = "red")

# With this function, the distance of each point from all roads is obtained
dist <- distance(coord_occ_dist, roads, unit ="m")

# Distance of each point from the nearest road
min_dist <- apply(dist, 1, min)
summary(min_dist)

# The distances are extracted and associated with the dataset of occurrences
coord_occ$distance <- min_dist

# With a cumulative density function, you can observe the distribution of points relative to roads
# To create the function, distances are extracted from the dataset of occurrences in the form of a dataframe
occ_hist <- coord_occ$distance
occ_hist <- as.data.frame(occ_hist)

# Empirical cumulative distribution function
ecdf_fun <- ecdf(occ_hist$occ_hist)

# Plot
plot(ecdf_fun, verticals = TRUE,
     main = "Cumulative Frequency",
     xlab = "Distance (m)",
     ylab = "Relative Frequency",
     xlim = c(0, 4000), pch = 20, col = "dodgerblue3")


## Sampling probability: since the probability of sampling decreases with distance from the road, we can imagine a decreasing logarithmic curve.
# Sampling probability formula
c <- 1
sampling_prob <- 1-(((log(c*min_dist))/(log(max(c*min_dist)))))

# Sampling probability curve
hist_prob <- bind_cols(sampling_prob = sampling_prob, distance = min_dist) %>%
  ggplot(aes(x = distance, y = sampling_prob)) +
  ylim(0, 1) +
  geom_point(color = "lightgreen") +
  geom_smooth(formula = y ~ log(x), method = "lm", color = "black") +
  labs(x = "distance", y = "sampling probability") +
  theme_bw()

plot(hist_prob)

# Summary distances
summary(min_dist)

# Istogram with number of points per distance
coord_occ %>% as.data.frame() %>% ggplot(., aes(x = distance)) +
  geom_histogram(binwidth = 200, fill = "darkgreen", color = "black") +
  labs(x = "distance", y = "number of points") +
  theme_bw()

# Based on statistics, the intervals' limits are established
lim <- c(0, 200, 500, 800, 1000, 2000, 4000)

# Each point is assigned to its respective class
distance_class <- bind_cols(sampling_prob = sampling_prob, coord_occ %>%
                              as.data.frame(geom="XY")) %>% add_column(class = cut(.$distance, breaks = lim, labels = F))
# Points per class
count(distance_class, class)

# By grouping into classes, for each class, you obtain the initial number of points, the average sampling probability, and the points actually sampled according to this probability
groups <- distance_class %>%
  group_by(class) %>%
  summarise(plotn = n(), meanprob = mean(sampling_prob)) %>%
  add_column(plot_sampl = .$plotn * .$meanprob)

# Union of the two dataframes through coordinates
whole_data <- merge(distance_class, groups, by="class")

## Resampling the original dataset through probability, which acts as a filter, results in a reduced dataset
# Points filtered by probability
by_prob <- whole_data %>%
  group_by(class) %>%
  group_split() %>%
  map(function(z){
    set.seed(1234)
    z %>%
      sample_frac(size = mean(z$sampling_prob))
  }) %>%
  do.call(bind_rows, .)

# In each group, multiplying the initial number of points by the average sampling probability gives a smaller number of points per class
count(by_prob, class)
# In red are the original points, in blue are those actually sampled
plot(roads)
points(whole_data %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
         terra::project(., "EPSG:3857"), col="blue")
points(by_prob %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
         terra::project(., "EPSG:3857"), col="red")


##  Niche analysis
# Let's choose to visualize the ecological niche generated by the points belonging to distance class 3, from 500 meters to 800 meters 
# How much does the niche of this subset differ from the realized niche?
# Filtering points of class 2
points_prob <- by_prob %>% filter(class==2) %>% terra::vect(., geom = c("x","y"),
                                                            crs="epsg:4326")
# Number of points
nrow(points_prob)

# Extract the coordinates of the points belonging to the chosen class
coordinate_prob <- points_prob %>% st_as_sf() %>% sf::st_coordinates()
coordinate_prob

# Create a dataframe
coordinate_prob <- data.frame(x = coordinate_prob[, "X"], y = coordinate_prob[,
                                                                              "Y"])
# Let's merge the two dataframes (original occurrences and points belonging to the chosen class) based on the 'x' and 'y' columns
new_points <- merge(filtered_occ, coordinate_prob, by = c("x", "y"))

# From dataframe to raster
raster_occ <- new_points %>% .[,-(8:9)] %>% rasterFromXYZ()

# With the getValues function, bioclimatic data is obtained in the form of a
# dataframe: the function is applied to both the realized niche and the subset
env_pa <- getValues(raster_pa)
env_occ <- getValues(raster_occ)

# Remove missing values
env_occ <- env_occ[complete.cases(env_occ), ]
env_pa <- env_pa[complete.cases(env_pa), ]

# Produce global environmental background data
globalEnvM <- rbind(env_pa, env_occ)

# PCA on the global data
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)
# Two-dimensional summary of the total environmental variability
global.scores <- pca.clim$li

# The observation data are mapped into that two-dimensional space using the suprow function
pa.scores <- suprow(pca.clim,
         data.frame(filtered_pa)[, colnames(globalEnvM)])$li
occ.scores <- suprow(pca.clim,
         data.frame(filtered_occ)[, colnames(globalEnvM)])$li
pa.scores1 <- suprow(pca.clim, env_pa)$li
occ.scores1<- suprow(pca.clim, env_occ)$li

# Calculate the Occurrence Density Grid
pagrid <- ecospat.grid.clim.dyn(global.scores,
                                pa.scores,
                                pa.scores1)
occgrid <- ecospat.grid.clim.dyn(global.scores,
                                 occ.scores,
                                 occ.scores1)
# Plot niche category
ecospat.plot.niche.dyn(pagrid, occgrid, quant = 0.1, interest = 2, name.axis1 =
                         "PC1", name.axis2 = "PC2")
# Calculate niche overlap
ecospat.niche.overlap(pagrid, occgrid, cor=T)

## Cube generation
# Iteration number
num_iter <- 3

# Empty lists
dataframe_list_before <- list()
dataframe_list_after <- list()

# For loop to generate species and calculate distances for each iteration
for (i in 1:num_iter) {
  set.seed(i)
  random.sp <- generateRandomSp(raster.stack = mydata,
                                convert.to.PA = FALSE,
                                species.type = "multiplicative",
                                approach = "random",
                                relations = "gaussian",
                                realistic.sp = TRUE,
                                plot = FALSE)
  
  new.pres <-convertToPA(random.sp,
                         beta = "random",
                         alpha = -0.05, plot = FALSE,
                         species.prevalence = 0.1)
  presence.points <- sampleOccurrences(new.pres,
                                       n = 100,
                                       type = "presence only",
                                       sample.prevalence = 0.9,
                                       error.probability = 0,
                                       detection.probability = 1,
                                       correct.by.suitability = TRUE,
                                       plot = FALSE)
  
  raster01 <- new.pres$pa.raster %>% raster()
  r1 <- mydata$mean.annual.T
  r2 <- mydata$annual.range.air.T
  r3 <- mydata$annual.precip
  r4 <- mydata$amount.of.prec..wettest.month
  r5 <- mydata$amount.prec..driest.month
  
  stack_pa <- brick(r1, r2, r3, r4, r5, raster01)
  values <- stack_pa %>% rasterToPoints() %>% as.data.frame()
  
  filtered_pa <- values %>% filter(., lyr.1 == 1) %>% as.data.frame()
  raster_pa <- filtered_pa %>% .[,-8] %>% rasterFromXYZ()
  
  raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]
  stack_occ <- brick(r1, r2, r3, r4, r5)
  values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
  
  filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))
  coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
  coord_occ_dist <- terra::project(coord_occ, "EPSG:3857")
  setwd("C:/chelsa")
  roads <- "strade_parma.shp" %>% st_read() %>% .$geometry %>% terra::vect()
  roads <- terra::project(roads, "EPSG:3857")
  
  dist <- distance(coord_occ_dist, roads, unit ="m")
  min_dist <- apply(dist, 1, min)
  coord_occ$distance <- min_dist
  
  c <- 1
  sampling_prob <- 1-(((log(c*min_dist))/(log(max(c*min_dist)))))
  lim <- c(0, 500, 1000, 2000, 5000, 10000, 15000)
  
  distance_class <- bind_cols(sampling_prob = sampling_prob, coord_occ %>%
                                as.data.frame(geom="XY")) %>%
    add_column(class = cut(.$distance, breaks = lim, labels = F))
  
  groups <- distance_class %>%
    group_by(class) %>%
    summarise(plotn = n(), meanprob = mean(sampling_prob)) %>%
    add_column(plot_sampl = .$plotn * .$meanprob)
  whole_data <- merge(distance_class, groups, by="class")
  
  # Time
  set.seed(i)
  whole_data$time <- sample(1980:2010, nrow(whole_data), replace = TRUE)
  
  by_prob <- whole_data %>%
    group_by(class) %>%
    group_split() %>%
    map(function(z){
      set.seed(1234)
      z %>%
        sample_frac(size = mean(z$sampling_prob))
    }) %>%
    do.call(bind_rows, .)
  
  coord_occ <- as.data.frame(whole_data)
  dataframe_before <- whole_data[,-(1:9)]
  dataframe_before <- dataframe_before[,-(4:6)]
  
  # Add species
  by_prob <- as.data.frame(by_prob)
  dataframe_after <- by_prob[,-(1:9)]
  dataframe_after <- dataframe_after[,-(4:6)]
  
  dataframe_before$specie <- paste("specie", i)
  dataframe_after$specie <- paste("specie", i)
  
  dataframe_list_before[[i]] <- dataframe_before
  dataframe_list_after[[i]] <- dataframe_after
}

# Combine all the dataframes obtained into a single dataframe
dataframe_before <- do.call(rbind, dataframe_list_before)
dataframe_after <- do.call(rbind, dataframe_list_after)

# Sort columns
dataframe_before <- dataframe_before[, c("specie", "x", "y", "time", "distance")]
dataframe_after <- dataframe_after[, c("specie", "x", "y", "time", "distance")]

dataframe_after

setwd("C:/chelsa")

# Save the cube!
write.csv(dataframe_before, file = "cube_before.csv", row.names = FALSE)
write.csv(dataframe_after, file = "cube_after.csv", row.names = FALSE)
