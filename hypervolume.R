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

# git add 

## Upload input data
# Set working directory

getwd()
setwd("/media/r_projects/phd_rocio/hackathon_project_8/shapefiles/")

# Upload shapefile
aoi <- st_read("parma.shp")
aoi <- aoi$geometry
plot(aoi)

setwd("/media/r_projects/phd_rocio/hackathon_project_8")

getwd()
# String containing the names of raster files
rastlist <- list.files(path ="bioclim_rasters", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual range air T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# Plot all data
plot(mydata)

# bio1
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))



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
# passi da 200, tutte le occorrenze di deh
presence.points <- sampleOccurrences(new.pres,
                                     n = 500,
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

# Aggiungi l'ID direttamente al dataset filtered_occ
filtered_occ <- filtered_occ %>%
  mutate(ID = row_number())



# Warning messages:
#1: In hypervolume(filtered_pa) : 
#  Consider removing some axes.
# 2: In hypervolume(filtered_pa) :
#  Log number of observations (6.96) is less than or equal to the number of dimensions (7).
# You may not have enough data to accurately estimate a hypervolume with this dimensionality.
# Consider reducing the dimensionality of the analysis.

# Ordina il dataset in base alla colonna 'ID'
filtered_occ <- filtered_occ[order(filtered_occ$ID), ]
filtered_occ


pippo <- function(x, no){
  fx <- x %>% 
    sample_n(size = 1) 
  ipervolumi <- 0
  num_occurrences <- 0
  
  for (i in 1:10000){
    fx <- x %>% 
      sample_n(size = no) %>% 
      bind_rows(fx) %>% 
      distinct()
    
    # Calcola l'ipervolume per il subset di occorrenze
    hv <- calcola_ipervolume(fx)
    
    # Salva l'ipervolume e il numero di occorrenze
    ipervolumi <- c(ipervolumi, hv)
    num_occurrences <- c(num_occurrences, nrow(fx))
    
    if(nrow(fx) == nrow(x)) break   
  }
  bind_cols(iperv = ipervolumi, n_occ = num_occurrences)
  
  
}

filtered_occ %>% distinct() %>% nrow()     

pluto <- pippo(filtered_occ, 200)

pluto
# Inizializza un vettore per salvare gli ipervolumi
ipervolumi <- numeric()
num_occurrences <- numeric()

# Funzione per calcolare l'ipervolume
calcola_ipervolume <- function(data) {
  hv_occ <- hypervolume(data)
  return(hv_occ@Volume)
}
end <- as.numeric(nrow(filtered_occ))
end
# Ciclo for per incrementare il numero di occorrenze da 10 a 100
for (i in seq(10, end, 10)) {
  # Seleziona le prime 'i' occorrenze
  subset_occ <- filtered_occ[1:i, ]
  
  # Calcola l'ipervolume per il subset di occorrenze
  hv <- calcola_ipervolume(subset_occ)
  
  # Salva l'ipervolume e il numero di occorrenze
  ipervolumi <- c(ipervolumi, hv)
  num_occurrences <- c(num_occurrences, i)
}

# Crea il dataframe
df <- data.frame(n_occurrences = num_occurrences, hypervolume = ipervolumi)
pluto
# Stampare il dataframe
print(df)


# Creazione del plot con ggplot
plot <- ggplot(df, aes(x = n_occurrences, y = hypervolume)) +
  geom_line() +  # Linea spezzata
  geom_point() + # Punti
  labs(x = "Numero di Occorrenze", y = "Ipervolume") +  # Etichette degli assi
  ggtitle("Ipervolume all'aumentare delle occorrenze") +   # Titolo del grafico
  theme_minimal() +  # Stile del tema, minimalista 
  theme_bw()

# Visualizza il plot
print(plot)


### Roadside
# The occurrences (data.frame) are transformed into a SpatVector object
coord_occ <- terra::vect(filtered_occ, geom = c("x","y"), crs="epsg:4326")
coord_occ
# For calculating distances in meters, it is necessary to change the CRS
# from 4326 to 3857. Keeping WGS84 will result in distances in degrees
coord_occ_dist <- terra::project(coord_occ, "EPSG:3857")
coord_occ_dist

# Upload the shapefile containing the road network and convert it into a
# SpatVector object
getwd()

roads <- "shapefiles/strade_parma.shp" %>% st_read() %>% .$geometry %>% terra::vect()

# Here, the CRS also needs to be changed
roads <- terra::project(roads, "EPSG:3857")

# Visualization
plot(roads)
points(coord_occ_dist, col = "red", pch = 16)

# With this function, the distance of each point from all roads is obtained
dist <- distance(coord_occ_dist, roads, unit ="m")

# Distance of each point from the nearest road
min_dist <- apply(dist, 1, min)

summary(min_dist)

# The distances are extracted and associated with the dataset of occurrences
coord_occ$distance <- min_dist

# Based on statistics, the intervals' limits are established
lim <- c(50000, 120000, 140000, 150000, 170000)

## Sampling probability: since the probability of sampling decreases with distance from the road, we can imagine a decreasing logarithmic curve.
# Sampling probability formula
c <- 0.0001
sampling_prob <- 1-(((log(c*min_dist))/(log(max(c*min_dist)))))

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
whole_data

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

points_closer_to_roads <- by_prob %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
  terra::project(., "EPSG:3857")

points_closer_to_roads

# Converti l'SpatVector in un dataframe
df_points_closer <- as.data.frame(points_closer_to_roads)

# Stampare il dataframe
print(df_points_closer)
df_points_closer <- df_points_closer[,-1]
df_points_closer <- df_points_closer[,-(7:11)]

df_points_closer
# Inizializza un vettore per salvare gli ipervolumi
ipervolumi_roads <- numeric()
num_occurrences_roads <- numeric()

# Funzione per calcolare l'ipervolume
calcola_ipervolume <- function(data) {
  hv_occ <- hypervolume(data)
  return(hv_occ@Volume)
}


end_roads <- as.numeric(nrow(df_points_closer))

filtered_occ

df_points_closer

ipervolumi_roads <- c()
num_occurrences_roads <- c()
# Ciclo for per incrementare il numero di occorrenze da 10 a 100
# Partire da 10
for (j in seq(10, end_roads, 15)) {
  # Seleziona le prime 'i' occorrenze
  subset_occ_roads <- df_points_closer[1:j, ]
  
  # Calcola l'ipervolume per il subset di occorrenze
  hv_roads <- calcola_ipervolume(subset_occ_roads)
  
  # Salva l'ipervolume e il numero di occorrenze
  ipervolumi_roads <- c(ipervolumi_roads, hv_roads)
  num_occurrences_roads <- c(num_occurrences_roads, j)
}

# Crea il dataframe
df_closer <- data.frame(n_occurrences_roads = num_occurrences_roads, hypervolume_roads = ipervolumi_roads)

# Stampare il dataframe
print(df_closer)

# Creazione del plot con ggplot
plot_1 <- ggplot(df_closer, aes(x = n_occurrences_roads, y = hypervolume_roads)) +
  geom_line() +  # Linea spezzata
  geom_point() + # Punti
  labs(x = "Numero di Occorrenze", y = "Ipervolume") +  # Etichette degli assi
  ggtitle("Ipervolume all'aumentare delle occorrenze") +   # Titolo del grafico
  theme_minimal() +  # Stile del tema, minimalista 
  theme_bw()

# Visualizza il plot
print(plot_1)

# Creazione del plot con ggplot sovrapposto
plot_combined <- ggplot() +
  geom_line(data = df, aes(x = n_occurrences, y = hypervolume), color = "blue", linetype = "solid") +
  geom_point(data = df, aes(x = n_occurrences, y = hypervolume), color = "blue") +
  geom_line(data = df_closer, aes(x = n_occurrences_roads, y = hypervolume_roads), color = "red", linetype = "dashed") +
  geom_point(data = df_closer, aes(x = n_occurrences_roads, y = hypervolume_roads), color = "red") +
  labs(x = "Numero di Occorrenze", y = "Ipervolume") +  
  ggtitle("Ipervolume all'aumentare delle occorrenze") +  
  theme_minimal() + 
  theme_bw()

# Visualizzazione del plot combinato
print(plot_combined)
