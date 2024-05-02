# Loading required packages
install.packages("devtools")
install.packages("gdata")
library(sf)
library(ClimDatDownloadR)
# Alternatively, a direct link with GitHub can be created
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
install.packages("hypervolume")
devtools::install_github("bblonder/hypervolume")

getwd()
setwd("/media/r_projects/phd_rocio/hackathon_project_8/shapefiles")

# Upload shapefile
aoi <- st_read("parma.shp")
aoi <- aoi$geometry
plot(aoi)


setwd("/media/r_projects/phd_rocio/hackathon_project_8")

# String containing the names of raster files
rastlist <- list.files(path ="/media/r_projects/phd_rocio/hackathon_project_8/raster_parma", pattern = "CHELSA", full.names = TRUE)

# Using the list of names, all the files are imported into a single raster package
mydata <- stack(rastlist)

# Change data names
names(mydata) <- c("mean annual T", "annual precip", "amount of prec. wettest month", "amount prec. driest month")

# bio1
plot(mydata[[1]], col = magma(500, alpha = 1, begin = 0, end = 1, direction = 1))

## Correlation Matrix
# Subsample 10% of pixels and calculate pairwise correlations
r1 <- mydata$mean.annual.T
cor <- cor(sampleRandom(mydata, size= ncell(r1) * 0.90 ), method = "pearson")

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
# passi da 200, tutte le occorrenze di deh
presence.points <- sampleOccurrences(new.pres,
                                     n =1000,
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


### Z transform
for (i in 1:nlayers(mydata)){
  mydata[[i]] <- (mydata[[i]] - cellStats(mydata[[i]], 'mean')) / cellStats(mydata[[i]], 'sd') 
}

# This step is about re-calibrate the raster values: crs is still the same
mydata

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

# The values are extracted from the stack, converted into a dataframe, and only the presence pixels (1) are retained
values <- stack_pa %>% rasterToPoints() %>% as.data.frame()
filtered_pa <- values %>% filter(., lyr.1 == 1) %>% as.data.frame()
filtered_pa <- filtered_pa[,-7]

# Nicchia realizzata con rispettivi valori climatici (step 2)
realized_niche_values <- filtered_pa[,-(1:2)]
realized_niche_values

# 1050 cells represent the realized niche values
# Lets calculate the hypervolume of the realized niche first
nrow(realized_niche_values)
?distinct


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

############################# Realized Niche Test ########### ####################
lista_output <- list()
valori <- c(100, 200, 300, 400, 600, 800, 1000)

# Loop for the function 
for (i in seq_along(valori_n_occ)) {
  pluto <- pippo(realized_niche_values, valori_n_occ[i])
  lista_output[[i]] <- pluto
}

lista_output[[1]]

# Plot con ggplot2
ggplot(lista_output[[3]], aes(x = n_occ, y = iperv)) +
  geom_point() + # Punti
  geom_smooth(method = "loess", se = TRUE) +  # Linea di interpolazione
  labs(x = "n_occ", y = "iperv", title = "graph") +
  theme_minimal()

?geom_smooth

### Tutti insieme
dataframe_outputs <- do.call(rbind, lista_output)

############################# Occurrences ###################
# Occurrences
# Conversion from a dataframe to a raster (to plot)
raster_pa <- filtered_pa %>% .[,-7] %>% rasterFromXYZ()
plot(raster_pa)


# The raster of occurrences is transformed into a dataset, from which the rows satisfying both conditions Real = 1 and Observed = 1 are preserved
raster_occurences <- presence.points$sample.points %>% as.data.frame() %>% .[.$Real == 1 & .$Observed == 1, ]
raster_occurences


# The environmental variables are associated with the occurrences using their coordinates
stack_occ <- brick(r1, r3, r4, r5)
values_occ <- stack_occ %>% rasterToPoints() %>% as.data.frame()
values_occ

filtered_occ <- merge(values_occ, raster_occurences, by = c("x", "y"))

occurrences_values <- filtered_occ[,-c(1:2, 7:8)]
# occurrences_values <- occurrences_values[,-(1:2)]

## Climatic values of the occurrences (step 3)
nrow(occurrences_values)
occurrences_values %>%  sample_n(size=1)

# List of the occurrences we want to test
valori_n_occ <- c(200, 500, 700, 990)

lista_output_occ <- list()

# Loop for the function 
for (i in seq_along(valori_n_occ)) {
  pluto <- pippo(occurrences_values, valori_n_occ[i])
  lista_output_occ[[i]] <- pluto
}

lista_output_occ
length(lista_output_occ)

### Tutti insieme
# df_total_occ <- do.call(rbind, lista_output_occ)

# Plot
ggplot(lista_output_occ[[3]], aes(x = n_occ, y = iperv)) +
  geom_point() + # Punti
  geom_smooth(method = "loess", se = TRUE) +  # Linea di interpolazione
  labs(x = "n_occ", y = "iperv", title = "graph") +
  theme_minimal()

##################### Roadside ##################################
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
# roads <- (roads - cellStats(roads, 'mean')) / cellStats(roads, 'sd')
# Here, the CRS also needs to be changed
roads <- terra::project(roads, "EPSG:3857")

# Visualization
plot(roads)
same.crs(roads, coord_occ_dist)
points(coord_occ_dist, col = "red", pch = 16)

# With this function, the distance of each point from all roads is obtained
dist <- distance(coord_occ_dist, roads, unit ="m")

# Distance of each point from the nearest road
min_dist <- apply(dist, 1, min)

summary(min_dist)

# The distances are extracted and associated with the dataset of occurrences
coord_occ$distance <- min_dist

# Based on statistics, the intervals' limits are established
lim <- c(0, 100, 1000, 2000, 5000, 10000, 20000, 27000)

## Sampling probability: since the probability of sampling decreases with distance from the road, we can imagine a decreasing logarithmic curve.
# Sampling probability formula
c <- 1
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

# In blue are the original points, in red are those actually sampled
plot(roads)
points(whole_data %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
         terra::project(., "EPSG:3857"), col="blue")
points(by_prob %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
         terra::project(., "EPSG:3857"), col="red")
dev.off()
occurrences_values
biased_df <- by_prob %>% terra::vect(., geom = c("x","y"), crs="epsg:4326") %>%
  terra::project(., "EPSG:3857") %>% as.data.frame()

biased_df <- biased_df[,-c(1:2, 7:12)]
biased_df

nrow(biased_df)

# List of the occurrences we want to test
valori_n_occ_biased <- c(50, 100, 150, 250)

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

dev.off()
################### salva separatamente ######################
# Inizializza una lista vuota per salvare i risultati
# lista_output <- list()

# Definisci i valori di n_occ desiderati
# valori_n_occ <- c(80, 100, 200, 300, 400)

# Esegui un ciclo for per calcolare pippo per ciascun valore di n_occ
#for (i in seq_along(valori_n_occ)) {
#  pluto <- pippo(occurrences_values, valori_n_occ[i])
#  nome_variabile <- paste0("pluto", i)
#  lista_output[[nome_variabile]] <- pluto
# }


####################################################################
# Inizializza un vettore per salvare gli ipervolumi
# ipervolumi <- numeric()
# num_occurrences <- numeric()

# Funzione per calcolare l'ipervolume
# calcola_ipervolume <- function(data) {
#  hv_occ <- hypervolume_gaussian(data)
#  return(hv_occ@Volume)
# }
# end <- as.numeric(nrow(occurrences_values))
# end

# Ciclo for per incrementare il numero di occorrenze
# for (i in seq(10, end, 100)) {
  # Seleziona le prime 'i' occorrenze
#  subset_occ <- occurrences_values[1:i, ]
  
  # Calcola l'ipervolume per il subset di occorrenze
 #  hv <- calcola_ipervolume(subset_occ)
  
  # Salva l'ipervolume e il numero di occorrenze
  # ipervolumi <- c(ipervolumi, hv)
  # num_occurrences <- c(num_occurrences, i)
# }

# Dataframe
# df <- data.frame(n_occurrences = num_occurrences, hypervolume = ipervolumi)
# Stampare il dataframe
# print(df)


# Output List 
## Bootstrap con incremento
## Regressione polinomiale 
