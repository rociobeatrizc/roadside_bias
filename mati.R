# report of the data in area of applicability and DI

library(CAST)
library(caret)
library(sf)
library(devtools)
library(raster)
library(viridis)
library(ggplot2)
library(tidyverse)
library(broom)
library(terra)
library(lme4)
library(performance)
library(mgcv)
library(ggeffects)

setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo/virtualspecies_output")

# function for pixel stat
analyze_AOA <- function(raster_stack, id) {
  
  # raster AOA
  aoa_null <- raster_stack$AOA_null
  aoa_biased <- raster_stack$AOA_biased
  aoa_all <- raster_stack$AOA_all
  
  # pixels 
  common_pixels <- sum(aoa_null[] == 1 & aoa_biased[] == 1, na.rm = TRUE)
  exclusive_null <- sum(aoa_null[] == 1 & aoa_biased[] == 0, na.rm = TRUE)
  exclusive_biased <- sum(aoa_null[] == 0 & aoa_biased[] == 1, na.rm = TRUE)
  count_all <- sum(aoa_all[] == 1, na.rm = TRUE)
  
  # data frame
  result <- data.frame(
    ID = id,
    Common_Pixels = common_pixels,
    Exclusive_Null = exclusive_null,
    Exclusive_Biased = exclusive_biased,
    Count_All = count_all
  )
  
  return(result)
}

analyze_DI <- function(raster_stack, species_id) {
  
  # DI null 
  # DI biased 
  # DI all
  di_null <- raster_stack$DI_null
  di_biased <- raster_stack$DI_biased
  di_all <- raster_stack$DI_all
  
  # total number of pixels in di_all
  total_pixels_di_all <- sum(!is.na(values(di_all)))
  
  # difference
  difference_DI <- di_null - di_biased
  difference_DI_all_biased <- di_all - di_biased
  difference_DI_all_null <- di_all - di_null
  
  # mean difference
  mean_difference <- mean(values(difference_DI), na.rm = TRUE)
  mean_difference_DI_all_biased <- mean(values(difference_DI_all_biased), na.rm = TRUE)
  mean_difference_DI_all_null <- mean(values(difference_DI_all_null), na.rm = TRUE)
  
  # na.remove
  total_pixels <- sum(!is.na(values(difference_DI)))
  
  # DI null > DI biased
  pixels_null_greater <- sum(values(di_null) > values(di_biased), na.rm = TRUE)
  
  # DI_null < DI_biased
  pixels_biased_greater <- sum(values(di_null) < values(di_biased), na.rm = TRUE)
  
  # DI_null == DI_biased
  unchanged_pixels <- sum(values(di_null) == values(di_biased), na.rm = TRUE)
  
  # mean and sd
  mean_null <- mean(values(di_null), na.rm = TRUE)
  mean_biased <- mean(values(di_biased), na.rm = TRUE)
  mean_all <- mean(values(di_all), na.rm = TRUE)
  sd_null <- sd(values(di_null), na.rm = TRUE)
  sd_biased <- sd(values(di_biased), na.rm = TRUE)
  sd_all <- sd(values(di_all), na.rm = TRUE)
  
  # Calculate ddt and ddb
  ddt <- di_all - di_null
  ddb <- di_all - di_biased
  
  # Count pixels where ddt and ddb are equal to zero
  pixels_ddt_zero <- sum(values(ddt) == 0, na.rm = TRUE)
  pixels_ddb_zero <- sum(values(ddb) == 0, na.rm = TRUE)
  
  # Count pixels where ddt and ddb are positive
  pixels_ddt_positive <- sum(values(ddt) >  0, na.rm = TRUE)
  pixels_ddb_positive <- sum(values(ddb) > 0, na.rm = TRUE)
  
  # DI all less than null/biased
  all_great_than_null_2 <- all(values(ddt) <  0, na.rm = TRUE)
  all_great_than_biased_2 <- all(values(ddb) < 0, na.rm = TRUE)
  

  # df
  result <- data.frame(
    Species_ID = species_id,
    Total_pixels_DI_all = total_pixels_di_all,
    Mean_DI_null = mean_null,
    SD_DI_null = sd_null,
    Mean_DI_biased = mean_biased,
    SD_DI_biased = sd_biased,
    Mean_DI_all = mean_all,
    SD_DI_all = sd_all,
    Pixels_null_greater = pixels_null_greater,
    Pixels_biased_greater = pixels_biased_greater,
    Unchanged_pixels = unchanged_pixels,
    Mean_difference_DI = mean_difference,
    Mean_difference_DI_all_biased = mean_difference_DI_all_biased,
    Mean_difference_DI_all_null = mean_difference_DI_all_null,
    Pixels_ddt_zero = pixels_ddt_zero,
    Pixels_ddb_zero = pixels_ddb_zero,
    Pixels_ddb_positive = pixels_ddb_positive,
    Pixels_ddt_positive = pixels_ddt_positive,
    all_great_than_null = all_great_than_null_2,
    all_great_than_biased = all_great_than_biased_2
    
  )
  
  return(result)
}


# process all the rasters in the folder with AOA and DI stats
process_rasters <- function(folder = ".") {
  # files with _aoa
  raster_files <- list.files(folder, pattern = "_aoa\\.tif$", full.names = TRUE)
  
  # df
  final_results <- data.frame()
  
  for (file in raster_files) {
    # ID species
    base_name <- tools::file_path_sans_ext(basename(file))
    id <- sub("_aoa$", "", base_name)
    
    # species prevalence and number of occurrences 
    matches <- regmatches(id, regexec("sp_prevalence_([0-9.]+).*n_occ_([0-9]+)", id))
    if (length(matches[[1]]) > 2) {
      sp_prev <- as.numeric(matches[[1]][2])
      n_occ <- as.integer(matches[[1]][3])
    } else {
      sp_prev <- NA
      n_occ <- NA
    }
    
    # upload stack
    raster_stack <- rast(file)
    
    # AOA
    aoa_results <- analyze_AOA(raster_stack, id)
    aoa_results$Species_Prevalence <- sp_prev
    aoa_results$N_Occ <- n_occ
    
    # DI
    di_results <- analyze_DI(raster_stack, id)
    di_results$Species_Prevalence <- sp_prev
    di_results$N_Occ <- n_occ
    
    # same df
    combined_results <- cbind(aoa_results, di_results[, -1])
    final_results <- rbind(final_results, combined_results)
  }
  
  # save in csv
  write.csv(final_results, file = "aoa_di.csv", row.names = FALSE)
  return(final_results)
}

# use the function
csvs <- process_rasters()
head(csvs)


# upload vs and get the number of occurrences (unbiased, biased) from there
merge_with_species <- function() {
  # aoa and di
  aoa_data <- read.csv("aoa_di.csv")
  
  # original csvs with the species
  csv_files <- list.files(pattern = "^species_.*\\.csv$")
  
  # unbiased and biased points 
  calculate_occurrences <- function(csv_file) {
    # ID
    file_id <- tools::file_path_sans_ext(basename(csv_file))
    
    # read file
    species_data <- read.csv(csv_file)
    
    # split between biased and unbiased
    biased <- species_data[species_data$BIASED == TRUE, ]
    unbiased_random <- species_data[species_data$UNBIASED == TRUE & species_data$BIASED == FALSE, ]
    
    # unbiased + 20%
    n_random_to_add <- ceiling(0.2 * nrow(biased))
    set.seed(123)
    random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_random_to_add), ]
    
    # combine
    unbiased_20 <- rbind(biased, random_points)
    
    # return ID
    data.frame(
      ID = file_id,
      Occurrences_Biased = nrow(biased),
      Occurrences_Unbiased_20 = nrow(unbiased_20)
    )
  }
  
  # lapply for each species
  occurrence_data <- do.call(rbind, lapply(csv_files, calculate_occurrences))
  
  # merge occurrence data with aoa-di 
  final_data <- merge(aoa_data, occurrence_data, by.x = "ID", by.y = "ID", all.x = TRUE)
  
  # output
  write.csv(final_data, "output_final_csv_file.csv", row.names = FALSE)
  #  cat("File salvato: output_final_csv_file.csv\n")
  
  return(final_data)
}


# apply function
csvs_1 <- merge_with_species()
head(csvs_1)
names(csvs_1)

# useless columns
csvs_1$N_Occ.1 <- NULL
csvs_1$Species_Prevalence.1 <- NULL

# % biased occurrences
csvs_1$Occurrences_Biased_norm <- csvs_1$Occurrences_Biased / csvs_1$N_Occ

# % null occurrences
csvs_1$Occurrences_null_norm <- csvs_1$Occurrences_Unbiased_20 / csvs_1$N_Occ

# ddb zero norm
csvs_1$ddb_zero_norm <- csvs_1$Pixels_ddb_zero / csvs_1$Total_pixels_DI_all

# ddb positive norm
csvs_1$ddb_postive_norm <- csvs_1$Pixels_ddb_positive /csvs_1$Total_pixels_DI_all


################################## Models ######################################
# GLM
glm_ddb <- glm(ddb_zero_norm ~ Occurrences_Biased_norm * Species_Prevalence, 
               data = csvs_1, 
               family = quasibinomial)
summary(glm_ddb)

# GAM
gam_ddb <- gam((ddb_zero_norm ~ Occurrences_Biased_norm * Species_Prevalence),
               data = csvs_1)
summary(gam_ddb)


### residuals ###
# GAM
plot(resid(gam_ddb) ~ predict(gam_ddb))
abline(h = 0, col = "red")

# GLM: ok
plot(resid(glm_ddb) ~ predict(glm_ddb))
abline(h = 0, col = "red")


## histogram
# GAM
hist(resid(gam_ddb), 
     main = "gam", 
     xlab = "residuals", 
     breaks = 30)

# GLM
hist(resid(glm_ddb), 
     main = "lm", 
     xlab = "residuals", 
     breaks = 30)



## Q-Q Plot
# GAM
qqnorm(resid(gam_ddb), main = "Q-Q Plot (GAM)")
qqline(resid(gam_ddb), col = "red")

# GLM
qqnorm(resid(glm_ddb), main = "Q-Q Plot (GLM)")
qqline(resid(glm_ddb), col = "red")

install.packages("scales")
library(scales)


## ggpredict 1
# glm: smooth plot
ggpredict(glm_ddb, terms = c("Occurrences_Biased_norm [all]", "Species_Prevalence")) %>% 
  plot(., show_ci = TRUE, ci_style = "ribbon",
       show_residuals_line = TRUE,
       colors = c("purple", "green3", "orange"),
       data_labels = TRUE) + 
  labs(title = "GLM") +
  xlab("occurrences biased normalized") + 
  ylab("pixels with ddb = 0 normalized") +
  guides(color = guide_legend(title = "Species Prevalence")) +
  scale_y_continuous(labels = label_percent(accuracy = 0.001, scale = 1))
# gam
ggpredict(gam_ddb, terms = c("Occurrences_Biased_norm", "Species_Prevalence")) %>% 
  plot()

## ggpredict 2

# glm
ggpredict(glm_ddb, terms = c("Species_Prevalence", "Occurrences_Biased_norm")) %>% 
   plot() +
  scale_y_continuous(labels = label_percent(accuracy = 0.001, scale = 1))


## ggeffects 

# predict response: understand results
# test_predictions: significant results
# plot: figures



# predict response
result <- predict_response(glm_ddb,terms =  c("Occurrences_Biased_norm [all]", "Species_Prevalence")) 
print(result, collapse_table = TRUE, collapse_ci = TRUE, n = Inf)
plot(result)

test_predictions(result)

ggpredict(glm_ddb,terms =  c("Species_Prevalence", "Occurrences_Biased_norm"))

#########################
## validazione glm

# Prevedi probabilità
predictions <- predict(glm_ddb, type = "response")


# 2. Deviance e Pseudo-R² (McFadden's R²)
model_deviance <- deviance(glm_ddb)
null_deviance <- deviance(glm(ddb_zero_norm ~ 1, data = csvs_1, family = quasibinomial))
pseudo_r_squared <- 1 - model_deviance / null_deviance
print(paste("Deviance del modello:", model_deviance))
print(paste("Deviance del modello nullo:", null_deviance))
print(paste("McFadden's R²:", pseudo_r_squared))

# 3. Cross-validation (10-fold)
train_control <- trainControl(method = "cv", number = 10)
cv_model <- train(ddb_zero_norm ~ Occurrences_Biased_norm * Species_Prevalence, 
                  data = csvs_1, 
                  method = "glm", 
                  family = quasibinomial, 
                  trControl = train_control)
print(cv_model)
dispersion <- summary(glm_ddb)$dispersion
print(paste("Dispersion parameter:", dispersion))


# Calcolo dei residui
residuals_glm <- residuals(glm_ddb, type = "response")

# Calcolo del RMSE
rmse_residuals <- sqrt(mean(residuals_glm^2))

# Output del RMSE
print(paste("RMSE dei residui:", rmse_residuals))

##############################################
# species prevalence, occurrences biased
library(viridis)

ggplot(csvs_1, aes(x = as.factor(Species_Prevalence), y = Occurrences_Biased_norm, fill = as.factor(Species_Prevalence))) +
  geom_boxplot(alpha = 0.8, outlier.color = "red", outlier.shape = 16, color = "black") +
  scale_fill_viridis_d(option = "viridis") +  # Usa viridis per i colori
  labs(x = "species prevalence", 
       y = "normalized biased occurrences", 
       title = "Distribution of biased occurrences in relation to species prevalence") +
  theme_minimal(base_size = 10) +  # Stile pulito con testo più leggibile
  theme(legend.position = "none",  # Nasconde la legenda perché il colore è ridondante
        plot.title = element_text(face = "bold", hjust = 0.5))  # Centra e grassetta il titolo

library(car)
leveneTest(Occurrences_Biased_norm ~ as.factor(Species_Prevalence), data = csvs_1)


library(dplyr)
var_data <- csvs_1 %>%
  group_by(Species_Prevalence) %>%
  summarise(variance = var(Occurrences_Biased_norm))

# Modello lineare o non lineare sulla varianza
lm_var <- lm(variance ~ Species_Prevalence, data = var_data)
summary(lm_var)

ggpredict(lm_var) %>% 
  plot(., show_ci = TRUE, ci_style = "ribbon",
       show_residuals_line = TRUE,
       colors = c("darkgreen"),
       data_labels = TRUE) +
  labs(title = "variance of biased occurrences as the species prevalence increases") +
  xlab("species prevalence") + 
  ylab("variance of biased occurrences") 
