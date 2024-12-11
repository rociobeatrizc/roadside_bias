# report of the data in area of applicability and DI
# for all the stacks generated before (script3_process.R)

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

# DI stats
analyze_DI <- function(raster_stack, species_id) {
  
  # DI null 
  # DI biased 
  # DI all
  di_null <- raster_stack$DI_null
  di_biased <- raster_stack$DI_biased
  di_all <- raster_stack$DI_all
  
  # difference
  difference_DI <- di_null - di_biased
  
  # mean difference
  mean_difference <- mean(values(difference_DI), na.rm = TRUE)
  
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
  
  # relationship bewteen all, null, biased: TRUE/FALSE
  all_less_than_null <- mean_all < mean_null
  all_less_than_biased <- mean_all < mean_biased
  
  # df
  result <- data.frame(
    Species_ID = species_id,
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
    All_less_than_null = all_less_than_null,
    All_less_than_biased = all_less_than_biased
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

# how many times ALL is greater than NULL?
count_null_false <- sum(csvs$All_less_than_null == FALSE)

# how many times ALL is greater than BIASED
count_biased_false <- sum(csvs$All_less_than_biased == FALSE)

# print
cat("All_less_than_null FALSE count:", count_null_false, "\n")
cat("All_less_than_biased FALSE count:", count_biased_false, "\n")


##### georg lm
# coords <- difference_DI_df %>% dplyr::select(Longitude = x, Latitude = y)

# data.wgs <- SpatialPointsDataFrame(
#  coords = coords,
#  data = difference_DI_df,
#  proj4string = CRS("+proj=eck4")
# )

# mod0 = mgcv::gam(DI_null ~ 1 + s(Longitude, Latitude, bs = "sos"), method = "REML", data = data.wgs)
summary(mod0)

# https://jakubnowosad.com/posts/2024-10-13-spatcomp-bp1/
# from raster to vector: how to show the dissimilarity info

# mod = lm(values(di_biased, na.rm = TRUE) ~ values(di_null, na.rm = TRUE))
# summary(mod)


# ggplot() +
#  geom_point(aes(x = values(di_biased, na.rm = TRUE), y = values(di_null, na.rm = TRUE))) +
#  geom_smooth(aes(x = values(di_biased, na.rm = TRUE), y =  values(di_biased, na.rm = TRUE) - values(di_null, na.rm = TRUE)),
#              method = "lm", color = "darkgreen") +
#  geom_smooth(aes(x = values(di_biased, na.rm = TRUE), y = values(di_null, na.rm = TRUE)),
#              method = "lm") +
#  geom_smooth(aes(x = c(0:7), y = c(0:7)), color = "red", method = "lm") + 
#  theme_minimal() +
#  labs(x = "DI biased", y = "DI null")



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
