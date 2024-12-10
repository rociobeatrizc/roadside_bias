# report of the data in area of applicability and DI
# just one csv

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

setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo")

# upload
sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa <-  rast("sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa.tif")
plot(sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa)

# select aoa from the stack
sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa


### AOA
aoa_null   <- sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa$AOA_null
aoa_biased <- sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa$AOA_biased
aoa_all    <- sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa$AOA_all


# common pixels between null and biased
common_pixels <- sum(aoa_null[] == 1 & aoa_biased[] == 1, na.rm = TRUE)
dev.off()

par(mfrow=c(2,1))
plot(aoa_null)
plot(aoa_biased)

# exclusive null
exclusive_null <- sum(aoa_null[] == 1 & aoa_biased[] == 0, na.rm = TRUE)

# exclusive biased
exclusive_biased <- sum(aoa_null[] == 0 & aoa_biased[] == 1, na.rm = TRUE)

# pixels of all
count_all <- sum(aoa_all[] == 1, na.rm = TRUE)

# results
common_pixels
exclusive_null
exclusive_biased
count_all

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


# stats
aoa_analysis <- analyze_AOA(sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa, id = "sp2_0.3_0.9_350")
print(aoa_analysis)



### DI
di_null   <- sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa$DI_null
di_biased <- sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa$DI_biased

par(mfrow=c(2,1))
plot(di_null)
plot(di_biased)

# DI difference
difference_DI <- di_null - di_biased

# negative values: DI_biased > DI_null
summary(difference_DI)

# dataframe
difference_DI_df <- as.data.frame(difference_DI, xy = TRUE)

# plot
ggplot(difference_DI_df, aes(x = x, y = y, fill = DI_null)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", high = "blue", mid = "white", midpoint = 0,
    limits = c(min(difference_DI[]), max(difference_DI[])),
    name = "Difference in DI"
  ) +
  theme_minimal() +
  coord_equal() +
  labs(title = "Difference in DI between AOA_null and AOA_biased")


# DI stats
analyze_DI <- function(raster_stack, species_id) {
 
  # DI null DI biased
  di_null <- raster_stack$DI_null
  di_biased <- raster_stack$DI_biased
  
  # difference
  difference_DI <- di_null - di_biased
  
  # na.remove
  total_pixels <- sum(!is.na(values(difference_DI)))
  
  # DI null > DI biased
  pixels_null_greater <- sum(values(di_null) > values(di_biased), na.rm = TRUE)
  
  # DI_null < DI_biased
  pixels_biased_greater <- sum(values(di_null) < values(di_biased), na.rm = TRUE)
  
  # Pixel invariati (DI_null == DI_biased)
  unchanged_pixels <- sum(values(di_null) == values(di_biased), na.rm = TRUE)
  
  # %
  percent_null_greater <- (pixels_null_greater / total_pixels) * 100
  percent_biased_greater <- (pixels_biased_greater / total_pixels) * 100
  percent_unchanged <- (unchanged_pixels / total_pixels) * 100
  
  # mean and sd
  mean_null <- mean(values(di_null), na.rm = TRUE)
  mean_biased <- mean(values(di_biased), na.rm = TRUE)
  sd_null <- sd(values(di_null), na.rm = TRUE)
  sd_biased <- sd(values(di_biased), na.rm = TRUE)
  
  # data frame
  result <- data.frame(
    Species_ID = species_id,
    Mean_DI_null = mean_null,
    SD_DI_null = sd_null,
    Mean_DI_biased = mean_biased,
    SD_DI_biased = sd_biased,
    Pixels_null_greater = pixels_null_greater,
    Percent_null_greater = percent_null_greater,
    Pixels_biased_greater = pixels_biased_greater,
    Percent_biased_greater = percent_biased_greater,
    Unchanged_pixels = unchanged_pixels,
    Percent_unchanged = percent_unchanged
  )
  
  return(result)
}

# example
species_id <- "sp2_0.3_0.9_350"  # ID della specie
di_analysis <- analyze_DI(sp2_sp_prev0.3_sample_prev0.9_nocc350_aoa, species_id)
print(di_analysis)

# csv
write.csv(di_analysis, "DI_analysis_species.csv", row.names = FALSE)



# georg lm
coords <- difference_DI_df %>% dplyr::select(Longitude = x, Latitude = y)

data.wgs <- SpatialPointsDataFrame(
  coords = coords,
  data = difference_DI_df,
  proj4string = CRS("+proj=eck4")
)

mod0 = mgcv::gam(DI_null ~ 1 + s(Longitude, Latitude, bs = "sos"), method = "REML", data = data.wgs)
summary(mod0)xx


# https://jakubnowosad.com/posts/2024-10-13-spatcomp-bp1/
# from raster to vector: how to show the dissimilarity info

mod = lm(values(di_biased, na.rm = TRUE) ~ values(di_null, na.rm = TRUE))
summary(mod)

ggplot() +
 geom_point(aes(x = values(di_biased, na.rm = TRUE), y = values(di_null, na.rm = TRUE))) +
 geom_smooth(aes(x = values(di_biased, na.rm = TRUE), y =  values(di_biased, na.rm = TRUE) - values(di_null, na.rm = TRUE)),
             method = "lm", color = "darkgreen") +
 geom_smooth(aes(x = values(di_biased, na.rm = TRUE), y = values(di_null, na.rm = TRUE)),
              method = "lm") +
 geom_smooth(aes(x = c(0:7), y = c(0:7)), color = "red", method = "lm") + 
 theme_minimal() +
 labs(x = "DI biased", y = "DI null")
