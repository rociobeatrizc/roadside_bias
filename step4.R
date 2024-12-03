# report of the data in area of applicability and DI

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

# DI
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

coords <- difference_DI_df %>% dplyr::select(Longitude = x, Latitude = y)

data.wgs <- SpatialPointsDataFrame(
  coords = coords,
  data = difference_DI_df,
  proj4string = CRS("+proj=eck4")
)

mod0 = mgcv::gam(DI_null ~ 1 + s(Longitude, Latitude, bs = "sos"), method = "REML", data = data.wgs)
summary(mod0)


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

# what should I save from this script? 
