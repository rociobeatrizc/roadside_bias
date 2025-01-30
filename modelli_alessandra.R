library(maps)
library(terra)
library(ggplot2)
library(tidyverse)
library(ggeffects)
library(lme4)
library(performance)
library(mgcv)

setwd("C:/Users/rocio/Downloads")

d <- read.csv("all_species_combined_results.csv")
# %>% 
#  pivot_longer(., cols = c("D", "I"), names_to = "Metric", values_to = "Value")


# prima: relazione fra species prevalence e biased ratio
hist(d$BiasedRatio, breaks = 40)

# rapporto fra i due (prova lm, glm, gam)
lm_d <- glm(BiasedRatio ~  SpeciesPrevalence, 
            family = quasibinomial, 
            data = d)

summary(lm_d)           
plot(residuals(lm_d))

# RMSE
residuals_glm <- residuals(lm_d, type = "response")

# Calcolo del RMSE
rmse_residuals <- sqrt(mean(residuals_glm^2))

# Output del RMSE
print(paste("RMSE dei residui:", rmse_residuals))


# relazione fra delta d, biased ratio e species prevalence
# deltaDI, biased ratio e species prevalence
glm_d <- glm(DeltaD ~ BiasedRatio + SpeciesPrevalence, 
              family = quasibinomial, 
              data = d)

summary(glm_d)

plot(residuals(glm_d))


# ggeffects 
ggpredict(glm_d, terms = c("BiasedRatio[all]", "SpeciesPrevalence")) %>% 
  plot() 


ggpredict(glm_d, terms = c("SpeciesPrevalence", "BiasedRatio"))  %>% 
  plot()
