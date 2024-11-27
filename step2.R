# hypervolume as occurrences increase
# original sampling, unbiased and biased

library(sf)
library(raster)
library(virtualspecies)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(terra) 
library(geodata)
library(osmdata)
library(osmextract)
library(hypervolume)

setwd("C:/Users/rocio/Desktop/PHD/1 year/Abruzzo")


# let's start from the csv we created before (one of them)
# this is the species n. 2, with species prevalence 0.3, sample prevalence 0.9 and 200 random sampled occurrences
sp2_sp_prev0.3_sample_prev0.9_nocc100 <- read.csv("species_2_sp_prevalence_0.3_sample_prev_0.9_n_occ_100.csv",
                                                  row.names = NULL)

typeof(sp2_sp_prev0.3_sample_prev0.9_nocc100)
# csv check 
head(sp2_sp_prev0.3_sample_prev0.9_nocc100)
nrow(sp2_sp_prev0.3_sample_prev0.9_nocc100)

# split the dataset in unbiased and biased
unbiased_random <- sp2_sp_prev0.3_sample_prev0.9_nocc100[sp2_sp_prev0.3_sample_prev0.9_nocc100$UNBIASED == TRUE & sp2_sp_prev0.3_sample_prev0.9_nocc100$BIASED == FALSE, ]
biased <- sp2_sp_prev0.3_sample_prev0.9_nocc100[sp2_sp_prev0.3_sample_prev0.9_nocc100$BIASED == TRUE, ]

head(unbiased_random)
head(biased)

# now we will simulate the same sampling effort, i.e. the two datasets should have more or less the
# same number of points 
# let's say that, if the number of biased points is n, the number of unbiased points will be n + 0.2 * n 
# please note that the unbiased sample must contain the biased one

# n + 0.2 * n 
n_random_to_add <- ceiling(0.2 * nrow(biased))

# random extraction of points BIASED = FALSE (UNBIASED = TRUE), without overlaps with BIASED = TRUE
set.seed(123) 
random_points <- unbiased_random[sample(1:nrow(unbiased_random), n_random_to_add), ]

# combine datasets
unbiased_20 <- rbind(biased, random_points)

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


# for each of these datasets (unbiased, biased, total occurrences)
# just the environmental variables are needed
head(biased)
head(unbiased_20)
head(sp2_sp_prev0.3_sample_prev0.9_nocc100)

# drop useless variables
drops <- c("X", "Y", "distance", "ID", "probability", "UNBIASED", "BIASED")

# subset with environmental data
biased_env <- biased[ , !(names(biased) %in% drops)]
unbiased_20_env <- unbiased_20[ , !(names(unbiased_20) %in% drops)]
sp2_env <- sp2_sp_prev0.3_sample_prev0.9_nocc100[ , !(names(sp2_sp_prev0.3_sample_prev0.9_nocc100) %in% drops)]

# number
nrow(biased_env)
nrow(unbiased_20_env)
nrow(sp2_env)


# let's build the curve that describes the hypervolume as the occurrences accumulate and increase
# the following function sets the starting number of points and the steps for the hypervolume building

define_hyp_steps <- function(data) {
  # total points
  total_points <- nrow(data)
  
  # initial value (20% points)
  from <- ceiling(0.2 * total_points)
  
  # each steps is equal to the 25% of the points
  by <- ceiling(0.25 * total_points)
  
  # generate the sequence
  c(seq(from = from, to = total_points, by = by), total_points)
}

# x_seq is needed for the curve building with the LOESS method
x_seq <- 1:(nrow(sp2_env))

# number of simulations
num_sim <- 2

# hypervolume building
process_hypervolume <- function(env_data, x_seq, num_sim, acc_curve_func, hyp_steps = NULL) {
  
  if (is.null(hyp_steps)) {
    hyp_steps <- define_hyp_steps(env_data)
  }
  
  # simulations
  all_simulations <- lapply(1:num_sim, function(sim) {
    lapply(hyp_steps, function(step) {
      d_hyp <- acc_curve_func(env_data, step)
      d_hyp[[1]]
    })
  })
  
  # combine all results in one df
  combined_df <- do.call(rbind, lapply(seq_along(all_simulations), function(sim) {
    do.call(rbind, lapply(all_simulations[[sim]], function(df) {
      df$sim <- sim
      df
    }))
  }))
  
  # mean values for each value of x_seq with LOESS
  pred_mean <- lapply(x_seq, function(n) {
    preds <- sapply(all_simulations, function(lista) {
      loess_fit <- loess(iperv ~ n_occ, data = do.call(rbind, lista))
      predict(loess_fit, newdata = data.frame(n_occ = n))
    })
    data.frame(n_occ = n, iperv_mean = mean(preds, na.rm = TRUE))
  })
  
  # one df
  pred_mean_df <- do.call(rbind, pred_mean)
  
  # graphic
  plot <- ggplot() +
    geom_smooth(data = combined_df, aes(x = n_occ, y = iperv, group = sim), 
                method = "loess", se = FALSE, color = "grey", size = 0.5, alpha = 0.5) +
    geom_line(data = pred_mean_df, aes(x = n_occ, y = iperv_mean), 
              color = "sienna1", size = 1.2) +
    labs(title = "Mean Hypervolume", x = "Occurrences", y = "Hypervolume") +
    theme_minimal()
  
  # output: graph and df
  return(list(data = pred_mean_df, plot = plot))
}


# unbiased
result_unbiased <- process_hypervolume(
  unbiased_20_env, 
  x_seq, 
  num_sim = num_sim, 
  acc_curve_func = acc_curve 
)

# biased
result_biased <- process_hypervolume(
  biased_env, 
  x_seq, 
  num_sim = num_sim, 
  acc_curve_func = acc_curve, 
)

# all
result_all <- process_hypervolume(
  sp2_env, 
  x_seq, 
  num_sim = num_sim, 
  acc_curve_func = acc_curve
)


# results: dataframe and plot
# unbiased
print(result_unbiased$plot)
print(result_unbiased$data)

# biased
print(result_biased$plot)
print(result_biased$data)

# all
print(result_all$plot)
print(result_all$data)


# plot all together
data_unbiased <- result_unbiased$data
data_unbiased$type <- "Unbiased"

data_biased <- result_biased$data
data_biased$type <- "Biased"

data_all <- result_all$data
data_all$type <- "All"

# all dfs
combined_data <- result_unbiased$data %>%
  mutate(type = "Unbiased") %>%
  bind_rows(
    result_biased$data %>% mutate(type = "Biased"),
    result_all$data %>% mutate(type = "All")
  )



# plot data
plot_hyp <- ggplot(combined_data, aes(x = n_occ, y = iperv_mean, color = type)) +
  geom_line(size = 1) +
  labs(
    title = "comparison",
    x = "n_occ",
    y = "iperv_mean",
    color = "type"
  ) +
theme_minimal()


print(plot_hyp)

# wide
wide_data <- combined_data %>%
pivot_wider(names_from = type, values_from = iperv_mean, names_prefix = "iperv_mean_")


# original dataframe: we need to copy the infos in the new dataframe
dataset_name <- deparse(substitute(sp2_sp_prev0.3_sample_prev0.9_nocc100))

# as csv
output_filename <- paste0(dataset_name, "_hypervolumes.csv")
write.csv(wide_data, file = output_filename, row.names = FALSE)


# save graph
output_file_plot <- paste0(dataset_name, "_hypervolume_plot.jpeg")
ggsave(output_file_plot, plot = plot_hyp, device = "jpeg", width = 8, height = 6, dpi = 300)
