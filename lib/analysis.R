
## Libraries -----------------------------------------------------------

library(tidyverse)
source("lib/model/phenoBayes.R")
source("lib/model/get_parameter.R")

## Parameters ----------------------------------------------------------

foresttype <- "broadleaved"
#foresttype <- "coniferous"
sampsize <- 25

## Data ----------------------------------------------------------------

load(paste0("data/landsat_", foresttype, ".RData"))
load(paste0("data/climate_dwd_spring.RData"))

temporal_drivers <- as.matrix(climate_vars[, c("preseason_temp_anomaly", "chilling_anomaly")])
temporal_drivers <- cbind(temporal_drivers, interaction = temporal_drivers[, 1] * temporal_drivers[, 2])

## Sampling ------------------------------------------------------------

select_pixels <- seq(1, 250, sampsize)

input <- vector("list", length = length(select_pixels))

for (i in 1:length(select_pixels)) {
  
  dat_tmp <- dat %>% filter(pixel %in% select_pixels[i]:(select_pixels[i] + (sampsize - 1)))
  dat_tmp$pixel <- as.integer(factor(dat_tmp$pixel, labels = 1:sampsize))
  spatial_drivers_tmp <- spatial_drivers[unique(dat_tmp$pixel), ]
  spatial_drivers_tmp <- data.frame(elevation = spatial_drivers_tmp)
  
  input[[i]] <- list(pixel = dat_tmp$pixel,
                     index = dat_tmp$evi, 
                     doy = dat_tmp$doy, 
                     year = dat_tmp$year,
                     spatial_drivers = spatial_drivers_tmp)
}

## Models ----------------------------------------------------------------

mod <- vector("list", length = length(input))

for (i in 1:length(mod)) {
  mod[[i]] <- phenoBayes(pixel = input[[i]]$pixel,
                         index = input[[i]]$index, 
                         doy = input[[i]]$doy, 
                         year = input[[i]]$year,
                         spatial_drivers = NULL, 
                         temporal_drivers = temporal_drivers
                         #contro = list(adapt_delta = 0.99)
                         )
}

## Ensemble development of rho -----------------------------------------

rho_development <- vector("list", length(mod))

for (i in 1:length(mod)) {
  tmp <- apply(simplify2array(map(mod[1:i], ~ .$rho)), 1:2, mean)
  rho_development[[i]] <- get_parameter(list(rho = tmp), "rho", p_lower = 0.025, p_upper = 0.975)
}

rho_development <- rho_development %>%
  set_names(as.character(1:length(mod))) %>%
  bind_rows(.id = "permutations") %>%
  mutate(permutations = as.integer(permutations)) %>%
  mutate(predictor = factor(predictor, labels = c("Temperature", "Chilling", "Temperature x Chilling")))

write_csv(rho_development, paste0("results/rho_development_", foresttype, ".csv"))

## Extract ensemble estimate ------------------------------------------------

rho <-  apply(simplify2array(map(mod, ~ .$rho)), 1:2, mean)
lambda <- apply(do.call("cbind", map(mod, ~ .$lambda)), 1, mean)

phi1 <- apply(simplify2array(map(mod, ~ .$phi[[1]])), 1:2, mean)
phi2 <- apply(simplify2array(map(mod, ~ .$phi[[2]])), 1:2, mean)
phi3 <- apply(simplify2array(map(mod, ~ .$phi[[3]])), 1:2, mean)
phi4 <- apply(simplify2array(map(mod, ~ .$phi[[4]])), 1:2, mean)

beta1 <- apply(simplify2array(map(mod, ~ .$beta[[1]])), 1:2, mean)
beta2 <- apply(simplify2array(map(mod, ~ .$beta[[2]])), 1:2, mean)
beta3 <- apply(simplify2array(map(mod, ~ .$beta[[3]])), 1:2, mean)
beta4 <- apply(simplify2array(map(mod, ~ .$beta[[4]])), 1:2, mean)

cor_beta <- apply(simplify2array(map(mod, ~ .$cor_beta)), 1:2, mean)
cor_phi <- apply(simplify2array(map(mod, ~ .$cor_phi)), 1:2, mean)

mod_ensemble <- list(data = list(year = 1985:2015, pixel = 1:10), 
                     phi = list(phi1, phi2, phi3, phi4),
                     beta = list(beta1, beta2, beta3, beta4),
                     rho = rho,
                     lambda = lambda,
                     cor_beta = cor_beta,
                     cor_phi = cor_phi)

save(mod_ensemble, file = paste0("results/mod_ensemble_", foresttype, ".RData"))
load(file = paste0("results/mod_ensemble_", foresttype, ".RData"))

## Extract Phi ----------------------------------------------------------------

phi_ensemble_estimate <- get_parameter(mod_ensemble, "phi", p_lower = 0.025, p_upper = 0.975) %>%
  mutate(parameter = factor(parameter, labels = c("Minimum", "Maxmium", "Gree-up rate", "Start of season")))

write_csv(phi_ensemble_estimate, paste0("results/phi_ensemble_estimate_", foresttype, ".csv"))

## Validation --------------------------------------------------------------

iter <- nrow(mod_ensemble$rho)

### In-situ measurements

unfolding_dat <- readxl::read_excel("data/phenology/pheno_ipg_waldhaeuser_schoenbrunn.xlsx") %>%
  dplyr::filter(year %in% c(1985:2015) & pft == foresttype) %>%
  mutate(year = as.integer(year))

unfolding <- unfolding_dat %>%
  group_by(year) %>%
  summarize(unfolding_mean = mean(unfolding, na.rm = TRUE),
            unfolding_se = sd(unfolding, na.rm = TRUE) / sqrt(length(unfolding))) %>%
  ungroup() %>%
  mutate(unfolding_deviance_mean = unfolding_mean - mean(unfolding_mean, na.rm = TRUE),
         unfolding_deviance_lower = unfolding_deviance_mean - 2 * unfolding_se,
         unfolding_deviance_upper = unfolding_deviance_mean + 2 * unfolding_se)

### Phi

validation_posterior <- mod_ensemble$phi[[4]] %>%
  as_data_frame(.) %>%
  gather(key = year, value = sos) %>%
  mutate(year = as.integer(rep(1985:2015, each = iter)),
         iter = rep(1:iter, length(1985:2015))) %>%
  group_by(iter) %>%
  mutate(sos_deviance = sos - mean(sos)) %>%
  ungroup() %>%
  left_join(unfolding, by = "year")

### Prediction based on rho

pred <- vector("list", iter)

for (i in 1:iter) {
  pred[[i]] <- temporal_drivers %*% mod_ensemble$rho[i,]
}

validation_posterior <- t(do.call("cbind", pred)) %>%
  as_data_frame(.) %>%
  gather(key = year, value = sos_prediction) %>%
  mutate(year = rep(1985:2015, each = iter),
         iter = rep(1:iter, length(1985:2015))) %>%
  group_by(iter) %>%
  mutate(sos_prediction_deviance = sos_prediction - mean(sos_prediction)) %>%
  ungroup() %>%
  left_join(validation_posterior, by = c("year", "iter"))

### Calculate validation measures

validation_dat <- validation_posterior %>%
  group_by(year) %>%
  summarize(sos_deviance_median = median(sos_deviance),
            sos_deviance_lower = quantile(sos_deviance, 0.025),
            sos_deviance_upper = quantile(sos_deviance, 0.975),
            sos_prediction_median = median(sos_prediction),
            sos_prediction_lower = quantile(sos_prediction, 0.025),
            sos_prediction_upper = quantile(sos_prediction, 0.975),
            unfolding_deviance_mean = unique(unfolding_deviance_mean),
            unfolding_deviance_lower = unique(unfolding_deviance_lower),
            unfolding_deviance_upper = unique(unfolding_deviance_upper)) %>%
  ungroup()

validation <- validation_posterior %>%
  group_by(iter) %>%
  summarize(sos.r = cor(sos_deviance, unfolding_deviance_mean, use = "complete.obs"),
            sos.rmse = sqrt(mean((sos_deviance - unfolding_deviance_mean)^2, na.rm = TRUE)),
            sos_prediction.r = cor(sos_prediction, unfolding_deviance_mean, use = "complete.obs"),
            sos_prediction.rmse = sqrt(mean((sos_prediction - unfolding_deviance_mean)^2, na.rm = TRUE))) %>%
  ungroup() %>%
  gather(key = key, value = value, -iter) %>%
  separate("key", c("model", "measure"), "\\.") %>%
  spread(key = measure, value = value) %>%
  group_by(model) %>%
  summarize(r_median = median(r),
            r_lower = quantile(r, 0.025),
            r_upper = quantile(r, 0.975),
            rmse_median = median(rmse),
            rmse_lower = quantile(rmse, 0.025),
            rmse_upper = quantile(rmse, 0.975))

write_csv(validation, paste0("results/validation_", foresttype, ".csv"))
write_csv(validation_dat, paste0("results/validation_dat_", foresttype, ".csv"))

## Response temporal model -----------------------------------------------

iter <- nrow(mod_ensemble$rho)
pred <- vector("list", iter)
newdat <- expand.grid(seq(min(temporal_drivers[,1]), max(temporal_drivers[,1]), length.out = 100),
                      c(-2, 0, 2))
newdat[,3] <- newdat[,1] * newdat[,2]
colnames(newdat) <- colnames(temporal_drivers)
newdat <- as.matrix(newdat)

for (i in 1:iter) {
  pred[[i]] <- newdat %*% mod_ensemble$rho[i,]
}

pred_posterior <- t(do.call("cbind", pred))

response_curves <- as.data.frame(newdat)
response_curves$pred <- apply(pred_posterior, 2, median)
response_curves$lower <- apply(pred_posterior, 2, quantile, 0.025)
response_curves$upper <- apply(pred_posterior, 2, quantile, 0.975)

write_csv(response_curves, paste0("results/response_curves_", foresttype, ".csv"))

## Estimates temporal model -----------------------------------------------

rho_estimates <- mod_ensemble$rho %>%
  as_data_frame(.) %>%
  gather(key = year, value = value) %>%
  mutate(predictor = rep(c("Temperature", "Chilling", "Temperature x Chilling"), each = iter),
         iter = rep(1:iter, 3)) %>%
  group_by(predictor) %>%
  summarize(median = median(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(St_median = median / c(apply(climate_vars[, c("chilling", "preseason_temp")], 2, sd), 1),
         St_lower = lower / c(apply(climate_vars[, c("chilling", "preseason_temp")], 2, sd), 1),
         St_upper = upper / c(apply(climate_vars[, c("chilling", "preseason_temp")], 2, sd), 1))

write_csv(rho_estimates, paste0("results/rho_estimates_", foresttype, ".csv"))

# R hat & Posterior predictions ---------------------------------------------------

quantile(unlist(map(mod, ~ .$rhat)) , 0.99, na.rm = TRUE)

sim_summary <- vector("list", length(mod))

for (i in 1:length(mod)) {
  
  sim <- as_data_frame(t(mod[[i]]$posterior_simulations))
  sim$pixel <- mod[[i]]$data$pixel
  sim$index <- mod[[i]]$data$index
  sim$doy <- mod[[i]]$data$doy
  sim$year <- mod[[i]]$data$year
  sim <- sim %>%
    gather(key = iteration, value = value, -pixel, -doy, -year, -index) %>%
    mutate(iteration = as.integer(gsub("V", "", iteration)))
  
  sim_summary[[i]] <- sim %>%
    group_by(pixel, doy, year, index) %>%
    summarize(lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    ungroup() %>%
    mutate(within = ifelse(index > lower & index < upper, TRUE, FALSE),
           model = i)

}

sim_summary_df <- do.call("rbind", sim_summary)

write_csv(sim_summary_df, paste0("results/posterior_predictions_", foresttype, ".csv"))



