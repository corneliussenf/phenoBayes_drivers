
## Libraries -----------------------------------------------------------

library(tidyverse)

## Get data ------------------------------------------------------------

names <- c("Coniferous", "Broadleaved")
colors <- c("#009988", "#EE7733")

phi_ensemble_estimate <- list(read_csv("results/phi_ensemble_estimate_coniferous.csv"),
                              read_csv("results/phi_ensemble_estimate_broadleaved.csv")) %>% 
  set_names(names) %>%
  bind_rows(.id = "foresttype")

response_curves <- list(read_csv("results/response_curves_coniferous.csv"),
                        read_csv("results/response_curves_broadleaved.csv")) %>% 
  set_names(names) %>%
  bind_rows(.id = "foresttype")

rho_development <- list(read_csv("results/rho_development_coniferous.csv"),
                        read_csv("results/rho_development_broadleaved.csv")) %>% 
  set_names(names) %>%
  bind_rows(.id = "foresttype")

validation_dat <- list(read_csv("results/validation_dat_coniferous.csv"),
                       read_csv("results/validation_dat_broadleaved.csv")) %>% 
  set_names(names) %>%
  bind_rows(.id = "foresttype")

simulations_dat <- list(read_csv("results/posterior_predictions_coniferous.csv"),
                        read_csv("results/posterior_predictions_broadleaved.csv")) %>% 
  set_names(names) %>%
  bind_rows(.id = "foresttype")

## Example time series (Figure 1) ------------------------------------------

load(paste0("data/landsat/landsat_broadleaved.RData"))

getPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))
cols <- getPalette(max(dat$year_index)) 

#p_sel <- sample(unique(dat$pixel), 1)
p_sel <- 5

# Estimate transition

d <- subset(dat, dat$pixel == p_sel & dat$doy > 80)
d <- d[order(d$doy),] 

slope <- c()
ddoy <- c()
for(i in 21:(nrow(d) - 21)) {
  dd <- d[(i - 20):(i + 20), ]
  dd$t <- 1:nrow(dd)
  slope <- c(slope, coef(lm(dd$evi ~ dd$t))[2])
  ddoy <- c(ddoy, d[i, "doy"])
}

transition <- unlist(ddoy)[which(slope < 0 & dplyr::lead(slope) < 0 & dplyr::lead(slope, 2) < 0)[1]]
dat_spring <- subset(dat, pixel == p_sel & dat$doy <= transition & dat$doy > 80)

# Estimate parameters

nls_fit <- nls(evi ~ b1 + (b2 / (1 + exp(- b3 * (doy - b4)))),
               start = list(b1 = min(dat_spring$evi), 
                            b2 = max(dat_spring$evi), 
                            b3 = 0.15, 
                            b4 = round(mean(dat_spring[which(dat_spring$evi > median(dat_spring$evi)), "doy"]), 0)),
               data = dat_spring)
parameters <- as.data.frame(t(coef(nls_fit)))
ff <- function (doy) parameters$b1 + (parameters$b2 / (1 + exp(- parameters$b3 * (doy - parameters$b4))))
fit <- data_frame(doy = sort(unique(dat_spring$doy)), fit = ff(sort(unique(dat_spring$doy))))

p1 <- dat %>%
  filter(pixel == p_sel) %>%
  ggplot(., aes(x = lubridate::as_date(date), y = evi)) +
  geom_point(aes(col = year)) +
  theme_bw() +
  labs(x = "Year", y = "EVI", col = "Year") +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        legend.position = "none") +
  scale_colour_gradientn(colors = cols)

p2 <- dat %>%
  filter(pixel == p_sel) %>%
  ggplot(., aes(x = doy, y = evi)) +
  geom_point(col = "grey") +
  geom_point(data = dat_spring, aes(x = doy, y = evi, col = year)) +
  theme_bw() +
  labs(x = "Day of year", y = NULL, col = "Year") +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.8, "line")) +
  scale_colour_gradientn(colors = cols) +
  geom_line(data = fit, aes(x = doy, y = fit)) +
  geom_point(data = data.frame(x = parameters$b4,
                               y = ff(parameters$b4)),
             aes(x, y), size = 2)

p <- gridExtra::grid.arrange(p1, p2, widths = c(1.8, 1.2)) 

ggsave(paste0("example_pixel.pdf"), p, path = "figures/", width = 7.5, height = 2)

## Plot phi (Figure 3) --------------------------------------------------

p <- ggplot(phi_ensemble_estimate, aes(x = year, y = median,
                                       col = foresttype, shape = foresttype)) + 
  geom_point() +
  geom_line(alpha = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.3) +
  facet_wrap(~parameter, ncol = 1, strip.position = "right", scales = "free_y") +
  theme_bw() +
  labs(x = "Year", y = "Inter-annual variability", col = NULL, shape = NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        #legend.position = "right",
        legend.position = c(1, 1.015),
        legend.justification = c(1, 1),
        legend.background = element_blank()) +
  geom_hline(data = phi_ensemble_estimate %>% group_by(parameter, foresttype) %>% summarize(mean = mean(median)),
             aes(yintercept = mean, col = foresttype), linetype = "dashed", alpha = 0.3) +
  scale_color_manual(values = colors) +
  guides(col = guide_legend(ncol = 1, 
                            keywidth = 0.1,
                            keyheight = 0.1,
                            default.unit = "inch"))

ggsave(paste0("inter-annual_deviance_phi.pdf"), p, path = "figures/", width = 3.5, height = 5)

## Response curves (Figure 4) -------------------------------------------

load(paste0("data/dwd_climate/climate_dwd_spring.RData"))
temporal_drivers <- as.matrix(climate_vars[, c("preseason_temp_anomaly", "chilling_anomaly")])
temporal_drivers <- cbind(temporal_drivers, interaction = temporal_drivers[, 1] * temporal_drivers[, 2])
temporal_drivers <- as_data_frame(temporal_drivers) %>%
  mutate(year = 1985:2015)

rawdat <- phi_ensemble_estimate %>% 
  filter(parameter == "Start of season") %>%
  group_by(foresttype) %>%
  mutate(median_deviance = median - mean(median),
         lower_deviance = lower - mean(median),
         upper_deviance = upper - mean(median)) %>%
  ungroup() %>%
  left_join(temporal_drivers, by = "year")

p <- ggplot(response_curves, aes(x = preseason_temp_anomaly, y = pred)) +
  geom_point(data = rawdat, 
             aes(x = preseason_temp_anomaly, y = median_deviance, col = chilling_anomaly),
             shape = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(chilling_anomaly)), alpha = 0.25) +
  geom_line(aes(col = chilling_anomaly, group = chilling_anomaly)) +
  theme_bw() +
  facet_wrap(~foresttype, ncol = 2) +
  labs(x = "Temperature", y = "Predicted variability\nin start of season", fill = "Chilling:", col = "Chilling:") +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "line")) +
  scale_color_gradient2(mid = "darkgrey",
                        breaks = c(-2, 0, 2),
                        guide = guide_legend(direction = "horizontal", 
                                             title.position = "top")) +
  scale_fill_manual(values = c(scales::muted("red"), 
                               "darkgrey", 
                               scales::muted("blue")), guide = "none")

ggsave(paste0("response.pdf"), p, path = "figures/", width = 3.5, height = 2)

## Validation (Figure 5) ------------------------------------------------

mod <- validation_dat %>%
  split(.$foresttype) %>%
  map(., ~ lmodel2::lmodel2(sos_prediction_median ~ unfolding_deviance_mean, data = ., "interval", "interval", 99)) %>%
  map(., ~ .$regression.results) %>%
  bind_rows(.id = "foresttype")

p <- ggplot(validation_dat, aes(x = unfolding_deviance_mean, y = sos_prediction_median)) +
  geom_point(aes(col = foresttype, shape = foresttype)) +
  geom_errorbar(aes(ymin = sos_prediction_lower, ymax = sos_prediction_upper, col = foresttype), width = 0.25, alpha = 0.5) +
  geom_errorbarh(aes(xmin = unfolding_deviance_lower, xmax = unfolding_deviance_upper, col = foresttype), height = 0.25, alpha = 0.5) +
  ylim(-20, 28) +
  geom_abline(intercept = 0, slope = 1, col = "grey", linetype = "dashed") +
  theme_bw() +
  scale_x_continuous(breaks = c(-20, -10, 0, 10, 20), limits = c(-20, 28)) +
  labs(x = "Deviance in leaf-out", 
       y = "Predicted variability\nin start of season", 
       col = NULL, shape = NULL) +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.background = element_blank()) +
  scale_color_manual(values = colors) +
  guides(col = guide_legend(keywidth = 0.15,
                            keyheight = 0.15,
                            default.unit = "inch")) +
  geom_abline(data = filter(mod, Method == "SMA"), aes(intercept = Intercept, slope = Slope, colour = foresttype), alpha = 0.75) +
  facet_wrap(~foresttype)

ggsave(paste0("validation.pdf"), p, path = "figures/", width = 3.5, height = 2)

## Plot ensemble development of rho (Figure S1) -------------------------

p <- ggplot(rho_development, aes(x = permutations, y = median, 
                                 col = foresttype, shape = foresttype)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.3) +
  facet_wrap(~predictor, ncol = 3) +
  #facet_grid(predictor ~ foresttype, scales = "free_y") +
  theme_bw() +
  labs(x = "Number of samples (with n = 10)", 
       y = "Estimated effect on variability\nin start of season",
       col = NULL, shape = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "darkgrey") +
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_blank()) +
  scale_color_manual(values = colors) +
  guides(col = guide_legend(ncol = 1, 
                            keywidth = 0.1,
                            keyheight = 0.1,
                            default.unit = "inch"))

ggsave(paste0("ensemble_development.pdf"), p, path = "figures/", width = 7.5, height = 2.5)

## Posterior predictions (Figure S2) ---------------------------------------------------

p <- simulations_dat %>%
  split(list(.$model, .$foresttype)) %>%
  map2(.y = c(1:10, 1:10), 
       ~ filter(., pixel %in% .y)) %>%
  bind_rows() %>%
  mutate(foresttype = ifelse(foresttype == 1, "Coniferous", "Broadleaved"),
         pixel = paste0("Pixel: ", pixel)) %>%
  ggplot(., aes(x = doy, y = index)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "darkgrey", width = 0) + 
  facet_wrap(foresttype~pixel, scales = "free") + 
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "EVI", col = "Year") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7),
        legend.position = "none")

ggsave(paste0("posterior_simulations.pdf"), p, path = "figures/", width = 7.5, height = 7.5)

