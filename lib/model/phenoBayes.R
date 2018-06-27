
#### Run a phenoBayes model

phenoBayes <- function(pixel,
                       index, 
                       doy, 
                       year,
                       parameter = "start",
                       temporal_drivers,
                       cut_winter = 80,
                       inital_iterations = 25,
                       n_cores = 4, 
                       ...
                       ) {
  
  #### Combine data #### ----------------------------------------------------------------------------------------------------------------------------
  
  dat <- data.frame(pixel,
                    index,
                    doy,
                    year)
  
  #### Estimate initial parameters #### ----------------------------------------------------------------------------------------------------------------------------
  
  ### Spring/autumn split
  
  transition <- c()
  pixel_id <- c()
  
  for (p in unique(dat$pixel)) {
    
    d <- subset(dat, dat$pixel == p & dat$doy > cut_winter)
    d <- d[order(d$doy),] 
    
    slope <- c()
    ddoy <- c()
    for(i in 21:(nrow(d) - 21)) {
      dd <- d[(i - 20):(i + 20), ]
      dd$t <- 1:nrow(dd)
      slope <- c(slope, coef(lm(dd$index ~ dd$t))[2])
      ddoy <- c(ddoy, d[i, "doy"])
    }
    
    transition <- c(transition, unlist(ddoy)[which(slope < 0 & 
                                                   dplyr::lead(slope) < 0 & 
                                                   dplyr::lead(slope, 2) < 0)[1]])
    pixel_id <- c(pixel_id, p)
  }
  
  transition <- data.frame(pixel = pixel_id,
                           transition = transition)
  
  dat <- merge(dat, transition, by = "pixel")
  
  dat <- subset(dat, dat$doy <= dat$transition & dat$doy > cut_winter)
  
  ### Initial parameters ----------------------------------------------------------------------------------------------------------------------------
  
  initial_parameters <- vector("list", inital_iterations)
  
  k <- 0
  
  while (k < inital_iterations) {
    dat_sample <- subset(dat, dat$pixel == sample(unique(dat$pixel), 1))
    nls_fit <- tryCatch(nls(index ~ b1 + (b2 / (1 + exp(- b3 * (doy - b4)))),
                   start = list(b1 = min(index), 
                                b2 = max(index), 
                                b3 = 0.15, 
                                b4 = round(mean(dat_sample[which(dat_sample$index > median(dat_sample$index)), "doy"]), 0)),
                   data = dat_sample), error = function(e) return(NA))
    if (class(nls_fit) == "nls") {
      k <- k + 1
      initial_parameters[[k]] <- as.data.frame(t(coef(nls_fit)))
      # print(paste("Find parameter, iteration", k, "out of 10."))
    }
  }
  
  initial_parameters <- apply(do.call("rbind", initial_parameters), 2, mean)

  #### Run pheoBayes model #### ----------------------------------------------------------------------------------------------------------------------------
  
  # Set priors
  mu_beta <- mu_beta <- matrix(initial_parameters[c(1:4)],
                               ncol = max(dat$pixel),
                               nrow = 4)
  
  # Create data list for stan
  dat$year_index <- dat$year - min(dat$year) + 1
  datalist <- list(N = nrow(dat),
                   NY = max(dat$year_index),
                   NP = max(dat$pixel),
                   doy = dat$doy,
                   vi = dat$index,
                   year = dat$year_index,
                   pixel = dat$pixel,
                   mu_beta = mu_beta,
                   scale_sigma = c(0.1, 0.1, 0.1, 1),
                   N_spatial = ncol(spatial_drivers),
                   spatial_drivers = spatial_drivers,
                   N_temporal = ncol(temporal_drivers),
                   temporal_drivers = temporal_drivers)
  
  # Run Stan model
  options(mc.cores = n_cores)
  rstan::rstan_options(auto_write = TRUE)
  
  library(rstan)
  
  model <- "lib/model/spring_temporal_drivers_sos.stan"
  
  fit <- rstan::stan(file = model, 
                     data = datalist, 
                     init = 0,
                     iter = 4000,
                     chains = 4,
                     ...)
  
  #### Extract model diagnostics #### ----------------------------------------------------------------------------------------------------------------------------
  
  # Rhat
  rhat <- bayesplot::rhat(fit)
  print(paste0(round(mean(rhat < 1.05, na.rm = TRUE) * 100, 2), "% of the Rhat values were below 1.05. Percentages >95% indicate good convergence!"))
  
  neff <- bayesplot::neff_ratio(fit)
  print(paste0(round(mean(neff > 0.5, na.rm = TRUE) * 100, 2), "% of the Neff values are 'good'."))

  # Posterior predictions and p-values
  draws <- as.matrix(fit)
  sim <- draws[, grep(glob2rx("sim*"), colnames(draws))]
    
  #### Extract results #### ----------------------------------------------------------------------------------------------------------------------------
  
  if (!exists("draws")) draws <- as.matrix(fit)
  
  phi1 <- draws[, grep(glob2rx("phi[1*"), colnames(draws))] + initial_parameters[1]
  phi2 <- draws[, grep(glob2rx("phi[2*"), colnames(draws))] + initial_parameters[2]
  phi3 <- draws[, grep(glob2rx("phi[3*"), colnames(draws))] + initial_parameters[3]
  phi4 <- draws[, grep(glob2rx("phi[4*"), colnames(draws))] + initial_parameters[4]
  
  beta1 <- draws[, grep(glob2rx("beta[1,*"), colnames(draws))]
  beta2 <- draws[, grep(glob2rx("beta[2,*"), colnames(draws))]
  beta3 <- draws[, grep(glob2rx("beta[3,*"), colnames(draws))]
  beta4 <- draws[, grep(glob2rx("beta[4,*"), colnames(draws))]
  
  rho <- draws[, grep(glob2rx("rho[*"), colnames(draws))]
  lambda <- draws[, grep(glob2rx("lambda[*"), colnames(draws))]
  
  cor_beta <- draws[, grep(glob2rx("cor_beta[*"), colnames(draws))]
  cor_phi <- draws[, grep(glob2rx("cor_phi[*"), colnames(draws))]
  
  #### Return results ##### ----------------------------------------------------------------------------------------------------------------------------
  
  return(list(data = dat, 
              phi = list(phi1, phi2, phi3, phi4),
              beta = list(beta1, beta2, beta3, beta4),
              rho = rho,
              lambda = lambda,
              cor_beta = cor_beta,
              cor_phi = cor_phi,
              rhat = rhat,
              posterior_simulations = sim))

}