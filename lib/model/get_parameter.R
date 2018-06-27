
#### Function for extracting parameters from a phenoBayes model

get_parameter <- function (
  mod,
  parameter, # Currebtly implemented: phi, beta, rho
  p_lower = 0.2,
  p_upper = 0.8
) {
  
  if (parameter == "phi") {
  
    iter <- nrow(mod$phi[[1]])
    years <- sort(unique(mod$data$year))
    
    get_phi <- function (x) {
      x <- tidyr::gather(as.data.frame(x), key = year, value = value)
      x$year <- rep(years, each = iter)
      x <- dplyr::summarize(dplyr::group_by(x, year),
                            median = median(value),
                            lower = quantile(value, p_lower),
                            upper = quantile(value, p_upper))
      return(dplyr::ungroup(x))
    }
    
    phi1 <- get_phi(mod$phi[[1]])
    phi2 <- get_phi(mod$phi[[2]])
    phi3 <- get_phi(mod$phi[[3]])
    phi4 <- get_phi(mod$phi[[4]])
  
    phi <- dplyr::bind_rows(list(phi1 = phi1, phi2 = phi2, phi3 = phi3, phi4 = phi4), .id = "parameter")
    
    return(phi)
    
  } else if (parameter == "beta") {
    
    iter <- nrow(mod$beta[[1]])
    pixels <- unique(mod$data$pixel)
    
    get_beta <- function (x) {
      x <- tidyr::gather(as.data.frame(x), key = year, value = value)
      x$pixel <- rep(pixels, each = iter)
      x <- dplyr::summarize(dplyr::group_by(x, pixel),
                            median = median(value),
                            lower = quantile(value, p_lower),
                            upper = quantile(value, p_upper))
      return(dplyr::ungroup(x))
    }
    
    beta1 <- get_beta(mod$beta[[1]])
    beta2 <- get_beta(mod$beta[[2]])
    beta3 <- get_beta(mod$beta[[3]])
    beta4 <- get_beta(mod$beta[[4]])
    
    beta <- dplyr::bind_rows(list(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4), .id = "parameter")
    
    return(beta)
    
  } else if (parameter == "rho") {
    
    iter <- nrow(mod$rho)
    names <- paste0("rho_", 1:ncol(mod$rho))
    
    rho <- tidyr::gather(as.data.frame(mod$rho), key = predictor, value = value)
    rho$predictor <- rep(names, each = iter)
    rho <- dplyr::summarize(dplyr::group_by(rho, predictor),
                          median = median(value),
                          lower = quantile(value, p_lower),
                          upper = quantile(value, p_upper))
    
    return(rho)
    
  }
  
}



