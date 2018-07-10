data {
  int<lower=0> N; // Number of observations
  int<lower=0> NY; // Number of years
  int<lower=0> NP; // Number of pixels
  real doy[N]; // Day of year
  real vi[N]; // Vegetation index
  int<lower=0, upper=NY> year[N]; // Observation year (consecutive number 1, ..., NY)
  int<lower=0, upper=NP> pixel[N]; // Observation pixel (consecutive number 1, ..., NP)
  matrix[4, NP] mu_beta; // Centering of beta
  real scale_sigma[4]; // Scaling of sigma phi/beta
  int<lower=0> N_temporal; // Number of temporal drivers
  matrix[NY, N_temporal] temporal_drivers; // Model matrix of temporal drivers
}

parameters {
  
  // Matrix for beta (centered)
  matrix[4, NP] z_beta;
  
  // Correlation matrix for beta with Cholesky factorization
  cholesky_factor_corr[4] cor_beta_L;
  
  // Scale parameters for beta
  vector<lower=0>[4] sigma_beta_raw;
  
  // Vector for phi (centered)
  matrix[4, NY] z_phi;
  
  // Correlation matrix for phi with Cholesky factorization
  cholesky_factor_corr[4] cor_phi_L;

  // Scale parameters for phi
  vector<lower=0>[4] sigma_phi_raw;
  
  // Temporal model parameters
  vector[N_temporal] rho;
  
  // Prediction error scale
  real<lower=0> sigma;
  
}

transformed parameters {
  
  // Matrix for beta (non-centered)
  matrix[4, NP] beta;
  
  // Correlation matrix for beta
  matrix[4, 4] cor_beta;
  
  // Matrix for phi (non-centered)
  matrix[4, NY] phi;
  
  // Correlation matrix for phi
  matrix[4, 4] cor_phi;
  
  // Matrix for temporal model
  matrix[4, NY] temporal_model;
  
  // Scale parameters for beta
  vector<lower=0>[4] sigma_beta;
  
  // Scale parameters for phi
  vector<lower=0>[4] sigma_phi;
  
  // Phi mean estimate (temporal drivers)
  temporal_model = rep_matrix(0, 4, NY);
  temporal_model[4] = (temporal_drivers * rho)';
  
  // Re-scale sigma phi/beta
  for (i in 1:4)
    sigma_beta[i] = sigma_beta_raw[i] * scale_sigma[i];
  
  for (i in 1:4)
    sigma_phi[i] = sigma_phi_raw[i] * scale_sigma[i];
  
  // Cholesky factorization of beta (incl. centering of beta)
  beta = diag_pre_multiply(sigma_beta, cor_beta_L) * z_beta + mu_beta;
  
  // Backtransform Cholesky correlation matrix (LL') of beta 
  cor_beta = tcrossprod(cor_beta_L);
  
  // Cholesky factorization of phi (incl. centering of phi)
  phi = diag_pre_multiply(sigma_phi, cor_phi_L) * z_phi + temporal_model;
  
  // Backtransform Cholesky correlation matrix (LL') of beta 
  cor_phi = tcrossprod(cor_phi_L);
  
}

model {
  
  vector[N] y;
  
  // Prior for beta and phi
  to_vector(z_beta) ~ normal(0, 1);
  to_vector(z_phi) ~ normal(0, 1);
  
  // Prior for correlation matrix of beta and phi
  cor_beta_L ~ lkj_corr_cholesky(2.0);
  cor_phi_L ~ lkj_corr_cholesky(2.0);
  
  // Prior for scales of beta and phi
  sigma_beta_raw ~ cauchy(0, 1);
  sigma_phi_raw ~ cauchy(0, 1);
  
  // Priors for temporal and spatial drivers
  rho ~ cauchy(0, 1);
  
  // Prior prediction error scale
  sigma ~ cauchy(0, 1);
  
  // Likelihood  
  for(i in 1:N)
    y[i] = (beta[1, pixel[i]] + phi[1, year[i]]) + ((beta[2, pixel[i]] + phi[2, year[i]]) / (1 + exp(-(beta[3, pixel[i]] + phi[3, year[i]]) * (doy[i] - (beta[4, pixel[i]] + phi[4, year[i]])))));
  
  vi ~ normal(y, sigma);
  
}

generated quantities {
  
  //vector[N] y;
  //vector[N] log_lik;
  vector[N] sim;

  // Log-likelihood
  //for (n in 1:N) 
  //  y[n] = (beta[1, pixel[n]]) + ((beta[2, pixel[n]]) / (1 + exp(-(beta[3, pixel[n]]) * (doy[n] - (beta[4, pixel[n]] + phi[year[n]])))));
  //for (n in 1:N) log_lik[n] = normal_lpdf(vi[n] | y[n], sigma);
  
  // Posterior simulations
  for (i in 1:N) 
   sim[i] = normal_rng((beta[1, pixel[i]] + phi[1, year[i]]) + ((beta[2, pixel[i]] + phi[2, year[i]]) / (1 + exp(-(beta[3, pixel[i]] + phi[3, year[i]]) * (doy[i] - (beta[4, pixel[i]] + phi[4, year[i]]))))), sigma);
  
}
