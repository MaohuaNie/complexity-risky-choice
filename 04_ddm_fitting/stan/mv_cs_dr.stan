data {
  int<lower=1> N;                          // number of data items
  int<lower=1> L;                          // number of participants
  array[N] int<lower=1, upper=L> participant;    // participant ID
  array[N] int<lower=-1, upper=1> cho;           // choice: complex choice 1, easy choice -1
  array[N] real<lower=0> rt;                     // response times
  array[N] real skew;                             // skewness parameter
  array[N] real evd;                              // expected value difference
  array[N] real sdd;                              // standard deviation difference
  real<lower=0, upper=1> starting_point;    // starting point for DDM (not estimated)
}

parameters {
  vector[6] mu;                             // population-level means for beta, theta, threshold, ndt, eta, zeta
  vector<lower=0>[6] sigma;                 // population-level standard deviations
  cholesky_factor_corr[6] L_corr;           // Cholesky factor of correlation matrix
  matrix[6, L] z;                           // unscaled participant-level deviations
}

transformed parameters {
  matrix[6, L] participant_params;          // participant-specific parameters (beta, theta, threshold, ndt, eta, zeta)
  vector[N] drift_ll;                       // trial-by-trial drift rate for likelihood
  vector[N] drift_t;                        // trial-by-trial drift rate for predictions
  vector<lower=0>[N] threshold_t;           // trial-by-trial threshold
  vector<lower=0>[N] ndt_t;                 // trial-by-trial non-decision time
  
  // Non-centered parameterization for participant-level parameters
  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] = participant_params[, l] + mu;  // add the population-level means to each participant's parameters
  }
  


  for (n in 1:N) {
    real beta = participant_params[1, participant[n]];
    real theta = log(1 + exp(participant_params[2, participant[n]]));  // softplus transformation for choice consistency
    real threshold = log(1 + exp(participant_params[3, participant[n]]));
    real ndt = log(1 + exp(participant_params[4, participant[n]]));
    real eta = participant_params[5, participant[n]];
    real zeta = participant_params[6, participant[n]];

    drift_t[n] = theta * (evd[n] + beta * sdd[n] + eta * skew[n]  + zeta);
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
  }
}

model {
  // Priors for population-level parameters
  mu[1] ~ normal(0, 3);  // Prior for beta's population mean
  mu[2] ~ normal(-2, 3);  // Prior for theta's population mean
  mu[3] ~ normal(5, 3);  // Prior for threshold's population mean
  mu[4] ~ normal(-1, 1);  // Prior for ndt's population mean
  mu[5] ~ normal(0, 3);  // Prior for eta's population mean
  mu[6] ~ normal(0, 5);  // Prior for zeta's population mean
  
  sigma ~ exponential(1);                             // weakly informative priors for standard deviations
  L_corr ~ lkj_corr_cholesky(1);                   // LKJ prior for correlation matrix

  // Priors for participant-level deviations (non-centered)
  to_vector(z) ~ normal(0, 1);

  // Likelihood using the Wiener distribution for response times
  rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}

generated quantities {
  matrix[6, 6] Omega;                       // correlation matrix
  vector[N] log_lik;                        // log-likelihood for each observation
  
  // Recover the full correlation matrix from the Cholesky factor
  Omega = multiply_lower_tri_self_transpose(L_corr);

  // Calculate log-likelihood for each observation
  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
  }
}

