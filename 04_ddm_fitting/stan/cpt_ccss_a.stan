functions {
  real pweight(real prop, real gamma) {
    return exp(-pow(-log(prop), gamma));
  }
}

data {
  int<lower=1> N;                                   // number of data items
  int<lower=1> L;                                   // number of participants
  array[N] int<lower=1, upper=L> participant;             // level (participant)
  array[N] int<lower=-1, upper=1> cho;                    // accuracy (1, -1)
  array[N] real<lower=0> rt;                              // response times
  array[N, 2] real o_risky;                                    // outcomes for option risky 
  array[N, 2] real o_safe;                                    // outcomes for option safe 
  array[N, 2] real p_risky;                                    // probabilities for option risky
  array[N, 2] real p_safe;                                    // probabilities for option safe
  real<lower=0, upper=1> starting_point;            // starting point for diffusion model, fixed value
  array[N] int<lower=-1, upper=1> con;                    // condition index (-1 = simple vs. simple, 1 = complex vs. complex)
}

parameters {
  vector[6] mu;                                     // population-level means for beta, theta, threshold, ndt, gamma, delta_threshold
  vector<lower=0>[6] sigma;                         // population-level standard deviations
  cholesky_factor_corr[6] L_corr;                   // Cholesky factor of correlation matrix
  matrix[6, L] z;                                   // unscaled participant-level deviations
}

transformed parameters {
  matrix[6, L] participant_params;                  // participant-specific parameters (beta, theta, threshold, ndt, gamma, delta_threshold)
  vector[L] beta_t;                                 // transformed beta for each participant
  vector[L] theta_t;                                // transformed theta for each participant
  vector[L] threshold_base_t;                       // transformed base threshold for each participant
  vector[L] ndt_t_p;                                // transformed ndt for each participant
  vector[L] gamma_t;                                // transformed gamma for each participant
  vector[L] delta_threshold_t;                      // delta_threshold for each participant (not transformed)
  vector[N] drift_ll;                               // trial-by-trial drift rate for likelihood
  vector[N] drift_t;                                // trial-by-trial drift rate for predictions
  vector<lower=0>[N] threshold_t;                   // trial-by-trial threshold
  vector<lower=0>[N] ndt_t;                         // trial-by-trial non-decision time
  array[N] real u_risky;                                       // utility for option risky
  array[N] real u_safe;                                       // utility for option safe

  // Non-centered parameterization for participant-level parameters
  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] = participant_params[, l] + mu;  // add the population-level means to each participant's parameters
  }

  // Transform parameters once per participant (not per trial)
  for (l in 1:L) {
    beta_t[l] = log1p_exp(participant_params[1, l]);
    theta_t[l] = log1p_exp(participant_params[2, l]);
    threshold_base_t[l] = participant_params[3, l];  // Keep untransformed for now
    ndt_t_p[l] = log1p_exp(participant_params[4, l]);
    gamma_t[l] = log1p_exp(participant_params[5, l]);
    delta_threshold_t[l] = participant_params[6, l];
  }

  // Loop over trials to compute utilities and drift rates
  for (n in 1:N) {
    
    real beta = beta_t[participant[n]];
    real theta = theta_t[participant[n]];
    real threshold = log1p_exp(threshold_base_t[participant[n]] + delta_threshold_t[participant[n]] * con[n]);
    real ndt = ndt_t_p[participant[n]];
    real gamma = gamma_t[participant[n]];
    
  
    // Compute utilities for options A and B
    u_risky[n] = pweight(p_risky[n,1], gamma) * pow(o_risky[n,1], beta) + (1 - pweight(p_risky[n,1], gamma)) * pow(o_risky[n,2], beta);
    u_safe[n] = pweight(p_safe[n,1], gamma) * pow(o_safe[n,1], beta) + (1 - pweight(p_safe[n,1], gamma)) * pow(o_safe[n,2], beta);

    // Compute trial-by-trial drift rates
    drift_t[n] = theta * (pow(u_risky[n], 1 / beta) - pow(u_safe[n], 1 / beta));
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
  }
}

model {
  // Priors for population-level parameters
  mu[1] ~ normal(0.5414, 3);  // Prior for beta's population mean
  mu[2] ~ normal(-1, 3);  // Prior for theta's population mean
  mu[3] ~ normal(5, 3);  // Prior for threshold's population mean
  mu[4] ~ normal(-1, 1);  // Prior for ndt's population mean
  mu[5] ~ normal(0.5414, 3);  // Prior for gamma's population mean
  mu[6] ~ normal(0, 3);  // Prior for delta_threshold's population mean
  
  sigma ~ exponential(1);                             // weakly informative priors for standard deviations
  L_corr ~ lkj_corr_cholesky(1);                   // LKJ prior for correlation matrix

  // Priors for participant-level deviations (non-centered)
  to_vector(z) ~ normal(0, 1);

  // Likelihood using the Wiener distribution for response times
  rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}

generated quantities {
  matrix[6, 6] Omega;                               // correlation matrix
  vector[N] log_lik;                                // log-likelihood for each observation

  // Recover the full correlation matrix from the Cholesky factor
  Omega = multiply_lower_tri_self_transpose(L_corr);

  // Calculate log-likelihood for each observation
  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
  }
}

