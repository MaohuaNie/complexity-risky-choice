functions {
  real pweight(real prop, real gamma) {
    return exp(-pow(-log(prop), gamma));
  }
}


data {
  int<lower=1> N;                                   // number of data items
  int<lower=1> L;                                   // number of participants
  array[N] int<lower=1, upper=L> participant;             // level (participant)
  array[N] int<lower=0, upper=1> accuracy_flipped;        // flipped accuracy (1, 0)
  array[N] int<lower=-1, upper=1> cho;                    // accuracy (-1, 1)
  array[N] real<lower=0> rt;                              // response times
  array[N, 2] real o_complex;                             // outcomes for option complex 
  array[N, 2] real o_simple;                              // outcomes for option simple 
  array[N, 2] real p_complex;                             // probabilities for option complex
  array[N, 2] real p_simple;                              // probabilities for option simple

}

parameters {
  vector[7] mu;                                     // population-level means for beta, theta, threshold, ndt, gamma, sp, zeta
  vector<lower=0>[7] sigma;                         // population-level standard deviations
  cholesky_factor_corr[7] L_corr;                   // Cholesky factor of correlation matrix
  matrix[7, L] z;                                   // unscaled participant-level deviations
}

transformed parameters {
  matrix[7, L] participant_params;                  // participant-specific parameters (beta, theta, threshold, ndt, gamma, sp, zeta)
  vector[L] beta_t;                                 // transformed beta for each participant
  vector[L] theta_t;                                // transformed theta for each participant
  vector[L] threshold_t_p;                          // transformed threshold for each participant
  vector[L] ndt_t_p;                                // transformed ndt for each participant
  vector[L] gamma_t;                                // transformed gamma for each participant
  vector<lower=0, upper=1>[L] rel_sp_t_p;          // transformed starting point for each participant
  vector[L] zeta_t;                                 // transformed zeta for each participant
  vector[N] drift_ll;                               // trial-by-trial drift rate for likelihood
  vector[N] drift_t;                                // trial-by-trial drift rate for predictions
  vector<lower=0>[N] threshold_t;                   // trial-by-trial threshold
  vector<lower=0>[N] ndt_t;                         // trial-by-trial non-decision time
  array[N] real u_complex;                                       // utility for option complex
  array[N] real u_simple;                                       // utility for option simple
  vector<lower=0, upper=1>[N] rel_sp_ll;    // trial-by-trial relative starting point for likelihood 
  vector<lower=0, upper=1>[N] rel_sp_t;     // trial-by-trial relative starting point

  // Non-centered parameterization for participant-level parameters
  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] = participant_params[, l] + mu;  // add the population-level means to each participant's parameters
  }

  // Transform parameters once per participant (not per trial)
  for (l in 1:L) {
    beta_t[l] = log1p_exp(participant_params[1, l]);
    theta_t[l] = log1p_exp(participant_params[2, l]);
    threshold_t_p[l] = log1p_exp(participant_params[3, l]);
    ndt_t_p[l] = log1p_exp(participant_params[4, l]);
    gamma_t[l] = log1p_exp(participant_params[5, l]);
    rel_sp_t_p[l] = Phi(participant_params[6, l]);
    zeta_t[l] = participant_params[7, l];
  }
  
  // Loop over trials to compute utilities and drift rates
  for (n in 1:N) {
    
    real beta = beta_t[participant[n]];
    real theta = theta_t[participant[n]];
    real threshold = threshold_t_p[participant[n]];
    real ndt = ndt_t_p[participant[n]];
    real gamma = gamma_t[participant[n]];
    real zeta = zeta_t[participant[n]];


    // Compute utilities for options A and B
    u_complex[n] = pweight(p_complex[n,1], gamma) * pow(o_complex[n,1], beta) + (1 - pweight(p_complex[n,1], gamma)) * pow(o_complex[n,2], beta);
    u_simple[n] = pweight(p_simple[n,1], gamma) * pow(o_simple[n,1], beta) + (1 - pweight(p_simple[n,1], gamma)) * pow(o_simple[n,2], beta);

    // Compute trial-by-trial drift rates
    drift_t[n] = theta * (pow(u_complex[n], 1 / beta) - pow(u_simple[n], 1 / beta) + zeta);
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
    
    
    rel_sp_t[n] = rel_sp_t_p[participant[n]];
    rel_sp_ll[n] = accuracy_flipped[n] + cho[n] * rel_sp_t[n];
  }
}

model {
  // Priors for population-level parameters
  mu[1] ~ normal(0.5414, 3);  // Prior for beta's population mean
  mu[2] ~ normal(-1, 3);  // Prior for theta's population mean
  mu[3] ~ normal(5, 3);  // Prior for threshold's population mean
  mu[4] ~ normal(-1, 1);  // Prior for ndt's population mean
  mu[5] ~ normal(0.5414, 3);  // Prior for gamma's population mean
  mu[6] ~ normal(0, 3);  // Prior for sp's population mean
  mu[7] ~ normal(0, 3);  // Prior for zeta's population mean

  
  sigma ~ exponential(1);                             // weakly informative priors for standard deviations
  L_corr ~ lkj_corr_cholesky(1);                   // LKJ prior for correlation matrix

  // Priors for participant-level deviations (non-centered)
  to_vector(z) ~ normal(0, 1);

  // Likelihood using the Wiener distribution for response times
  rt ~ wiener(threshold_t, ndt_t, rel_sp_ll, drift_ll);
}

generated quantities {
  matrix[7, 7] Omega;                               // correlation matrix
  vector[N] log_lik;                                // log-likelihood for each observation


  // Recover the full correlation matrix from the Cholesky factor
  Omega = multiply_lower_tri_self_transpose(L_corr);

  // Calculate log-likelihood for each observation
  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], rel_sp_ll[n], drift_ll[n]);
  }
}

