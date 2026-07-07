data {
  int<lower=1> N;
  int<lower=1> L;
  array[N] int<lower=1, upper=L> participant;
  array[N] int<lower=-1, upper=1> cho;
  array[N] real<lower=0> rt;
  array[N] real skew;
  array[N] real evd;
  array[N] real sdd;
  real<lower=0, upper=1> starting_point;
  array[N] int<lower=-1, upper=1> con;
}

parameters {
  vector[7] mu; // (beta, theta, threshold, ndt, eta, delta_theta, delta_eta)
  vector<lower=0>[7] sigma;
  cholesky_factor_corr[7] L_corr;
  matrix[7, L] z;
}

transformed parameters {
  matrix[7, L] participant_params;  
  vector[N] drift_ll;
  vector[N] drift_t;
  vector<lower=0>[N] threshold_t;
  vector<lower=0>[N] ndt_t;


  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] += mu;
  }

  for (n in 1:N) {
    real beta = participant_params[1, participant[n]];
    real theta = log(1 + exp(participant_params[2, participant[n]] + participant_params[6, participant[n]]*con[n]));
    real threshold = log(1 + exp(participant_params[3, participant[n]]));
    real ndt = log(1 + exp(participant_params[4, participant[n]]));
    real eta = participant_params[5, participant[n]] + participant_params[7, participant[n]]*con[n];
    
    
    drift_t[n] = theta * (evd[n] + beta * sdd[n] + eta * skew[n]);
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
  }
}

model {
  mu[1] ~ normal(0, 3);    
  mu[2] ~ normal(-2, 3);   
  mu[3] ~ normal(5, 3);    
  mu[4] ~ normal(-1, 1);  
  mu[5] ~ normal(0, 3);    
  mu[6] ~ normal(0, 3);    // delta_theta
  mu[7] ~ normal(0, 3);    // delta_eta

  sigma ~ exponential(1);
  L_corr ~ lkj_corr_cholesky(1);

  to_vector(z) ~ normal(0, 1);

  rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}

generated quantities {
  matrix[7,7] Omega = multiply_lower_tri_self_transpose(L_corr);
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
  }
}


