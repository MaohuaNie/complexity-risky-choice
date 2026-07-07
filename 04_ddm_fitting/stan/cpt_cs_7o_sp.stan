functions {
  real pweight(real prop, real gamma) {
    return exp(-pow(-log(prop), gamma));
  }
}

data {
  int<lower=1> N;
  int<lower=1> L;
  array[N] int<lower=1, upper=L> participant;
  array[N] int<lower=-1, upper=1> cho;
  array[N] real<lower=0> rt;
  array[N, 7] real o_complex;
  array[N, 2] real o_simple;
  array[N, 7] real p_complex;
  array[N, 2] real p_simple;
  real<lower=0, upper=1> starting_point;
  array[N] int<lower=0, upper=1> accuracy_flipped;
}

parameters {
  vector[6] mu;
  vector<lower=0>[6] sigma;
  cholesky_factor_corr[6] L_corr;
  matrix[6, L] z;
}

transformed parameters {
  matrix[6, L] participant_params;
  vector[L] beta_t;
  vector[L] theta_t;
  vector[L] threshold_t_p;
  vector[L] ndt_t_p;
  vector[L] gamma_t;
  vector<lower=0, upper=1>[L] rel_sp_t_p;
  vector[N] drift_ll;
  vector[N] drift_t;
  vector<lower=0>[N] threshold_t;
  vector<lower=0>[N] ndt_t;
  array[N] real u_complex;
  array[N] real u_simple;
  vector<lower=0, upper=1>[N] rel_sp_ll;
  vector<lower=0, upper=1>[N] rel_sp_t;

  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] = participant_params[, l] + mu;
  }

  for (l in 1:L) {
    beta_t[l] = log1p_exp(participant_params[1, l]);
    theta_t[l] = log1p_exp(participant_params[2, l]);
    threshold_t_p[l] = log1p_exp(participant_params[3, l]);
    ndt_t_p[l] = log1p_exp(participant_params[4, l]);
    gamma_t[l] = log1p_exp(participant_params[5, l]);
    rel_sp_t_p[l] = Phi(participant_params[6, l]);
  }

  for (n in 1:N) {
    real beta = beta_t[participant[n]];
    real theta = theta_t[participant[n]];
    real threshold = threshold_t_p[participant[n]];
    real ndt = ndt_t_p[participant[n]];
    real gamma = gamma_t[participant[n]];

    // 7-outcome rank-dependent CPT for complex option
    real w1c = pweight(p_complex[n, 1], gamma);
    real w2c = pweight(p_complex[n, 1] + p_complex[n, 2], gamma) - w1c;
    real w3c = pweight(p_complex[n, 1] + p_complex[n, 2] + p_complex[n, 3], gamma) - w2c - w1c;
    real w4c = pweight(p_complex[n, 1] + p_complex[n, 2] + p_complex[n, 3] + p_complex[n, 4], gamma) - w3c - w2c - w1c;
    real w5c = pweight(p_complex[n, 1] + p_complex[n, 2] + p_complex[n, 3] + p_complex[n, 4] + p_complex[n, 5], gamma) - w4c - w3c - w2c - w1c;
    real w6c = pweight(p_complex[n, 1] + p_complex[n, 2] + p_complex[n, 3] + p_complex[n, 4] + p_complex[n, 5] + p_complex[n, 6], gamma) - w5c - w4c - w3c - w2c - w1c;
    real w7c = 1 - (w6c + w5c + w4c + w3c + w2c + w1c);

    u_complex[n] = w1c * pow(o_complex[n, 1], beta) + w2c * pow(o_complex[n, 2], beta) +
                   w3c * pow(o_complex[n, 3], beta) + w4c * pow(o_complex[n, 4], beta) +
                   w5c * pow(o_complex[n, 5], beta) + w6c * pow(o_complex[n, 6], beta) +
                   w7c * pow(o_complex[n, 7], beta);

    // 2-outcome CPT for simple option
    u_simple[n] = pweight(p_simple[n, 1], gamma) * pow(o_simple[n, 1], beta) +
                  (1 - pweight(p_simple[n, 1], gamma)) * pow(o_simple[n, 2], beta);

    drift_t[n] = theta * (pow(u_complex[n], 1 / beta) - pow(u_simple[n], 1 / beta));
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
    rel_sp_t[n] = rel_sp_t_p[participant[n]];
    rel_sp_ll[n] = accuracy_flipped[n] + cho[n] * rel_sp_t[n];
  }
}

model {
  mu[1] ~ normal(0.5414, 3);
  mu[2] ~ normal(-1, 3);
  mu[3] ~ normal(5, 3);
  mu[4] ~ normal(-1, 1);
  mu[5] ~ normal(0.5414, 3);
  mu[6] ~ normal(0, 3);

  sigma ~ exponential(1);
  L_corr ~ lkj_corr_cholesky(1);
  to_vector(z) ~ normal(0, 1);

  rt ~ wiener(threshold_t, ndt_t, rel_sp_ll, drift_ll);
}

generated quantities {
  matrix[6, 6] Omega;
  vector[N] log_lik;

  Omega = multiply_lower_tri_self_transpose(L_corr);

  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], rel_sp_ll[n], drift_ll[n]);
  }
}
