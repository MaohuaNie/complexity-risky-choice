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
  array[N, 9] real o_risky;
  array[N, 9] real o_safe;
  array[N, 9] real p_risky;
  array[N, 9] real p_safe;
  real<lower=0, upper=1> starting_point;
  array[N] int<lower=-1, upper=1> con;
}

parameters {
  vector[6] mu;
  vector<lower=0>[6] sigma;
  cholesky_factor_corr[6] L_corr;
  matrix[6, L] z;
}

transformed parameters {
  matrix[6, L] participant_params;
  vector[L] beta_base_t;
  vector[L] theta_t;
  vector[L] threshold_t_p;
  vector[L] gamma_t;
  vector[L] ndt_t_p;
  vector[L] delta_beta_t;
  vector[N] drift_ll;
  vector[N] drift_t;
  vector<lower=0>[N] threshold_t;
  vector<lower=0>[N] ndt_t;
  array[N] real u_risky;
  array[N] real u_safe;

  participant_params = diag_pre_multiply(sigma, L_corr) * z;
  for (l in 1:L) {
    participant_params[, l] = participant_params[, l] + mu;
  }

  for (l in 1:L) {
    beta_base_t[l] = participant_params[1, l];
    theta_t[l] = log1p_exp(participant_params[2, l]);
    threshold_t_p[l] = log1p_exp(participant_params[3, l]);
    gamma_t[l] = log1p_exp(participant_params[5, l]);
    ndt_t_p[l] = log1p_exp(participant_params[4, l]);
    delta_beta_t[l] = participant_params[6, l];
  }

  for (n in 1:N) {
    real beta = log1p_exp(beta_base_t[participant[n]] + delta_beta_t[participant[n]] * con[n]);
    real theta = theta_t[participant[n]];
    real threshold = threshold_t_p[participant[n]];
    real ndt = ndt_t_p[participant[n]];
    real gamma = gamma_t[participant[n]];

    if (con[n] == 1) {
      // 7-outcome rank-dependent CPT (cols 3..9)
      real w1r = pweight(p_risky[n, 3], gamma);
      real w2r = pweight(p_risky[n, 3] + p_risky[n, 4], gamma) - w1r;
      real w3r = pweight(p_risky[n, 3] + p_risky[n, 4] + p_risky[n, 5], gamma) - w2r - w1r;
      real w4r = pweight(p_risky[n, 3] + p_risky[n, 4] + p_risky[n, 5] + p_risky[n, 6], gamma) - w3r - w2r - w1r;
      real w5r = pweight(p_risky[n, 3] + p_risky[n, 4] + p_risky[n, 5] + p_risky[n, 6] + p_risky[n, 7], gamma) - w4r - w3r - w2r - w1r;
      real w6r = pweight(p_risky[n, 3] + p_risky[n, 4] + p_risky[n, 5] + p_risky[n, 6] + p_risky[n, 7] + p_risky[n, 8], gamma) - w5r - w4r - w3r - w2r - w1r;
      real w7r = 1 - (w6r + w5r + w4r + w3r + w2r + w1r);

      u_risky[n] = w1r * pow(o_risky[n, 3], beta) + w2r * pow(o_risky[n, 4], beta) +
                   w3r * pow(o_risky[n, 5], beta) + w4r * pow(o_risky[n, 6], beta) +
                   w5r * pow(o_risky[n, 7], beta) + w6r * pow(o_risky[n, 8], beta) +
                   w7r * pow(o_risky[n, 9], beta);

      real w1s = pweight(p_safe[n, 3], gamma);
      real w2s = pweight(p_safe[n, 3] + p_safe[n, 4], gamma) - w1s;
      real w3s = pweight(p_safe[n, 3] + p_safe[n, 4] + p_safe[n, 5], gamma) - w2s - w1s;
      real w4s = pweight(p_safe[n, 3] + p_safe[n, 4] + p_safe[n, 5] + p_safe[n, 6], gamma) - w3s - w2s - w1s;
      real w5s = pweight(p_safe[n, 3] + p_safe[n, 4] + p_safe[n, 5] + p_safe[n, 6] + p_safe[n, 7], gamma) - w4s - w3s - w2s - w1s;
      real w6s = pweight(p_safe[n, 3] + p_safe[n, 4] + p_safe[n, 5] + p_safe[n, 6] + p_safe[n, 7] + p_safe[n, 8], gamma) - w5s - w4s - w3s - w2s - w1s;
      real w7s = 1 - (w6s + w5s + w4s + w3s + w2s + w1s);

      u_safe[n] = w1s * pow(o_safe[n, 3], beta) + w2s * pow(o_safe[n, 4], beta) +
                  w3s * pow(o_safe[n, 5], beta) + w4s * pow(o_safe[n, 6], beta) +
                  w5s * pow(o_safe[n, 7], beta) + w6s * pow(o_safe[n, 8], beta) +
                  w7s * pow(o_safe[n, 9], beta);
    } else {
      // 2-outcome CPT (cols 1..2)
      u_risky[n] = pweight(p_risky[n, 1], gamma) * pow(o_risky[n, 1], beta) +
                   (1 - pweight(p_risky[n, 1], gamma)) * pow(o_risky[n, 2], beta);
      u_safe[n] = pweight(p_safe[n, 1], gamma) * pow(o_safe[n, 1], beta) +
                  (1 - pweight(p_safe[n, 1], gamma)) * pow(o_safe[n, 2], beta);
    }

    drift_t[n] = theta * (pow(u_risky[n], 1 / beta) - pow(u_safe[n], 1 / beta));
    drift_ll[n] = drift_t[n] * cho[n];
    threshold_t[n] = threshold;
    ndt_t[n] = ndt;
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

  rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}

generated quantities {
  matrix[6, 6] Omega;
  vector[N] log_lik;

  Omega = multiply_lower_tri_self_transpose(L_corr);

  for (n in 1:N) {
    log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
  }
}
