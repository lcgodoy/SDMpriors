functions {
}
data {
  // Observational (Field) Data
  int<lower=1> N_obs;
  array[N_obs] int<lower = 0, upper = 1> y_obs;
  // Experimental (Literature) Data
  int<lower = 1> K;
  matrix[N_obs, K] X;
  vector[K] pri_beta;
}
transformed data {
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower = 0> sigma;
}
transformed parameters {
  vector[N_obs] logit_mu;
  logit_mu = alpha + X * beta;
}
model {
  target += cauchy_lpdf(sigma | 0, 1) -
    cauchy_lccdf(sigma | 0, 1);
  for (k in 1:K)
    target += normal_lpdf(beta[k] | pri_beta[k], sigma);
  y_obs ~ bernoulli_logit(logit_mu);
}
