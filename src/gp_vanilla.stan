functions {
}
data {
  // Observational (Field) Data
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  array[N_obs] int<lower = 0, upper = 1> y_obs;
}
parameters {
  real<lower = 0> rho;
  real<lower = 0> alpha;
  vector[N_obs] eta;  // For non-centered parameterization
}
transformed parameters {
  vector[N_obs] f_obs;
  {
    matrix[N_obs, N_obs] K =
      add_diag(gp_matern52_cov(x_obs, alpha, rho),
               1e-10);
    matrix[N_obs, N_obs] L = cholesky_decompose(K);
    f_obs = L * eta;
  }
}
model {
  // Priors
  rho ~ inv_gamma(2, 1);
  alpha ~ std_normal();
  eta ~ std_normal();
  // Likelihood (Only evaluated on the real field data!)
  y_obs ~ bernoulli_logit(f_obs);
}
