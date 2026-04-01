functions {
#include utils/conditionals.stan
}
data {
  // Observational (Field) Data
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  array[N_obs] int<lower = 0, upper = 1> y_obs;
  // Experimental (Literature) Data
  int<lower = 1> N_exp;
  array[N_exp] real x_exp;
  vector[N_exp] f_exp;  // We pass the fixed logit values here (e.g., -5.0)
  real<lower = 0, upper = 1> p_rho;
  real<lower = 0> rho_thresh;
}
transformed data {
  real beta_rho;
  beta_rho = - log(p_rho) / rho_thresh;
}
parameters {
  real<lower = 0> rho;
  real<lower = 0> alpha;
  real<lower = 0> kappa;
  vector[N_obs] eta;  // For non-centered parameterization
}
transformed parameters {
  vector[N_obs] f_obs;
  {
    matrix[N_obs, N_obs] K_obs_cond;
    matrix[N_obs, N_obs] L_K_cond;
    // Build the sub-matrices
    matrix[N_exp, N_exp] K_exp_exp =
      add_diag(gp_matern52_cov(x_exp, kappa * alpha, rho),
               1e-10);
    matrix[N_obs, N_obs] K_obs_obs =
      add_diag(gp_matern52_cov(x_obs, alpha, rho),
               1e-10);
    matrix[N_exp, N_obs] K_exp_obs =
      gp_matern52_cov(x_exp, x_obs, sqrt(kappa) * alpha, rho);
    // conditional parameters for GP
    tuple(vector[N_obs], matrix[N_obs, N_obs]) musigma_cond;
    musigma_cond = gp_cond_pars(f_exp, K_exp_exp, K_obs_obs, K_exp_obs);
    L_K_cond = cholesky_decompose(musigma_cond.2);
    // Construct the final latent GP for the observations
    f_obs = musigma_cond.1 + L_K_cond * eta;
  }
}
model {
  // Priors
  // rho ~ inv_gamma(2, 1);
  rho ~ exponential(beta_rho);
  alpha ~ std_normal();
  kappa ~ std_normal();
  eta ~ std_normal();
  // Likelihood (Only evaluated on the real field data!)
  y_obs ~ bernoulli_logit(f_obs);
}
