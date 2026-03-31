functions {
#include utils/conditionals.stan
}
data {
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  array[N_obs] int<lower = 0, upper = 1> y_obs;
  int<lower = 1> N_exp;
  array[N_exp] real x_exp;
  vector[N_exp] f_exp;
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}
transformed data {
  int N_full;
  N_full = N_obs + N_exp;
}
parameters {
  // Posterior samples (fitted parameters)
  real<lower = 0> rho;
  real<lower = 0> alpha;
  real<lower = 0> kappa;
  vector[N_obs] f_obs;  // Posterior samples of latent function
}
generated quantities {
  vector[N_pred] f_pred;
  vector[N_pred] y_pred_prob;
  {
    vector[N_full] f_full;
    f_full = append_row(f_obs, f_exp);
    matrix[N_full, N_full] K_full;
    // observational
    K_full[1:N_obs, 1:N_obs] =
      add_diag(gp_matern52_cov(x_obs, alpha, rho), 1e-10);
    // experimental
    K_full[(N_obs + 1):(N_full),
           (N_obs + 1):(N_full)] =
      add_diag(gp_matern52_cov(x_exp, kappa * alpha, rho), 1e-10);
    // cross
    K_full[1:N_obs, (N_obs + 1):(N_full)] =
      gp_matern52_cov(x_obs, x_exp, sqrt(kappa) * alpha, rho);
    K_full[(N_obs + 1):(N_full), 1:N_obs] =
      K_full[1:N_obs, (N_obs + 1):(N_full)]';
    // pred (marginal)
    matrix[N_pred, N_pred] K_pred_pred =
      add_diag(gp_matern52_cov(x_pred, alpha, rho), 1e-10);
    // pred cross
    matrix[N_full, N_pred] K_full_pred;
    K_full_pred[1:N_obs, 1:N_pred] =
      gp_matern52_cov(x_obs, x_pred, alpha, rho);
    K_full_pred[(N_obs + 1):(N_full), 1:N_pred] =
      gp_matern52_cov(x_exp, x_pred, sqrt(kappa) * alpha, rho);
    tuple(vector[N_pred], matrix[N_pred, N_pred]) musigma_cond;
    musigma_cond = gp_cond_pars(f_full, K_full, K_pred_pred, K_full_pred);
    // Sample predictions
    f_pred = multi_normal_rng(musigma_cond.1, musigma_cond.2);
    y_pred_prob = inv_logit(f_pred);
  }
}
