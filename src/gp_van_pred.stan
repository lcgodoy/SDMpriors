functions {
#include utils/conditionals.stan
}
data {
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}
transformed data {
}
parameters {
  // Posterior samples (fitted parameters)
  real<lower = 0> rho;
  real<lower = 0> alpha;
  vector[N_obs] f_obs;  // Posterior samples of latent function
}
generated quantities {
  vector[N_pred] f_pred;
  vector[N_pred] y_pred_prob;
  {
    matrix[N_obs, N_obs] K_obs;
    // observational
    K_obs = add_diag(gp_matern52_cov(x_obs, alpha, rho), 1e-10);
    // pred (marginal)
    matrix[N_pred, N_pred] K_pred =
      add_diag(gp_matern52_cov(x_pred, alpha, rho), 1e-10);
    // pred cross
    matrix[N_obs, N_pred] K_cross;
    K_cross = gp_matern52_cov(x_obs, x_pred, alpha, rho);
    tuple(vector[N_pred], matrix[N_pred, N_pred]) musigma_cond;
    musigma_cond = gp_cond_pars(f_obs, K_obs, K_pred, K_cross);
    // Sample predictions
    f_pred = multi_normal_rng(musigma_cond.1, musigma_cond.2);
    y_pred_prob = inv_logit(f_pred);
  }
}
