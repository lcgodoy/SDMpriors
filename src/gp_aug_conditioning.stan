data {
  // Observational (Field) Data
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  array[N_obs] int<lower = 0, upper = 1> y_obs;
  // Experimental (Literature) Data
  int<lower = 1> N_exp;
  array[N_exp] real x_exp;
  vector[N_exp] f_exp;  // We pass the fixed logit values here (e.g., -5.0)
}
parameters {
  real<lower = 0> rho;
  real<lower = 0> alpha;
  vector[N_obs] eta;  // For non-centered parameterization
}
transformed parameters {
  vector[N_obs] f_obs;
  {
    matrix[N_obs, N_obs] K_obs_cond;
    matrix[N_obs, N_obs] L_K_cond;
    vector[N_obs] mu_cond;
    // 1. Build the sub-matrices
    matrix[N_exp, N_exp] K_exp_exp = add_diag(gp_matern52_cov(x_exp, alpha, rho),
                                              1e-16);
    matrix[N_obs, N_obs] K_obs_obs = add_diag(gp_matern52_cov(x_obs, alpha, rho),
                                              1e-16);
    matrix[N_obs, N_exp] K_obs_exp = gp_matern32_cov(x_obs, x_exp, alpha, rho);
    // 2. Perform the conditioning using Cholesky decomposition (faster/stabler
    // than inverse)
    matrix[N_exp, N_exp] L_exp = cholesky_decompose(K_exp_exp);
    matrix[N_exp, N_obs] v = mdivide_left_tri_low(L_exp, K_obs_exp');
    // Calculate Conditional Mean
    mu_cond = K_obs_exp *
      mdivide_left_tri_low(L_exp', mdivide_left_tri_low(L_exp, f_exp));
    // Calculate Conditional Covariance
    K_obs_cond = K_obs_obs - crossprod(v);
    L_K_cond = cholesky_decompose(K_obs_cond);
    // 3. Construct the final latent GP for the observations
    f_obs = mu_cond + L_K_cond * eta;
  }
}
model {
  // Priors
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  eta ~ std_normal();
  // Likelihood (Only evaluated on the real field data!)
  y_obs ~ bernoulli_logit(f_obs);
}
