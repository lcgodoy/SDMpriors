data {
  int<lower=1> N;                   // Total observations (Field + Pseudo)
  array[N] real x;                  // Temperatures
  array[N] int<lower=0, upper=1> y; // Occurrences (0 or 1)
  vector<lower=0>[N] weights;
}
parameters {
  // GP Hyperparameters
  real<lower=0> rho;               // Length-scale (wiggliness)
  real<lower=0> alpha;             // GP marginal standard deviation
  // Latent GP non-centered parameterization
  vector[N] eta;
  // A simple flat baseline mean (instead of a biological curve)
  real baseline_logit;             
}

transformed parameters {
  vector[N] f;                     // The latent GP function
  {
    matrix[N, N] K;                  // Covariance matrix
    matrix[N, N] L_K;                // Cholesky decomposition
    // Build the Covariance matrix
    K = add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-10);
    L_K = cholesky_decompose(K);
  // Construct the latent GP: Flat Baseline + GP variation
  f = rep_vector(baseline_logit, N) + L_K * eta;
  }
}
model {
  // --- Priors ---
  baseline_logit ~ normal(0, 2);   // Weakly informative prior for the overall mean 
  rho ~ inv_gamma(5, 5);           // Length-scale prior
  alpha ~ std_normal();            // Variance prior
  eta ~ std_normal();              // Non-centered trick
  // --- Likelihood ---
  for (n in 1:N) {
    target += weights[n] * bernoulli_logit_lpmf(y[n] | f[n]);
  }
}
