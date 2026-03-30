data {
  int<lower=1> N;                  // Number of observations
  array[N] real x;                 // Temperature data (predictor)
  array[N] int<lower=0, upper=1> y;// Occurrence data (0 or 1)
  // Literature-derived priors passed as data
  real lit_T_opt;                  // Literature value for optimal temp
  real lit_T_opt_sd;               // Uncertainty around the literature T_opt
}
parameters {
  // GP Hyperparameters
  real<lower=0> rho;               // Length-scale (wiggliness)
  real<lower=0> alpha;             // GP marginal standard deviation
  // Latent GP non-centered parameterization trick
  vector[N] eta;                   
  // Mean Function Parameters
  real T_opt;                      // Estimated thermal optimum
  real<lower=0> beta_width;        // Controls the width of the thermal niche
  real alpha_max;                  // Max logit-probability at T_opt
}
transformed parameters {
  vector[N] f;                     // The final latent GP function
  {
    vector[N] mu;                    // The biological mean function
    matrix[N, N] K;                  // Covariance matrix
    matrix[N, N] L_K;                // Cholesky decomposition of K
    // 1. Construct the biological mean function
    for (i in 1:N) {
      mu[i] = alpha_max - beta_width * square(x[i] - T_opt);
    }
    // 2. Build the GP Covariance matrix using a Squared Exponential kernel
    K = add_diag(gp_exp_quad_cov(x, alpha, rho), 1e-10);
    L_K = cholesky_decompose(K);
    // 3. Construct the latent GP (Mean function + GP variation)
    f = mu + L_K * eta;
  }
}
model {
  // --- Priors for the Mean Function ---
  // Informative prior tying the peak of the curve to your literature value
  T_opt ~ normal(lit_T_opt, lit_T_opt_sd); 
  
  // Prior for how wide the niche is. 
  // You will need to tune this based on the distance between T_min and T_max
  beta_width ~ exponential(1);       
  
  // Prior for peak occurrence probability (e.g., normal(0, 1.5) on logit scale)
  alpha_max ~ normal(0, 1.5);        

  // --- Priors for GP Hyperparameters ---
  rho ~ inv_gamma(5, 5);           // Keeps the GP from becoming too wiggly
  alpha ~ std_normal();            // Keeps GP variance reasonable
  eta ~ std_normal();              // Required for the non-centered parameterization
  
  // --- Likelihood ---
  // Connect the latent GP to the binary occurrence data
  y ~ bernoulli_logit(f);          
}
