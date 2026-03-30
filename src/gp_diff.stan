functions {
  // Build Joint Covariance (Functions + Derivatives)
  matrix build_joint_K(vector x, vector v, real alpha, real rho, real jitter) {
    int N_f = rows(x);
    int N_d = rows(v);
    int N_joint = N_f + N_d;
    matrix[N_joint, N_joint] K;
    real rho_sq = rho^2;
    real alpha_sq = alpha^2;

    for (i in 1:N_f) {
      for (j in i:N_f) {
        K[i, j] = alpha_sq * exp(-0.5 * square(x[i] - x[j]) / rho_sq);
        K[j, i] = K[i, j];
      }
      K[i, i] += jitter;
    }

    for (i in 1:N_d) {
      for (j in i:N_d) {
        real dist_sq = square(v[i] - v[j]);
        K[N_f + i, N_f + j] = (alpha_sq / rho_sq) * exp(-0.5 * dist_sq / rho_sq) * (1.0 - dist_sq / rho_sq);
        K[N_f + j, N_f + i] = K[N_f + i, N_f + j];
      }
      K[N_f + i, N_f + i] += jitter;
    }

    for (i in 1:N_f) {
      for (j in 1:N_d) {
        real diff = x[i] - v[j];
        K[i, N_f + j] = (alpha_sq / rho_sq) * diff * exp(-0.5 * square(diff) / rho_sq);
        K[N_f + j, i] = K[i, N_f + j]; 
      }
    }
    return K;
  }
}
data {
  int<lower=1> N_obs;
  vector[N_obs] x_obs;
  array[N_obs] int<lower=0, upper=1> y_obs;

  int<lower=1> N_deriv;
  vector[N_deriv] x_deriv;  
  real<lower=0> nu;         
  real<lower=0> kappa;      
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real baseline_logit;
  real<lower=10, upper=30> T_opt; 
  vector[N_obs + N_deriv] eta_joint; 
}
transformed parameters {
  vector[N_obs] f_obs;
  vector[N_deriv] d_val;
  {
    matrix[N_obs + N_deriv, N_obs + N_deriv] K_joint;
    matrix[N_obs + N_deriv, N_obs + N_deriv] L_joint;
    vector[N_obs + N_deriv] f_joint;
    K_joint = build_joint_K(x_obs, x_deriv, alpha, rho, 1e-8);
    L_joint = cholesky_decompose(K_joint);
    f_joint = L_joint * eta_joint;
    f_obs = rep_vector(baseline_logit, N_obs) + f_joint[1:N_obs];
    d_val = f_joint[N_obs + 1 : N_obs + N_deriv];
  }
}
model {
  // Stiffer prior on rho to prevent excessive wiggling
  rho ~ inv_gamma(10, 50); 
  alpha ~ std_normal();
  baseline_logit ~ normal(0, 2);
  eta_joint ~ std_normal();
  T_opt ~ normal(20, 10); 
  y_obs ~ bernoulli_logit(f_obs);
  for (i in 1:N_deriv) {
    real dynamic_sign = tanh(kappa * (T_opt - x_deriv[i]));
    target += std_normal_lcdf((dynamic_sign * d_val[i]) / nu);
  }
}
