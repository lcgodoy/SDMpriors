functions {
  // 1. We need the base joint function to reconstruct the Cholesky factor locally
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

  // 2. The Cross-Covariance function for the new grid
  matrix build_K_pred_joint(array[] real x_pred, vector x_obs,
                            vector x_deriv, real alpha, real rho) {
    int N_p = size(x_pred);
    int N_o = rows(x_obs);
    int N_d = rows(x_deriv);
    matrix[N_p, N_o + N_d] K;
    real rho_sq = rho^2;
    real alpha_sq = alpha^2;

    for (i in 1:N_p) {
      for (j in 1:N_o) {
        K[i, j] = alpha_sq * exp(-0.5 * square(x_pred[i] - x_obs[j]) / rho_sq);
      }
      for (j in 1:N_d) {
        real diff = x_pred[i] - x_deriv[j];
        K[i, N_o + j] = (alpha_sq / rho_sq) * diff * exp(-0.5 * square(diff) / rho_sq);
      }
    }
    return K;
  }
}
data {
  // We need the original field data to condition the predictions
  int<lower=1> N_obs;
  vector[N_obs] x_obs;
  int<lower=1> N_deriv;
  vector[N_deriv] x_deriv;  

  // The new high-res prediction grid
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}
parameters {
  // MUST MATCH THE BASE MODEL EXACTLY
  real<lower=0> rho;
  real<lower=0> alpha;
  real baseline_logit;
  real<lower=10, upper=30> T_opt; 
  vector[N_obs + N_deriv] eta_joint; 
}
generated quantities {
  vector[N_pred] f_pred;
  vector[N_pred] p_pred;

  { 
    // 1. Rebuild the base model's Cholesky factor locally for this iteration
    matrix[N_obs + N_deriv, N_obs + N_deriv] K_joint = build_joint_K(x_obs, x_deriv, alpha, rho, 1e-8);
    matrix[N_obs + N_deriv, N_obs + N_deriv] L_joint = cholesky_decompose(K_joint);

    // 2. Build prediction matrices
    matrix[N_pred, N_obs + N_deriv] K_pred_joint = build_K_pred_joint(x_pred, x_obs, x_deriv, alpha, rho);
    matrix[N_pred, N_pred] K_pred_pred = gp_exp_quad_cov(x_pred, alpha, rho);
    for (i in 1:N_pred) K_pred_pred[i, i] += 1e-8; 

    // 3. Condition the multivariate normal
    matrix[N_obs + N_deriv, N_pred] v = mdivide_left_tri_low(L_joint, K_pred_joint');
    vector[N_pred] mu_pred = rep_vector(baseline_logit, N_pred) + v' * eta_joint;
    
    matrix[N_pred, N_pred] cov_pred = K_pred_pred - v' * v;
    for (i in 1:N_pred) cov_pred[i, i] += 1e-8; 
    // 4. Draw the smooth curve
    f_pred = multi_normal_cholesky_rng(mu_pred, cholesky_decompose(cov_pred));   
    for (i in 1:N_pred) {
      p_pred[i] = inv_logit(f_pred[i]);
    }
  }
}
