tuple(vector, matrix) gp_cond_pars(vector z,
                                   matrix sigma,
                                   matrix sigma_p,
                                   matrix sigma_c) { // mean pred
  int N = rows(sigma);
  int Np = rows(sigma_p);
  matrix[N, N] L = cholesky_decompose(sigma);
  matrix[N, Np] v_pred;
  matrix[Np, Np] sigma_cond;
  vector[Np] mu_cond;
  v_pred = mdivide_left_tri_low(L, sigma_c);
  mu_cond = mdivide_right_tri_low(v_pred', L) * z;
  sigma_cond = sigma_p - crossprod(v_pred);
  return (mu_cond, sigma_cond);
}
