// ---------------------------------------------------------------------------//
// Title:                                                                     //
// Author: Jan Ian Failenschmid                                               //
// Created Date: 08-04-2024                                                   //
// -----                                                                      //
// Last Modified: 30-01-2025                                                  //
// Modified By: Jan Ian Failenschmid                                          //
// -----                                                                      //
// Copyright (c) 2024 by Jan Ian Failenschmid                                 //
// E-mail: J.I.Failenschmid@tilburguniveristy.edu                             //
// -----                                                                      //
// License: GNU General Public License v3.0 or later                          //
// License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html          //
// ---------------------------------------------------------------------------//


data {
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  vector[N_obs] y_obs;
  real<lower=0> alpha;
}

transformed data {
  real xmean = mean(x_obs);
  real xsd = sd(x_obs);
  vector[N_obs] xs = (to_vector(x_obs) - xmean)/xsd;
  array[N_obs] real xn = to_array_1d(xs);
}

parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
}

model {
  rho ~ normal(0, 1);
  sigma ~ normal(0, 1);

  matrix[N_obs, N_obs] cov =  gp_exp_quad_cov(xn, alpha, rho)
  + diag_matrix(rep_vector(square(sigma), N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(cov);
  
  y_obs ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}

generated quantities {
  vector[N_obs] f_predict;
  real gcv_val;
  vector[N_obs] f_post_predict;
  // array[N_obs] real y_predict;
  {
    matrix[N_obs, N_obs] K_x1 = gp_exp_quad_cov(xn, alpha, rho);
    matrix[N_obs, N_obs] K = gp_exp_quad_cov(xn, alpha, rho)
      + diag_matrix(rep_vector(square(sigma), N_obs));
    matrix[N_obs, N_obs] A = mdivide_right_spd(K_x1', K);
    f_predict = A * y_obs;
    gcv_val = N_obs*sum((y_obs - f_predict)^2)/((N_obs - trace(A))^2);
    f_post_predict = multi_normal_rng(f_predict, K_x1 - A * K_x1
      + diag_matrix(rep_vector(1e-10, N_obs))); 
    // y_predict = normal_rng(f_post_predict, sigma);
  }
}
