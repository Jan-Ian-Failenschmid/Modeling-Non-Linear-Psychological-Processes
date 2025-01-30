// ---------------------------------------------------------------------------//
// Title:                                                                     //
// Author: Jan Ian Failenschmid                                               //
// Created Date: 18-08-2024                                                   //
// -----                                                                      //
// Last Modified: 29-01-2025                                                  //
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
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
//   real t_diff_min;
//   real t_diff_max;
}

transformed data {
  real xmean = mean(x_obs);
  real ymean = mean(y_obs);
  real xsd = sd(x_obs);
  real ysd = sd(y_obs);
  vector[N_obs] xs = (to_vector(x_obs) - xmean)/xsd;
  vector[N_obs] yn = (y_obs - ymean)/ysd;
  array[N_obs] real xn = to_array_1d(xs);
  array[N_pred] real xp = to_array_1d((to_vector(x_pred) - xmean)/xsd);
  // vector[2] par_guess = [log(10), log(20)]';
  // vector[2] theta = [t_diff_min, t_diff_max]';
  // vector[2] par;
  // array[0] real x_r;
  // array[0] int x_i;

  // par = algebra_solver(tail_delta, par_guess, theta, x_r, x_i);

  // print("a = ", exp(par[1]));
  // print("b = ", exp(par[2]));
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  // rho ~ inv_gamma(100, 250);
  rho ~ lognormal(log(700/xsd), 1);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 1);

  matrix[N_obs, N_obs] cov =  gp_exp_quad_cov(xn, alpha, rho)
  + diag_matrix(rep_vector(square(sigma), N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(cov);
  
  yn ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}

generated quantities {
  vector[N_pred] f_predict;
  real gcv_val;
  vector[N_pred] f_post_predict;
  // array[N_obs] real y_predict;
  {
    matrix[N_obs, N_pred] K_x1 = gp_exp_quad_cov(xn, xp, alpha, rho);
    matrix[N_obs, N_obs] K = gp_exp_quad_cov(xn, alpha, rho)
      + diag_matrix(rep_vector(square(sigma), N_obs));
    matrix[N_pred, N_obs] A = mdivide_right_spd(K_x1', K);
    f_predict = (A * yn)*ysd + ymean;
    // gcv_val = N_obs*sum((y_obs - f_predict)^2)/((N_obs - trace(A))^2);
    // f_post_predict = multi_normal_rng(f_predict, K_x1 - A * K_x1
    //   + diag_matrix(rep_vector(1e-10, N_obs))); 
    // y_predict = normal_rng(f_post_predict, sigma);
  }
}
