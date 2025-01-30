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

// functions {
//   vector tail_delta(vector y, vector theta, 
//                     array[] real x_r, array[] int x_i) {
//     vector[2] deltas;
//     deltas[1] = inv_gamma_cdf(theta[1]| exp(y[1]), exp(y[2])) - 0.01;
//     deltas[2] = 1 - inv_gamma_cdf(theta[2]| exp(y[1]), exp(y[2])) - 0.01;
//     return deltas;
//   }
// }

data {
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
}

transformed data {
  real xmean = mean(x_obs);
  real xsd = sd(x_obs);
  vector[N_obs] xs = (to_vector(x_obs) - xmean)/xsd;
  array[N_obs] real xn = to_array_1d(xs);
}

generated quantities {

  real<lower=0> rho = abs(normal_rng(0, 1));
  real<lower=0> alpha = abs(normal_rng(0, 5));

  matrix[N_obs, N_obs] cov = gp_exp_quad_cov(xn, alpha, rho)
  + diag_matrix(rep_vector(1e-6, N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(cov);
  vector[N_obs] f = multi_normal_cholesky_rng(rep_vector(0, N_obs), L_cov);
}
   
