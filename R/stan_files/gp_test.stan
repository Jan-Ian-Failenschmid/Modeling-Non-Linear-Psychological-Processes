// ---------------------------------------------------------------------------//
// Title:                                                                     //
// Author: Jan Ian Failenschmid                                               //
// Created Date: 14-06-2024                                                   //
// -----                                                                      //
// Last Modified: 25-06-2024                                                  //
// Modified By: Jan Ian Failenschmid                                          //
// -----                                                                      //
// Copyright (c) 2024 by Jan Ian Failenschmid                                 //
// E-mail: J.I.Failenschmid@tilburguniveristy.edu                             //
// -----                                                                      //
// License: GNU General Public License v3.0 or later                          //
// License URL: https://www.gnu.org/licenses/gpl-3.0-standalone.html          //
// ---------------------------------------------------------------------------//



functions {
  // Smoothing spline covariance basis
  real cov_fun_sp(real x1, real x2) {
    real dist = abs(x1 - x2);
    real v = min([ x1, x2 ]);
    real k_sp = (dist * v^2) / 2 + v^3 / 3;
    return k_sp;
  }

  // Function for creating a custom covariance matrix
  matrix cust_cov(array[] real x1, array[] real x2, int N1, int N2) {
    
    matrix[N1, N2] cov_mat;
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        cov_mat[i, j] = cov_fun_sp(x1[i], x2[j]);
      }
    }

    return cov_mat; 
  }
}

data {
  // Observations
  int<lower=1> N_obs;
  array[N_obs] real x_obs;
  vector[N_obs] y_obs;

  // Predictions
  int<lower=1> N_pred;
  array[N_pred] real x_pred;
}

transformed data {
   // Create H matrix for mean basis functions
   matrix[N_obs, 2] H = append_col(rep_vector(1, N_obs), to_vector(x_obs));
   matrix[N_obs, 2] H_pred = append_col(rep_vector(1, N_pred), to_vector(x_pred));
}

parameters {
  real<lower=0> sigma_f;
  real<lower=0> sigma_n;
  vector[2] beta;
}

model {
  // Create covariance matrix
  matrix[N_obs, N_obs] K = sigma_f * cust_cov(x_obs, x_obs, N_obs, N_obs)
   + diag_matrix(rep_vector(square(sigma_n), N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(K);
  // Likelihood
  y_obs ~ multi_normal_cholesky(H*beta, L_cov);
  // real log_p;
  // {
  //   matrix[N_obs, N_obs] K_inv = inverse_spd(K);
  //   matrix[2, 2] A = H' * K_inv * H;
  //   matrix[N_obs, N_obs] C = K_inv * H * inverse_spd(A) * H' * K_inv; 
  //   log_p = - 0.5 * y_obs' * K_inv * y_obs + 0.5 * y_obs' * C * y_obs
  //     - 0.5 * log_determinant(K_inv) - 0.5 * log_determinant(A) 
  //     - ((N_obs - 2) / 2.0) * log(2 * pi());
  // }
  // target += log_p;
}

generated quantities {
   // Generate predictions for f
   vector[N_pred] f;
   {  
      // Create component matrices
      // Covariance Matrix
      matrix[N_obs, N_obs] K = sigma_f * cust_cov(x_obs, x_obs, N_obs, N_obs)
        + diag_matrix(rep_vector(square(sigma_n), N_obs));
      matrix[N_obs, N_obs] L_cov = cholesky_decompose(K);
      // Mean function
      matrix[N_obs, N_pred] k_x1_x2 = sigma_f * cust_cov(x_obs, x_pred, N_obs, N_pred);
      matrix[N_pred, 2] R = H_pred - ((H' / K) * k_x1_x2)';
      vector[2] beta_hat = ((H' / K) * H) \ ((H' / K) * y_obs);

      // Invert K through the cholesky decomposition
      vector[N_obs] L_K_div_y1 = mdivide_left_tri_low(L_cov, y_obs);
      vector[N_obs] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_cov)';
      f = (R * beta_hat) + (k_x1_x2' * K_div_y1);
   } 
}
