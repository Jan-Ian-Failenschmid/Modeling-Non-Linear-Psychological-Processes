functions {
  #include gpbasisfun_functions.stan
}
data {
  int<lower=1> N;       // number of observations
  int<lower=1> N_eval;  // evaluation locations
  vector[N] x;          // univariate covariate
  vector[N_eval] x_eval;// covariate evaluation locations
  vector[N] y;          // target variable
        
  real<lower=0> c_f1;   // factor c to determine the boundary value L
  int<lower=1> M_f1;    // number of basis functions for smooth function
}
transformed data {

  // Basis functions for f1
  real L_f1 = c_f1*max(x);
  matrix[N, M_f1] PHI_f1 = PHI(N, M_f1, L_f1, x);
  matrix[N_eval, M_f1] PHI_eval = PHI(N_eval, M_f1, L_f1, x_eval);
}
parameters {
  vector[M_f1] beta_f1;         // the basis functions coefficients
  real<lower=0> lengthscale_f1; // lengthscale of f1
  real<lower=0> sigma_f1;       // scale of f1
  real<lower=0> sigma;          // residual scale
}
model {
  // spectral densities for f1
  vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
  // priors
  beta_f1 ~ normal(0, 1);
  lengthscale_f1 ~ inv_gamma(1.2213340, 0.1962925);
  sigma_f1 ~ normal(0, 2);
  sigma ~ normal(0, 1);
  // model
  y ~ normal_id_glm(PHI_f1, 0.0, diagSPD_f1 .* beta_f1, sigma); 
}
generated quantities {
  vector[N_eval] f;
  // vector[N] log_lik;
  {
    // spectral densities for f1
    vector[M_f1] diagSPD_f1 = diagSPD_EQ(sigma_f1, lengthscale_f1, L_f1, M_f1);
    // function scaled back to original scale
    f = (PHI_eval * (diagSPD_f1 .* beta_f1));
    // log_liks for loo
    // for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | f[n], sigma);
  }
}
