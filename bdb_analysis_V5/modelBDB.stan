
data {
  int<lower = 1> N;                    // number of patients
  int<lower = 1> M;                    // number of covariates + intercept
  int<lower = 1> S;                    // number of studies
  int<lower = 1> NN;                   // number of patients in prosperity
  matrix[N, M] x;                      // covariates
  real<lower = 0> time[N];             // survival times
  int<lower = 0, upper = 1> status[N]; // censoring indicators (1 = event, 0 = censored)
  int<lower = 1, upper = S> study[N];  // studies
  real<lower = 0> nu1;                 // prior: tau2_beta ~ inv_gamma(nu1, nu2)
  real<lower = 0> nu2;                 // prior: tau2_beta ~ inv_gamma(nu1, nu2)
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;

  // Standardize the covariates
  x_std[, 1] = x[, 1];  // intercept does not need standardization
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  vector[M] beta_std[S];
  vector[M] mu_std;
  vector<lower = 0>[M] tau2_beta;
  real<lower = 0> alpha;  // shape parameter for Weibull distribution
}

model {
  for (s in 1:S)
    beta_std[s] ~ normal(mu_std, sqrt(tau2_beta));
  mu_std ~ normal(0, 100);
  tau2_beta ~ inv_gamma(nu1, nu2);
  
  for (n in 1:N) {
    if (status[n] == 1)
      target += weibull_lpdf(time[n] | alpha, exp(x_std[n, ] * beta_std[study[n]]));
    else
      target += weibull_lccdf(time[n] | alpha, exp(x_std[n, ] * beta_std[study[n]]));
  }
  alpha ~ cauchy(0, 2);
}

generated quantities {
  vector[M] beta;
  vector[NN] log_lik;

  beta[1] = beta_std[1, 1];
  for (m in 2:M) {
    beta[m] = beta_std[1, m] / sd_x[m];
    beta[1] -= beta[m] * mu_x[m];
  }

  for (n in 1:NN) {
    if (status[n] == 1)
      log_lik[n] = weibull_lpdf(time[n] | alpha, exp(x_std[n, ] * beta));
    else
      log_lik[n] = weibull_lccdf(time[n] | alpha, exp(x_std[n, ] * beta));
  }
}

