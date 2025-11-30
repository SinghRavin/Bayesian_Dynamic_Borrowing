
data {
  int<lower = 1> N;                    // number of patients
  int<lower = 1> M;                    // number of covariates + intercept
  int<lower = 1> S;                    // number of studies
  int<lower = 1> NN;                   // number of patients in prosperity
  matrix[N, M] x;                      // covariates
  real<lower=0> time[N];               // survival times
  int<lower = 0, upper = 1> status[N]; // censoring indicator (1 = observed, 0 = censored)
  int<lower = 1, upper = 3> study[N];  // studies
}

transformed data {
  vector[M] mu_x;
  vector[M] sd_x;
  matrix[N, M] x_std;
  
  // x[, 1] is the intercept
  x_std[, 1] = x[, 1];
  for (m in 2:M) {
    mu_x[m] = mean(x[, m]);
    sd_x[m] = sd(x[, m]);
    x_std[, m] = (x[, m] - mu_x[m]) / sd_x[m];
  }
}

parameters {
  vector[M] beta_std;
  real<lower = 0> alpha;              // shape parameter for Weibull
}

model {
  vector[N] log_lambda;

  beta_std ~ normal(0, 100);
  alpha ~ cauchy(0, 2);
  
  log_lambda = x_std * beta_std;

  for (n in 1:N) {
    if (status[n] == 1) {
      target += weibull_lpdf(time[n] | alpha, exp(log_lambda[n]));
    } else {
      target += weibull_lccdf(time[n] | alpha, exp(log_lambda[n]));
    }
  }
}

generated quantities {
  vector[M] beta;
  vector[NN] log_lik;

  for (m in 1:M) {
    beta[m] = beta_std[m];
  }

  for (n in 1:NN) {
    if (status[n] == 1) {
      log_lik[n] = weibull_lpdf(time[n] | alpha, exp(x[n, ] * beta));
    } else {
      log_lik[n] = weibull_lccdf(time[n] | alpha, exp(x[n, ] * beta));
    }
  }
}

