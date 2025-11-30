
data {
  int<lower = 1> N;                     // number of patients
  int<lower = 1> M;                     // number of covariates + intercept
  int<lower = 1> S;                     // number of studies
  int<lower = 1> NN;                    // number of patients in prosperity
  matrix[N, M] x;                       // covariates
  real<lower=0> time[N];                // survival times
  int<lower = 1, upper = 3> study[N];   // studies
  vector<lower = 0, upper = 1>[S] a;    // weight for each study
  int<lower = 0, upper = 1> status[N];  // censoring indicator (0 for censored, 1 for observed)
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
  real<lower = 0> shape;  // Weibull shape parameter
  real<lower = 0> scale;  // Weibull scale parameter
}

model {
  beta_std ~ normal(0, 100);
  shape ~ cauchy(0, 2);
  scale ~ cauchy(0, 2);
  
  for (n in 1:N) {
    if (status[n] == 1) {
      target += a[study[n]] * weibull_lpdf(time[n] | shape, scale * exp(x_std[n, ] * beta_std));
    } else {
      target += a[study[n]] * weibull_lccdf(time[n] | shape, scale * exp(x_std[n, ] * beta_std));
    }
  }
}

generated quantities {
  vector[M] beta;
  vector[NN] log_lik;
  
  // Transform the coefficients back
  beta[1] = beta_std[1];
  for (m in 2:M) {
    beta[m] = beta_std[m] / sd_x[m];
    beta[1] -= beta[m] * mu_x[m];
  }
  
  for (n in 1:NN) {
    if (status[n] == 1) {
      log_lik[n] = weibull_lpdf(time[n] | shape, scale * exp(x[n, ] * beta));
    } else {
      log_lik[n] = weibull_lccdf(time[n] | shape, scale * exp(x[n, ] * beta));
    }
  }
}

