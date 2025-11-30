# R Libraries we need.

library(simcausal)
library(dplyr)
library(sas7bdat)
library(haven)
library(survival)
library(ggmcmc)
library(robustHD)
library(ggpubr)
library(MASS)
library(rstan)
library(loo)
library(rstudioapi)
library(doParallel)
library(parallel)
library(HDInterval)
library(progress)
library(pbapply)

# Some constants

rstan_iters <- 10000 # Number of posterior draws per chain (5,000 are burn draws)
n.chains <- 3 # Number of posterior markov chains.
m <- 7 # Number of covariates + intercept.

# Functions for the estimation bias/error/prediction/accuracy.

log.MSE <- function(y.true, y.pred){log((y.true-y.pred)^2)}
per.bias <- function(y.true, y.pred){(y.pred-y.true)/y.true}

extract.stan <- function(fit, M = m, coef, pool.var, total.obs) {
  est <- summary(fit)$summary
  beta <-  est[grep('^beta\\[', rownames(est)), ][1:M, ]
  
  # Variance of the posterior distribution of each effect.
  var <- beta[,"sd"]^2
  
  # 95% HPD credible interval.
  hpd_intervals_95 = hdi(extract(fit)$beta)
  
  # Calculate coverage for 95% HPD credible interval
  coverage_hpd_95 <- sapply(1:M, function(i) {
    coef[i] >= hpd_intervals_95["lower", i] && coef[i] <= hpd_intervals_95["upper", i]
  })
  
  list(looic = loo(fit)$estimates[2:3, 1],
       log.MSE = log.MSE(coef, beta[, 1]),
       per.bias = per.bias(coef, beta[, 1]),
       coverage_hpd_95 = coverage_hpd_95,
       width_hpd_95 = hpd_intervals_95 |> diff(),
       TESS = total.obs*(pool.var/var))
}

# Function to set the data into a stan format i.e., list.

data.stan <- function(dat_combined) {
  dat <- dat_combined[rank(dat_combined$Study, ties.method = 'first'), ]
  x <- cbind(Intercept = 1, dat[, 3:8])
  status <- dat$composite
  time <- dat$follow_up_time_composite_in_days_RS
  study <- dat$Study
  N <- nrow(x) # number of patients
  M <- ncol(x) # number of covariates + intercept
  S <- length(unique(study)) # number of studies
  NN <- sum(study == 3) # number of patients in prosperity
  list(N = N, M = M, S = S, NN = NN, x = x, status = status, time = time, study = study)
}

# Model string for complete pooling and non-informative borrowing methods.

modelstring <- '
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

  // Transform the coefficients back
  beta[1] = beta_std[1];
  for (m in 2:M) {
    beta[m] = beta_std[m] / sd_x[m];
    beta[1] -= beta[m] * mu_x[m];
  }

  for (n in 1:NN) {
    if (status[n] == 1) {
      log_lik[n] = weibull_lpdf(time[n] | alpha, exp(x[n, ] * beta));
    } else {
      log_lik[n] = weibull_lccdf(time[n] | alpha, exp(x[n, ] * beta));
    }
  }
}
'

# Function to suppress output
suppress_output <- function(expr) {
  capture.output(result <- expr, file = NULL)
  return(result)
}

pool <- function(dat_combined) {
  writeLines(modelstring, con = 'model_pool.stan')
  data <- data.stan(dat_combined)
  suppress_output(stan('model_pool.stan', data = data, iter = rstan_iters, chains = n.chains))
}


non_inf <- function(dat_combined) {
  writeLines(modelstring, con = 'model_non_inf.stan')
  data <- data.stan(subset(dat_combined, Study == 3))
  suppress_output(stan('model_non_inf.stan', data = data, iter = rstan_iters, chains = n.chains))
}

# Model string for static power priors method.

modelstring.pp <- '
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
'

pp <- function(dat_combined, a = 0.5) {
  writeLines(modelstring.pp, con = 'model_pp.stan')
  data <- data.stan(dat_combined)
  data$a <- c(rep(a, data$S - 1), 1) 
  suppress_output(stan('model_pp.stan', data = data, iter = rstan_iters, chains = n.chains))
}

# Model string for Bayesian dynamic borrowing methods.

modelstring.bdb <- '
data {
  int<lower = 1> N;                    // number of patients
  int<lower = 1> M;                    // number of covariates + intercept
  int<lower = 1> S;                    // number of studies
  int<lower = 1> NN;                   // number of patients in prosperity
  matrix[N, M] x;                      // covariates
  real<lower = 0> time[N];             // survival times
  int<lower = 0, upper = 1> status[N]; // censoring indicators (1 = event, 0 = censored)
  int<lower = 1, upper = 3> study[N];  // studies
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
'

bdb <- function(dat_combined, nu1 = 0.001, nu2 = 0.001) {
  writeLines(modelstring.bdb, con = 'model_bdb.stan')
  data <- data.stan(dat_combined)
  data$nu1 <- nu1
  data$nu2 <- nu2
  suppress_output(stan('model_bdb.stan', data = data, iter = rstan_iters, chains = n.chains))
}






