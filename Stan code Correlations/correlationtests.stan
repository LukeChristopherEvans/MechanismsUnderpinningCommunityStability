data {
    int<lower=1> N;  // number of observations
    vector[2] x[N];  // input data: rows are observations, columns are the two variables
}

parameters {
    vector[2] mu;
    real<lower=0> sigma[2];
    real<lower=-1, upper=1> rho;  // correlation coefficient
}

transformed parameters {
    // Covariance matrix
    cov_matrix[2] cov = [[      sigma[1] ^ 2       , sigma[1] * sigma[2] * rho],
                         [sigma[1] * sigma[2] * rho,       sigma[2] ^ 2       ]];
}

model {
  // Likelihood
  // using normal
  x ~ multi_normal(mu, cov);

  // standard priors
  sigma ~ exponential(1);
  mu ~ normal(0, 1);

}
