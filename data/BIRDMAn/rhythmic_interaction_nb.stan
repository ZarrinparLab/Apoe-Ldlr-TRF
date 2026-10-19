data {
    int<lower=1> N;                   // Number of samples
    array[N] int<lower=0> y;          // Count data
    vector[N] log_depths;             // Offset (log of sequencing depth)
    int<lower=1> p;                   // Number of covariates in design matrix
    matrix[N, p] x;                   // Design matrix
}

parameters {
    vector[p] beta_cov;               // All interaction & rhythmic coefficients
    real<lower=0> inv_disp;           // Overdispersion parameter
}

transformed parameters {
    // eta is the "linear predictor" (the log-mean)
    vector[N] eta = log_depths + x * beta_cov; 
    // phi is the precision (reciprocal of overdispersion)
    real phi = 1.0 / inv_disp;
}

model {
    beta_cov ~ normal(0, 2);          // Prior for all coefficients
    inv_disp ~ exponential(1);        // Prior for overdispersion
    
    // Using neg_binomial_2_log is more stable than neg_binomial_2(exp(eta), phi)
    y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
    vector[N] log_lik;
    for (n in 1:N) {
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta[n], phi);
    }
}