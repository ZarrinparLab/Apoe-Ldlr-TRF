data {
    int<lower=1> N;                   // Number of samples
    array[N] int<lower=0> y;          // Feature counts (one microbe)
    vector[N] log_depths;             // Offset (log of sequencing depth)
    vector[N] time;                   // Time points
    real<lower=0> period;             // Fixed period (e.g., 24 for circadian)
    
    // Additional covariates (e.g., Treatment, Age)
    int<lower=0> p;                   
    matrix[N, p] x;                   
}

parameters {
    real alpha;                       // Mesor (Mean value of rhythm)
    real beta_sin;                    // Sine coefficient
    real beta_cos;                    // Cosine coefficient
    vector[p] beta_cov;               // Other covariate coefficients
    real<lower=0> inv_disp;           // 1/phi (Inverse dispersion)
}

transformed parameters {
    vector[N] mu;
    real phi = 1.0 / inv_disp;
    
    for (n in 1:N) {
        // Log-link: log(mu) = offset + mesor + sin_term + cos_term + covariates
        mu[n] = exp(
            log_depths[n] + alpha + 
            beta_sin * sin(2 * pi() * time[n] / period) + 
            beta_cos * cos(2 * pi() * time[n] / period) + 
            x[n] * beta_cov
        );
    }
}

model {
    // Priors (Standard BIRDMAn/Weakly Informative)
    alpha ~ normal(0, 5);
    beta_sin ~ normal(0, 2);
    beta_cos ~ normal(0, 2);
    beta_cov ~ normal(0, 2);
    inv_disp ~ exponential(1);

    y ~ negative_binomial_2(mu, phi);
}

generated quantities {
    // Convert sin/cos coefficients to MetaCycle parameters
    real amplitude = sqrt(beta_sin^2 + beta_cos^2);
    
    // Phase in radians (atan2 ensures correct quadrant)
    real phase_rad = atan2(beta_sin, beta_cos); 
    
    // Phase in time units (e.g., hours)
    real phase_time = (phase_rad / (2 * pi())) * period;
    if (phase_time < 0) phase_time += period;

    // Log-likelihood for BIRDMAn diagnostics/LOO
    vector[N] log_lik;
    for (n in 1:N) {
        log_lik[n] = negative_binomial_2_lpmf(y[n] | mu[n], phi);
    }
}