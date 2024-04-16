functions {
  real beta_4p_lpdf(real y, real a_alpha, real b_alpha, real lower, real upper) {
    // Scale 4-parameter Beta RV to 2-parameter Beta RV
    real x = (y - lower) / (upper - lower);
    
    // Return scaled 2-parameter beta lpdf
    return beta_lpdf(x | a_alpha, b_alpha) - log(upper - lower);
  }
}

data {              // This is where you format the input data (Known)
 int<lower=1> N; // Sample size with restriction
 int<lower=1> K; // Number of Predictors
 vector[N] G;       // Response variable
 vector[N] S;
 matrix[N, K] X;
 vector[K] mu_0;
 matrix[K, K] sig20;
 real c;
 real d;
 real a_alpha;
 real b_alpha;
 real a_gamma;
 real gamma_max;
}

parameters {        // This is where you state your parameters (Unknown)
 real <lower=a_gamma, upper = gamma_max> gamma;    
 // real <lower=0, upper=max_tau2> tau2;
 real <lower=0> tau;
 real beta_0;
 vector[K] beta; 
 real <lower=a_alpha, upper=b_alpha> alpha;
}

transformed parameters{ // This is where you transform your parameters
}                       // Good for hierarchical modeling

model {                     // This is where you specify distributions.
 vector[N] mu;
 for(i in 1:N) {
    mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
 }

 // Priors:
 // gamma ~ gamma( c, d )T[,gamma_max];
 gamma ~ uniform( a_gamma, gamma_max );
 tau ~ uniform( 0, 100 );
 // tau ~ cauchy( 0, 1 )T[0,];
 beta_0 ~ normal( 0, 20 )T[0,];
 beta ~ multi_normal( mu_0, sig20 );
 alpha ~ uniform(a_alpha, b_alpha);
 // target += beta_4p_lpdf(alpha | a_alpha, b_alpha, 0, 5);
 
 // Likelihood:
 G ~  normal( mu, tau );

}

generated quantities {      // Good for prediction / diagnostic measures
 // array[N] real y_rep1;
 // array[N] real y_rep2;
 // vector[N] mu;
 // 
 // for(i in 1:N) {
 //   mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
 // }
 // 
 // y_rep1 = normal_rng(mu, tau);
 // y_rep2 = normal_rng(mu, tau);

}
