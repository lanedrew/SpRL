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
 real <lower=0> nu_beta;
 real <lower=0> nu_delta;
 real mu_delta;
 real <lower=0> sigma_delta;
 real h_1;
 real gamma_max;
}

parameters {        // This is where you state your parameters (Unknown)
 real <lower=200, upper=gamma_max> gamma;    
 real <lower=0, upper=100> tau;
 real <lower=0> beta_0;
 vector[K] beta; 
 real <lower=0.5, upper=5> alpha;
 real delta;
}

transformed parameters{ // This is where you transform your parameters
}                       // Good for hierarchical modeling

model {                     // This is where you specify distributions.
 vector[N] mu;
 
 for(i in 1:N) {
    mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
 }

 // Priors:
 // gamma ~ gamma( c, d ) T[,gamma_max];
 gamma ~ uniform( 200, gamma_max );
 tau ~ uniform( 0, 100 );
 beta_0 ~ normal( 0, 20 )T[0,];
 beta ~ multi_normal( mu_0, sig20 );
 alpha ~ uniform(0.5, 5);
 delta ~ normal( 0, 1 );
 
 G ~ skew_normal(mu, tau, delta);

}

generated quantities {      // Good for prediction / diagnostic measures
 array[N] real y_rep1;
 array[N] real y_rep2;
 vector[N] mu;

 for(i in 1:N) {
   mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
 }
 y_rep1 = skew_normal_rng(mu, tau, delta);
 y_rep2 = skew_normal_rng(mu, tau, delta);
 
}
