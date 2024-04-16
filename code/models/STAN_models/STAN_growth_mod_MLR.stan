data {              // This is where you format the input data (Known)
 int<lower=1> N; // Sample size with restriction
 int<lower=1> K; // Number of Predictors
 vector[N] G;       // Response variable
 matrix[N, K] X;
 vector[K] mu_0;
 matrix[K, K] sig20;
}

parameters {        // This is where you state your parameters (Unknown)
 real <lower=0> tau;
 real beta_0;
 vector[K] beta; 
}

transformed parameters{ // This is where you transform your parameters
}                       // Good for hierarchical modeling

model {                     // This is where you specify distributions.
 vector[N] mu;
 
 for(i in 1:N) {
    mu[i] = beta_0 + X[i] * beta;
 }

 // Priors:
 tau ~ uniform( 0, 100 );
 beta_0 ~ normal( 0, 20 );
 beta ~ multi_normal( mu_0, sig20 );
 
 G ~ normal(mu, tau);

}

generated quantities {      // Good for prediction / diagnostic measures
 array[N] real y_rep1;
 array[N] real y_rep2;
 vector[N] mu;

 for(i in 1:N) {
   mu[i] = beta_0 + X[i] * beta;
 }
 y_rep1 = normal_rng(mu, tau);
 y_rep2 = normal_rng(mu, tau);
 
}
