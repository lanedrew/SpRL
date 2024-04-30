functions {
  #include skew_generalized_t_helpers.stanfile
}
data {
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
parameters {
  real <lower=200, upper=gamma_max> gamma;    
  real <lower=0, upper=100> tau;
  real <lower=0> beta_0;
  vector[K] beta; 
  real <lower=0.5, upper=5> alpha;
  real <lower=-1, upper=1> lambda;
  real <lower=1> q; 
}
model {
  vector[N] mu;
  
  for(i in 1:N){
    mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
  }
  
  // tau ~ inv_gamma(.001, .001)T[,5];
  tau ~ uniform(0, 100);
  // gamma ~ gamma(.001, .001) T[,gamma_max];
  gamma ~ uniform(200, gamma_max);
  beta_0 ~ normal(0, 20)T[0,];
  beta ~ multi_normal( mu_0, sig20 );
  alpha ~ uniform(0.5, 5);
  lambda ~ normal(0, .25)T[-1,1];
  q ~ gamma(2, 0.2);
  G ~ skew_generalized_t(mu, tau, lambda, 2.0, q);
  
}

generated quantities {      // Good for prediction / diagnostic measures
 array[N] real y_rep1;
 array[N] real y_rep2;
 vector[N] mu;

 for(i in 1:N) {
   mu[i] = ((beta_0 + X[i] * beta) * S[i]^alpha) / (gamma^alpha + S[i]^alpha);
   y_rep1[i] = skew_generalized_t_rng(mu[i], tau, lambda, 2.0, q);
   y_rep2[i] = skew_generalized_t_rng(mu[i], tau, lambda, 2.0, q);
 }

}
