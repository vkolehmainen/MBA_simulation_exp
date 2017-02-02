// This model is for direct inference in full hierarchical model (FHM) in the normal example

data {
  int<lower=0> J;   // number of groups
  int<lower=0> n_j; // number of observations within a group
  vector[J] y[n_j] ;  // array of vectors
}

parameters {
  real<lower=-1, upper=1> glob_mu;
  real eta[J];
  real<lower=0, upper=10> tau;
  real<lower=0, upper=2> sigma[J];
}

transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] <- glob_mu + tau * eta[j];
}

model {
  // priors:
  glob_mu ~ normal(0, 1);
  //sigma ~ uniform(0, 10);
  eta ~ normal(0, 1);
  tau ~ uniform(0, 5);
  
  // likelihood:
  for (j in 1:J) {
      y[j] ~ normal(theta, 1); 
  }        
}
