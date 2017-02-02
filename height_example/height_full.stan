// This model is for direct inference in full hierarchical model (FHM) in the height example
// all N observations are stored in a single vector y and start_pos contains the indexes of first observations of each group in y
// e.g if data is:
// group  observations 
//   1    6 2 6 
//   2    3 3 7 1 5 9
//   3    5 7 4 2 
// then:
//         y = c(6, 2, 6, 3, 3, 7, 1, 5, 9, 5, 7, 4, 2)
// start_pos = c(1, 4, 10)

data {
  int N;  // total number of observations
  vector[N] y;  // observations
  int<lower=0> J;  // number of groups  
  int<lower=1> start_pos[J+1];  // the first position on observations for each group 
}

parameters {
  real<lower=0> glob_mu;
  real eta[J];
  real<lower=0, upper=20> tau;
  real<lower=0, upper=15> sigma[J];
}

transformed parameters {
  // here each theta is reparameterized as glob_mu + tau * eta, eta ~ normal(0, 1) to facilitate better convergence
  // see Betancourt & Girolami (2013) for more information
  real theta[J];
  for (j in 1:J)
    theta[j] <- glob_mu + tau * eta[j];
}

model {
  // priors:
  glob_mu ~ normal(173, 5);
  sigma ~ uniform(0, 15);
  eta ~ normal(0, 1);
  
  // likelihood:
  for (j in 1:J) {
    y[start_pos[j]:start_pos[j + 1] - 1]  ~ normal(theta[j], sigma[j]);
  }        
}