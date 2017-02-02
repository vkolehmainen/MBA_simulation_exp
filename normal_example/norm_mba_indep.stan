// this model is for acquiring first-stage individual posterior distributions in MBA normal example
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
  int<lower=0> J;   // number of groups
  int<lower=0> n_j; // number of observations within a group
  vector[J] y[n_j] ;  // array of vectors
} 

parameters { 
  real theta[J];
  real<lower=0, upper=10> sigma[J];
}

model { 
  // priors:
  sigma ~ uniform(0, 10);
  theta ~ normal(0, 5);
  
  // likelihood:
  for (j in 1:J) {
      y[j] ~ normal(theta, sigma);
  }
} 
