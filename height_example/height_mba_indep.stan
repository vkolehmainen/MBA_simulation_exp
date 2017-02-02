// This model is for acquiring first-stage individual posterior distributions in MBA height example
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

  real pr_mean[J];  // hyperparameters
  real pr_sd[J];    // hyperparameters
} 

parameters { 
  real theta[J];
  real<lower=0, upper=15> sigma[J];
}

model { 
  // priors:
  sigma ~ normal(7, 2);
  theta ~ normal(pr_mean, pr_sd);
  
  // likelihood:
  for (j in 1:J) {
      y[start_pos[j]:start_pos[j + 1] - 1]  ~ normal(theta[j], sigma[j]);
  }

  
} 