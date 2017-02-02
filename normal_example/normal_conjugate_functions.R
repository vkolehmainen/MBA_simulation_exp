draw_full_cond <- function(n, obs_mean, obs_var, prior_mean, prior_var) {
  # Draws one sample from the normal conjugate posterior with known variance
  # Source: https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf, Chapter 2
  # Args:
  #   n: amount of observations
  #   obs_mean: mean of observation distribution
  #   obs_var: variance of observation distribution
  #   prior_mean: upper level prior mean
  #   prior_var: upper level prior variance
  new_sigma_sq <- update_sigma_sq(n, obs_var, prior_var)
  new_mu <- update_mu(new_sigma_sq, prior_mean, prior_var, n, obs_mean, obs_var)
  return (rnorm(1, new_mu, sqrt(new_sigma_sq)))
}

update_sigma_sq <- function(n, obs_var, prior_var) {
  # bayesGauss Eq. (20)
  return (1 / (n/obs_var + 1/prior_var))
}

update_mu <- function(sigma_sq_n, prior_mu, prior_var, n, obs_mean, obs_var) {
  # bayesGauss Eq. (24)
  return (sigma_sq_n * (prior_mu/prior_var + n * obs_mean/obs_var))
}