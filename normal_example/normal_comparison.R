library(rstan)
library(xtable)
source("normal_conjugate_functions.R")

# This file runs a comparison of 1-stage (FHM) and 2-stage (MBA) Bayesian meta-analyses using simulated data.
# Read README.txt to get an overview of the files needed in the analysis.
# Note: this is very similar to height example, the difference is that here data is generated before the analysis and not read from a file
# height_comparison.R contains some more detailed comments
#
# Contents:
# 1. DATA INITIALIZATION
# 2. MBA
# 3. FHM
# 4. PLOTS

# 1. DATA INITIALIZATION---------------------------------------------------------
J <- 10           # number of groups
n_j <- 30         # observations per group
global_mu <- 0    # global mean
tau <- 1          # standard deviation of group mean distribution 

# generate vector of groups specific means normally distributed around the global mean with standard deviation std
mu_true <- rnorm(n=J, mean=global_mu, sd=tau)
names(mu_true) <- LETTERS[1:J]

# initialize data matrix y
y <- matrix(NA, nrow=n_j, ncol=J)

# generate each group's data vector to column j
for (j in 1:J) {
  y[,j] <- rnorm(n=n_j, mean=mu_true[j], sd=1)
}

m <- rbind(y, mu_true)

# 2. MBA-------------------------------------------------------------------------
# 1st stage:

# calculate individual posterior distributions in STAN
data_list_1 <- list(y = y, J = J, n_j = n_j)
fit_indep <- stan(file = 'norm_mba_indep.stan', data = data_list_1,
                  iter = 1000, chains = 4, control = list(stepsize=0.01, adapt_delta = 0.99), seed = 1)

indep_post <- as.data.frame(fit_indep)
indep_post$lp__ <- NULL   # remove unnecessary column

# calculate sample statistics (mean, sd) from independent posterior samples
mu_sigma <- matrix(NA, nrow=2, ncol=J, dimnames = list(c("mu", "sigma")))
for (j in 1:J) {
  mu_sigma[,j] <- c(mean(indep_post[,j]), sd(indep_post[,j]))
}

# 2nd stage:
# basic principle of Gibbs sampler taken from http://stats.stackexchange.com/questions/45946/generating-samples-from-gibbs-sampling (Answer of user10525)

I <- 20999      # total number of iterations
burnin <- 1000   # number of burn-in iterations
thinning <- 10   # select every kth element from posteriors

# initialize containers
psi_post <-  matrix(NA, nrow=I, ncol=J)
psi_post[1,] = rnorm(J,0,1) # draw first values from normal distribution
mu_post <- vector(length = I)
mu_post[1] <- 1

# set priors
mu_pr_mean <- 0
mu_pr_var <- 1
psi_pr_var <- 1

# sampling
for (i in 2:I) {  
  #  i = iteration, j = group, function definitions in "normal_conjugate_functions.R"
  mu_post[i] <- draw_full_cond(n = J, 
                               obs_mean = mean(psi_post[i-1,]), 
                               obs_var = var(psi_post[i-1,]),
                               prior_mean = mu_pr_mean, 
                               prior_var = mu_pr_var)
  
  for (j in 1:J) {
    psi_post[i, j] <- draw_full_cond(n = 1,
                                 obs_mean = mu_sigma[1,j], 
                                 obs_var = (mu_sigma[2,j])^2,
                                 prior_mean = mu_post[i], 
                                 prior_var = psi_pr_var)
  }
}

# remove burn-in 
psi_post <- psi_post[seq(from=burnin, to=I, by=thinning),]
mu_post <- mu_post[seq(from=burnin, to=I, by=thinning)]

# 3. FHM-------------------------------------------------------------------------

# calculate posterior distributions in STAN
data_list <- list(y = y, J = J, n_j = n_j)
fit_full <- stan(file = 'norm_full.stan', data = data_list,
                 iter = 1000, chains = 4, control = list(stepsize=0.01, adapt_delta = 0.99), seed = 1)

full_df <- as.data.frame(fit_full)

# 4. PLOTS-----------------------------------------------------------------------

for (i in 1:J) {
  plot(density(full_df[[(2+2*J)+i]]), main = paste("Yksivaiheinen päättely:", LETTERS[i]), xlim = c(-3, 3), ylim = c(0, 2), xlab = "lukuarvo", ylab = "tiheys")
  abline(v = mu_true[i])
  legend(x="topleft", legend="todellinen arvo", col="black", lty=1, cex=1)
  plot(density(psi_post[,i]), main = paste("Kaksivaiheinen päättely:", LETTERS[i]), xlim = c(-3, 3), ylim = c(0, 2), xlab = "lukuarvo", ylab = "tiheys")
  abline(v = mu_true[i])
  legend(x="topleft", legend="todellinen arvo", col="black", lty=1, cex=1)
}

abs_error <- matrix(NA, J+1, 2, dimnames = list(c(LETTERS[1:J], "Yht."), c("1-vaiheinen", "2-vaiheinen")))
for (i in 1:J) {
  abs_error[i,] <- c(abs(mean(full_df[[(2+2*J)+i]])-mu_true[i]), abs(mean(psi_post[,i])-mu_true[i]))
}
abs_error[J+1,] <- apply(abs_error, 2, sum, na.rm=TRUE)
abs_error

d <- data.frame(values = c(as.vector(psi_post),  # MBA group distr.
                           mu_post,   # MBA overall distr.
                           as.vector(as.matrix(full_df[,(3+2*J):(2+2*J+J)])),  # FHM group distr.
                           as.vector(full_df[,"glob_mu"])),  # FHM overall distr.
                group = rep(rep(c(LETTERS[1:J], "Keskiarvo"), each=2000), times=2), 
                method = rep(c("2-vaiheinen", "1-vaiheinen"), each=2000*(J+1)))

d$group <- factor(d$group, levels = c("Keskiarvo", names(sort(mu_true))), ordered = TRUE)

true_values <- data.frame(values = c(mu_true, 0),
                          group = c(LETTERS[1:J], "Keskiarvo"),
                          method = rep("Todellinen arvo", times = J+1)) 

ggplot(d, aes(y=values, x=group, fill = method)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("tomato1", "dodgerblue1", "black")) +
  geom_point(aes(shape = "tod.arvo"), data = true_values, size = 4, shape = 18, show.legend=TRUE) +
  coord_flip() +
  ggtitle("Laatikoiden reunat: 25.- ja 75. persentiilit, keskiviivat: mediaanit") +
  labs(x = "ryhmä", y = "lukuarvo", fill = "Menetelmä") +
  guides(fill = guide_legend(override.aes=list(shape=c(NA,NA,18), linetype=c(1,1,0), fill=c("tomato1","dodgerblue1","white"))))