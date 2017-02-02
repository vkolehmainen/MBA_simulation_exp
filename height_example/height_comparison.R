library(rstan)
library(xtable)
source("normal_conjugate_functions.R")

# This file runs a comparison of 1-stage (FHM) and 2-stage (MBA) Bayesian meta-analyses using simulated data.
# Read README.txt to get an overview of the files needed in the analysis.
#
# Contents:
# 1. DATA INITIALIZATION
# 2. MBA
# 3. FHM
# 4. PLOTS

# 1. DATA INITIALIZATION---------------------------------------------------------------------------------
# the data to be read has to be in the following format, (height_generator.R does this automatically):
#
# group   1       2       3       4  ...  J
#         172.33  180.10  173.45  168.55  182.11
#         NA      NA      175.25  166.53  178.44
#         NA      NA      NA      165.27  179.61
#         NA      NA      NA      170.35  176.54
#         NA      NA      NA      169.26  NA
# n_j_max NA      NA      NA      170.35  NA
# mu_true 175.24  180.77  177.20  168.35  179.35
#
# in short: 1. true values used in data generation should be in the last row
#           2. missing observations are marked as NA
#           3. last row before mu_true should have at least one observation

filename <- "heights_table_example3.txt"
m <- as.matrix(read.csv(filename, sep="\t"))

data <- m[-(dim(m)[1]-1) : -dim(m)[1], ] # raw data excludes the last row

# groups have a varying amount of observations
# -> all the data is combined into a single vector y to help Stan read it

y <- as.vector(data) # convert data to a single vector by each column
y <- y[!is.na(y)] # remove NA values
N <- length(y)  # total number of observations
J <- dim(data)[2]  # number of groups

# calculate start point indexes
start_pos <- rep(1, times=J+1)
for (i in 2:(J+1)) {
  start_pos[i] <- start_pos[i-1] + sum(!is.na(data[,i-1]))
}

# 2. MBA-------------------------------------------------------------------------------------------------
# 1st stage:

# assign prior mean and sd for each group
pr_mean <- as.vector(m[dim(m)[1],]) # true values are used here
pr_sd <- rep(7, times=J)

# calculate individual posterior distributions in STAN
data_list_indep <- list("N", "y", "J", "start_pos", "pr_mean", "pr_sd")
fit_indep <- stan(file = 'height_mba_indep.stan', data = data_list_indep,
                  iter = 1000, chains = 4, control = list(stepsize=0.01, adapt_delta = 0.99), seed = 1)

indep_post <- as.data.frame(fit_indep)
indep_post$lp__ <- NULL   # remove unnecessary column

# calculate sample statistics (mean, sd) from independent posterior samples
mu_sigma <- matrix(NA, nrow=2, ncol=J, dimnames = list(c("mu", "sigma")))
for (j in 1:J) {
  mu_sigma[,j] <- c(mean(indep_post[,j]), sd(indep_post[,j]))
}

# 2nd stage: 
# basic principle of a Gibbs sampler taken from http://stats.stackexchange.com/questions/45946/generating-samples-from-gibbs-sampling (Answer of user10525)
I <- 20999       # total number of iterations
burnin <- 1000   # number of burn-in iterations
thinning <- 10   # select every kth element from posteriors

# initialize containers to store Gibbs sampler's samples of hierarchical model parameters psi and mu
# psi: group specific mean parameter
# mu:  global mean parameter
psi_post <-  matrix(NA, nrow=I, ncol=J)
psi_post[1,] = rnorm(J, 173, 6) # draw initial values from normal distribution
mu_post <- vector(length = I)
mu_post[1] <- 173 # set initial value

# set priors
mu_pr_mean <- 173
mu_pr_var <- 5^2
psi_pr_var <- 6^2

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

# remove burn-in and thin
psi_post <- psi_post[seq(from=burnin, to=I, by=thinning),]
mu_post <- mu_post[seq(from=burnin, to=I, by=thinning)]

# 3. FHM------------------------------------------------------------------------------------------------

# calculate posterior distributions in Stan
data_list_fhm <- list("N", "y", "J", "start_pos")
fit_full <- stan(file = 'height_full.stan', data = data_list_fhm,
                 iter = 1000, chains = 4, control = list(stepsize=0.01, adapt_delta = 0.99), seed = 1)

full_df <- as.data.frame(fit_full)

# 4. PLOTS-----------------------------------------------------------------------------------------------
# note: the code in this section is quite ugly but gets the job done...
# inspect full_df to see the reasoning behind seemingly arbitrary indexes (2+2*J etc.), Stan orders parameters in the order they are introduced in the model
# TODO: clean and optimize this section

mu_true <- as.vector(m[dim(m)[1],]) # extract true mean values to be used in plots

# plot individual posteriors of each group and both methods
for (i in 1:J) {
  plot(density(full_df[[(2+2*J)+i]]), main = paste("1-vaiheinen päättely:", colnames(m)[i]), xlim = c(150, 190), ylim = c(0, 0.3), xlab = "pituus (cm)", ylab = "tiheys")
  abline(v = mu_true[i])
  legend(x="topleft", legend="todellinen arvo", col="black", lty=1, cex=1)
  plot(density(psi_post[,i]), main = paste("2-vaiheinen päättely:", colnames(m)[i]), xlim = c(150, 190), ylim = c(0, 0.3), xlab = "pituus (cm)", ylab = "tiheys")
  abline(v = mu_true[i])
  legend(x="topleft", legend="todellinen arvo", col="black", lty=1, cex=1)
}

# print table of absolute error between posterior mean estimates and true values
abs_error <- matrix(NA, J+1, 2, dimnames = list(c(colnames(m), "Yht."), c("1-vaiheinen", "2-vaiheinen")))
for (i in 1:J) {
  abs_error[i,] <- c(abs(mean(full_df[[(2+2*J)+i]])-mu_true[i]), abs(mean(psi_post[,i])-mu_true[i]))
}
abs_error[J+1,] <- apply(abs_error, 2, sum, na.rm=TRUE)
abs_error

# build a dataframe containing final posteriors of each group and both methods
d <- data.frame(values = c(as.vector(psi_post),  # MBA group distr.
                           mu_post,   # MBA overall distr.
                           as.vector(as.matrix(full_df[,(3+2*J):(2+2*J+J)])),  # FHM group distr.
                           as.vector(full_df[,"glob_mu"])),  # FHM overall distr.
                group = rep(rep(c(colnames(m), "Keskiarvo"), each=2000), times=2), 
                method = rep(c("2-vaiheinen", "1-vaiheinen"), each=2000*(J+1)))

# order the groups by their true values in ascending order
d$group <- factor(d$group, levels = c("Keskiarvo", names(sort(m[dim(m)[1],]))), ordered = TRUE)

true_values <- data.frame(values = c(mu_true, 173),
                          group = c(colnames(m), "Keskiarvo"),
                          method = rep("tod.arvo", times = dim(m)[2]+1)) 

ggplot(d, aes(y=values, x=group, fill = method)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("tomato1", "dodgerblue1", "black")) +
  geom_point(aes(shape = "tod.arvo"), data = true_values, size = 4, shape = 18, show.legend=TRUE) +
  coord_flip() +
  ggtitle("Laatikoiden reunat: 25.- ja 75. persentiilit, keskiviivat: mediaanit") +
  labs(x = "ryhmä", y = "pituus (cm)", fill = "Menetelmä") +
  guides(fill = guide_legend(override.aes=list(shape=c(NA,NA,18), linetype=c(1,1,0), fill=c("tomato1","dodgerblue1","white"))))
  