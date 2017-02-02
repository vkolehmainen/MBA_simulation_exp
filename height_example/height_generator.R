# This file is for generating hierarchical normal data that represents heights in different countries

# INITIALIZATION
J <- 15           # number of groups
global_mu <- 173  # global mean
tau <- 6          # standard deviation of group mean distribution

filename <- "heights_table_example5.txt"

# generate vector of groups specific means normally distributed around the global mean with standard deviation tau
mu_true <- rnorm(n=J, mean=global_mu, sd=tau)

# generate vector of length J to specify amount of datapoints in each group, here they are set manually but feel free to randomize them (e.g. with sample())
n_j <- c(1, 5, 9, 2, 3, 7, 2, 4, 30, 13)

# push values to matrix m
m <- matrix(NA, nrow=max(n_j)+1, ncol=J)
for (j in 1:J) {
  vals <- rnorm(n=n_j[j], mean=mu_true[j], sd=6)
  m[seq(along=vals), j] <- round(vals, digits=3)
}

m[max(n_j)+1,] <- mu_true # add true mean values to last row

# label examples, rename after generation for more realistic labels
colnames(m) <- c("Netherlands", "Vietnam", "Serbia", "Finland", "Thailand", "Germany", "Japan", "France", "Mexico", "Sweden")

write.table(m, file = filename,
            na = "NA",
            append = FALSE, sep = "\t")