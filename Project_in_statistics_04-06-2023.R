## Install and load required packages

# install.packages("momentuHMM")
# install.packages("tseries")
# install.packages("TruncatedNormal")
# install.packages("mvtnorm")
# install.packages("tmvtnorm")
library(momentuHMM)
library(tseries)
library(TruncatedNormal)
library(mvtnorm)
library(tmvtnorm)

## Simulate realisations from a Markov chain

set.seed(123)
N <- 5 # Number of states
T <- 5000 # Number of realisations
delta <- c(1 / N, 1 / N, 1 / N, 1 / N, 1 / N) # Initial distribution
Gamma <- matrix(c(0.8, 0.05, 0.05, 0.05, 0.05, 
                  0.05, 0.8, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.8, 0.05, 0.05,
                  0.05, 0.05, 0.05, 0.8, 0.05,
                  0.05, 0.05, 0.05, 0.05, 0.8), ncol = 5) # Transition probability matrix
s <- rep(NA, times = T) # State vector
s[1] <- sample(1:N, size = 1, prob = delta)
for(t in 2:T) {
  s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1],])
}

pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- rep(NA, times = T)
for(i in 1:N) {
  cols[s == i] <- pal[i]
}
plot(s, col = cols, xlab = "t", ylab = expression(s[t]))

## Simulate realizations from an HMM

y <- matrix(NA, nrow = T, ncol = 2) # Observation matrix
mu = list()
mu[[1]] <- c(0.1, 0) # Means of the state-dependent distributions
mu[[2]] <- c(0.5, 0.1)
mu[[3]] <- c(0.9, -0.1)
mu[[4]] <- c(1.3, 1)
mu[[5]] <- c(2.5, -1)
sigma <- list()
sigma[[1]] <- matrix(c(0.2, -0.2, -0.2, 0.4), ncol = 2) # Covariance matrices of the state-dependent distributions
sigma[[2]] <- matrix(c(0.2, -0.2, -0.2, 0.4), ncol = 2)
sigma[[3]] <- matrix(c(0.2, -0.2, -0.2, 0.4), ncol = 2)
sigma[[4]] <- matrix(c(0.2, -0.2, -0.2, 0.4), ncol = 2)
sigma[[5]] <- matrix(c(0.2, -0.2, -0.2, 0.4), ncol = 2)

for(t in 1:T) {
  y[t, ] <- TruncatedNormal::rtmvnorm(n = 1, 
                                      mu = mu[[s[t]]], 
                                      sigma = sigma[[s[t]]],
                                      lb = c(0, -Inf),
                                      ub = c(Inf, Inf))
}

par(mfrow = c(2, 1))
plot(x = y[, 1], y = y[, 2], col = cols, xlab = "Horizontal step", ylab = "Vertical step")
plot(s, col = cols, xlab = "t", ylab = expression(s[t]))

hist(y[, 1])
hist(y[, 2])

## Function that computes minus the log-likelihood

par0 <- c(mu0_horizontal_step,
          mu0_vertical_step,
          sigma0_horizontal_step,
          sigma0_vertical_step,
          sigma_horizontal_vertical_step,
          wGamma0, 
          wDelta0)


negative_log_likelihood <- function(Horizontal_step, Vertical_step, par) {
  T <- length(Horizontal_step)
  mu_horizontal_step <- par[1:5] # Unpack parameters
  mu_vertical_step <- par[6:10]
  sigma_horizontal_step <- par[11:15]
  sigma_vertical_step <- par[16:20]
  sigma_horizontal_vertical_step <- par[21:25] # Covariance
  cov_matrix <- list() # Create one covariance matrix for each state
  for(i in 1:5) {
    cov_matrix[[i]] <- matrix(c(sigma_horizontal_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_vertical_step[i]), ncol = 2)
  }
  Gamma <- diag(5) # Diagonal of ones
  Gamma[!Gamma] <- par[26:45] # Fill non-diagonal entries 
  Gamma <- Gamma / rowSums(Gamma) # Divide by row sums
  delta <- c(par[46], par[47], par[48], par[49], 1)
  delta <- delta / sum(delta)
  all_probs <- matrix(1, nrow = T, ncol = 5) # Probabilities of observations conditional on state
  ind <- which(!is.na(Horizontal_step) & !is.na(Vertical_step))
  for(i in 1:5) {
    all_probs[ind, i] <- tryCatch(
      TruncatedNormal::dtmvnorm(x = cbind(Horizontal_step[ind], Vertical_step[ind]), 
                                                         mu = c(mu_horizontal_step[i], mu_vertical_step[i]), 
                                                         sigma = cov_matrix[[i]],
                                                         lb = c(0, -Inf),
                                                         ub = c(Inf, Inf)),
                               error = function(e) NA)  

        #print(cov_matrix[[i]])
  }
  # Forward algorithm (computes the log-likelihood)
  foo <- delta * all_probs[1,] 
  llk <- 0
  for(t in 2:T) {
    foo <- foo %*% Gamma * all_probs[t, ]
    llk <- llk + log(sum(foo))
    foo <- foo / sum(foo) # Scaling to avoid numerical problems
  }
  return(-llk) # Return the negative log-likelihood (because we use a numerical minimizer)
}


Horizontal_step <- prep$step
Vertical_step <- prep$vertical_step

mu_horizontal_step <- par0[1:5] # Unpack parameters
mu_vertical_step <- par0[6:10]
sigma_horizontal_step <- par0[11:15]
sigma_vertical_step <- par0[16:20]
sigma_horizontal_vertical_step <- par0[21:25] # Covariance
cov_matrix <- list() # Create one covariance matrix for each state
for(i in 1:5) {
  cov_matrix[[i]] <- matrix(c(sigma_horizontal_step[i], 
                              sigma_horizontal_vertical_step[i], 
                              sigma_horizontal_vertical_step[i], 
                              sigma_vertical_step[i]), ncol = 2)
}

Horizontal_step[1]
Vertical_step[1]

TruncatedNormal::dtmvnorm(x = cbind(Horizontal_step[1], Vertical_step[1]), 
                          mu = c(mu_horizontal_step[1], mu_vertical_step[1]), 
                          sigma = cov_matrix[[5]],
                          lb = c(0, -Inf),
                          ub = c(Inf, Inf))
set.seed(123)
#############################
# combine the data into a matrix
Horizontal_step <- prep$step
Vertical_step <- prep$vertical_step

data <- cbind(Horizontal_step, Vertical_step/1000)

data[is.na(data)] <- 0.1

# set the number of clusters (states)
k <- 5

# run k-means clustering
km <- kmeans(data, centers = k)

# get the cluster assignments
clusters <- km$cluster

# calculate the mean and standard deviation for each cluster (state) and variable
mu_horizontal_step <- tapply(Horizontal_step, clusters, mean)
mu_vertical_step <- tapply(Vertical_step, clusters, mean)
sigma_horizontal_step <- tapply(Horizontal_step, clusters, sd)
sigma_vertical_step <- tapply(Vertical_step, clusters, sd)

# print the results
print(mu_horizontal_step)
print(mu_vertical_step)
print(sigma_horizontal_step)
print(sigma_vertical_step)
#############################

# timos 


mu0_horizontal_step <- c(0.1, 2, 1, 1.25, 0.75)
mu0_vertical_step <- c(-0.3, 1, 0, -1, 0.5)
sigma0_horizontal_step <- c(0.15, 0.15, 0.3, 0.2, 0.3) ^ 2
sigma0_vertical_step <- c(0.25, 0.35, 0.5, 0.35, 0.3) ^ 2
sigma_horizontal_vertical_step <- c(0, 0, 0, 0, 0)


#mu0_horizontal_step <- c(0.1, 0.3, 0.8, 1.5, 2.0) # Initial values for the mean
#mu0_vertical_step <- c(0.0, 0.1, -0.1, 1, -1)
#mu0_horizontal_step <- c(0.1, 100, 500, 1000, 2000) # Initial values for the mean
#mu0_vertical_step <- c(0.0, 100, -100, 500, -500)

mu0_horizontal_step <- mu_horizontal_step
mu0_vertical_step <- mu_vertical_step
sigma0_horizontal_step <- sigma_horizontal_step
sigma0_vertical_step <- sigma_vertical_step

#sigma0_horizontal_step <- 31*c(0.2, 0.2, 0.2, 0.2, 0.2) # Initial values for the variance
#sigma0_vertical_step <- 31*c(0.4, 0.4, 0.4, 0.4, 0.4)
sigma_horizontal_vertical_step <- c(-0.2, -0.2, -0.2, -0.2, -0.2)# Initial values for the covariance
sigma_horizontal_vertical_step <- rep(0.0, 5)


Gamma0 <- Gamma
wGamma0 <- Gamma0 / diag(Gamma0) # Transform Gamma0 and delta0 to working scale
wGamma0 <- Gamma0[!diag(5)]
delta0 <- delta
wDelta0 <- delta0[-5] / sum(delta0[-5])

par0 <- c(mu0_horizontal_step,
          mu0_vertical_step,
          sigma0_horizontal_step,
          sigma0_vertical_step,
          sigma_horizontal_vertical_step,
          wGamma0,
          wDelta0)

## Check if the negative_log_likelihood() function returns some value

negative_log_likelihood(Horizontal_step = prep$step, 
                        Vertical_step = prep$vertical_step, 
                        par = par0)

## Use nlminb() to minimise the negative log-likelihood (this is equivalent to
## maximising the likelihood)


#mod <- nlminb(start = par0, 
#              objective = negative_log_likelihood, 
#              Horizontal_step = y[, 1],
#              Vertical_step = y[, 2],
#              lower = c(rep(-Inf, 10), rep(0, 10), rep(-Inf, 10), rep(0, 10)),
#              upper = c(rep(Inf, 10), rep(Inf, 10), rep(Inf, 10), rep(Inf, 10)),
#              control = list(eval.max = 1000, 
#                             iter.max = 20, 
#                             trace = 1))

# load data
EagleData <- read.csv('/Users/victorbaekgaard/Desktop/iCloud Drive (Archive)/Documents/Matematik/Kandidat/Project in Stat 23/Eagle_data_with_GPS_positions.csv')

# Rename for prepData. 
colnames(EagleData)[colnames(EagleData) == "Segment_ID"] <- "ID" 

# calculate the vertical steps
EagleData$vertical_step <- unlist(tapply(EagleData$Altitude, EagleData$ID, function(x) c(NA, diff(x))))

EagleData$vertical_step <- EagleData$vertical_step/1000

prep <- prepData(data = EagleData,
                 type = "LL",
                 coordNames = c("Longitude", "Latitude")
)
prep$step[is.na(prep$step)] <- 0.1
prep$vertical_step[is.na(prep$vertical_step)] <- 0.1

missing_values <- which(is.na(prep$step) | is.na(prep$vertical_step))
print(missing_values)

#prep$step <- prep$step / 1000
#prep$vertical_step <- prep$vertical_step / 1000

hist(prep$step)
hist(prep$vertical_step)

mod <- nlminb(start = par0, 
              objective = negative_log_likelihood, 
              Horizontal_step = prep$step,
              Vertical_step = prep$vertical_step,
              lower = c(rep(-Inf, 10), rep(0, 10), rep(-Inf, 10), rep(0, 8)),
              upper = c(rep(Inf, 10), rep(Inf, 10), rep(Inf, 10), rep(Inf, 8)),
              control = list(eval.max = 1000, 
                             iter.max = 1000, 
                             trace = 1))
## Unpack maximum likelihood estimates

mod$par[1:5]
mod$par[6:10]
mod$par[11:15]
mod$par[16:20]
mod$par[21:25]
Gamma_MLE <- diag(5)
Gamma_MLE[!Gamma_MLE] <- mod$par[26:45]
Gamma_MLE <- Gamma_MLE / apply(Gamma_MLE, 1, sum)
Gamma_MLE
delta_MLE <- c(mod$par[46], mod$par[47], mod$par[48], mod$par[49], 1)
delta_MLE <- delta_MLE / sum(delta_MLE)
delta_MLE


sum(Gamma_MLE[1,])
sum(Gamma_MLE[2,])
sum(Gamma_MLE[3,])
sum(Gamma_MLE[4,])
sum(Gamma_MLE[5,])


save(mod, file = "mod.RData")

##########################################################

df_sim_trunc <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())

N <- 5 # Number of states
n_eagles <- 600

for (n in 1:n_eagles) {
  T <- sample(50:100, size = 1)  # Randomly select a value for T between 50 and 500
  
  s <- rep(NA, times = T) # State vector
  
  rows <- list()
  
  s[1] <- sample(1:N, size = 1, prob = delta)
  for (t in 2:T) {
    s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1], ])
  }
  
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- rep(NA, times = T)
  for (i in 1:N) {
    cols[s == i] <- pal[i]
  }
  
  y_new <- matrix(NA, nrow = T, ncol = 2) # Observation vector
 
  for(t in 1:T) {
    y_new[t, ] <- TruncatedNormal::rtmvnorm(n = 1, 
                                            mu = mu[[s[t]]], 
                                            sigma = sigma[[s[t]]],
                                            lb = c(0, -Inf),
                                            ub = c(Inf, Inf))
    
    rows[[t]] <- data.frame("ID" = n, "step" = y_new[t, 1], "vertical_step" = y_new[t, 2], "state" = s[t])
    df_sim_trunc <- rbind(df_sim_trunc, rows[[t]])
  }
}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim_trunc$step, xlab="step length (simulated)", main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim_trunc$vertical_step, xlab="Vertical Step length (simulated)", main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")


Decoded_states <- viterbi(m = mod) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == df_sim_cor1$state) / nrow(simulations_df) # Classification accuracy
sum(Decoded_states != df_sim_cor1$state) / nrow(simulations_df) # Misclassification rate



########## decode with viterbi #########
## Decode the states using the Viterbi algorithm

Viterbi <- function(Horizontal_step, Vertical_step, par) {
  T <- length(Horizontal_step)
  mu_horizontal_step <- par[1:5] # Unpack parameters
  mu_vertical_step <- par[6:10]
  sigma_horizontal_step <- par[11:15]
  sigma_vertical_step <- par[16:20]
  sigma_horizontal_vertical_step <- par[21:25] # Covariance
  Gamma <- diag(5) # Diagonal of ones
  Gamma[!Gamma] <- par[26:45] # Fill non-diagonal entries 
  Gamma <- Gamma / apply(Gamma, 1, sum) # Divide by row sums
  delta <- c(par[46], par[47], par[48], par[49], 1)
  delta <- delta / sum(delta)
  all_probs <- matrix(1, nrow = T, ncol = 5) 
  cov_matrix <- list() # Create one covariance matrix for each state
  for(i in 1:5) {
    cov_matrix[[i]] <- matrix(c(sigma_horizontal_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_vertical_step[i]), ncol = 2)
  }
  ind <- which(!is.na(Horizontal_step) & !is.na(Vertical_step))
  for(i in 1:5) {
    all_probs[ind, i] <- TruncatedNormal::dtmvnorm(x = cbind(Horizontal_step[ind], Vertical_step[ind]), 
                                                   mu = c(mu_horizontal_step[i], mu_vertical_step[i]), 
                                                   sigma = cov_matrix[[i]],
                                                   lb = c(0, -Inf),
                                                   ub = c(Inf, Inf))
  }
  xi  <-  matrix(0, nrow = T, ncol = 5)
  v <- delta * all_probs[1, ]
  xi[1, ] <- v / sum(v)
  for(t in 2:T) {
    v <- apply(xi[t - 1, ] * Gamma, 2, max) * all_probs[t, ] 
    xi[t, ] <- v / sum(v)
  }
  stSeq <- rep(NA, T)
  stSeq[T] <- 4#which.max(xi[T, ])
  
  for (t in (T - 1):1) {
    stSeq[t] <- which.max(Gamma[, stSeq[t + 1]] * xi[t,])
  }
  return(stSeq) 
}

states <- Viterbi(Horizontal_step = df_sim_trunc$step, 
                  Vertical_step = df_sim_trunc$vertical_step,
                  par = mod$par)
states

#############################################
## Compute the pseudo-residuals



pseudo_residuals <- function(Horizontal_step, Vertical_step, par) {
  T <- length(Horizontal_step)
  mu_horizontal_step <- par[1:5] # Unpack parameters
  mu_vertical_step <- par[6:10]
  sigma_horizontal_step <- par[11:15]
  sigma_vertical_step <- par[16:20]
  sigma_horizontal_vertical_step <- par[21:25] # Covariance
  cov_matrix <- list()
  for(i in 1:5) {
    cov_matrix[[i]] <- matrix(c(sigma_horizontal_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_horizontal_vertical_step[i], 
                                sigma_vertical_step[i]), ncol = 2)
  }
  Gamma <- diag(5) # Diagonal of ones
  Gamma[!Gamma] <- par[26:45] # Fill non-diagonal entries 
  Gamma <- Gamma / apply(Gamma, 1, sum) # Divide by row sums
  delta <- c(par[46], par[47], par[48], par[49], 1)
  delta <- delta / sum(delta)
  foo_horizontal_step <- matrix(NA, nrow = T, ncol = 5)
  foo_vertical_step <- matrix(NA, nrow = T, ncol = 5)
  pseudo_residuals_horizontal_step <- rep(NA, T)
  pseudo_residuals_vertical_step <- rep(NA, T)
  log_alpha <- matrix(NA, nrow = T, ncol = 5)
  all_probs <- matrix(1, nrow = T, ncol = 5) 
  
  ind <- which(!is.na(Horizontal_step) & !is.na(Vertical_step))
  for(i in 1:5) {
    all_probs[ind, i] <- TruncatedNormal::dtmvnorm(x = cbind(Horizontal_step[ind], Vertical_step[ind]), 
                                                   mu = c(mu_horizontal_step[i], mu_vertical_step[i]), 
                                                   sigma = cov_matrix[[i]],
                                                   lb = c(0, -Inf),
                                                   ub = c(Inf, Inf))
  }
  lscale <- 0
  v <- delta * all_probs[1, ]
  log_alpha[1, ] <- log(v)
  for(t in 2:T) {
    v <- v %*% Gamma * all_probs[t, ] 
    lscale <- lscale + log(sum(v)) 
    v <- v / sum(v)
    log_alpha[t,] <- log(v) + lscale
  }
  for(i in 1:5) {
    for(t in ind) {
      foo_horizontal_step[t, i] <- tmvtnorm::ptmvnorm.marginal(xn = c(Horizontal_step[t], 0),
                                                               n = 1,
                                                               mean = c(mu_horizontal_step[i], mu_vertical_step[i]), 
                                                               sigma = cov_matrix[[i]],
                                                               lower = c(0, -Inf),
                                                               upper = c(Inf, Inf))[1]
      foo_vertical_step[t, i] <- tmvtnorm::ptmvnorm.marginal(xn = c(0, Vertical_step[t]),
                                                             n = 2,
                                                             mean = c(mu_horizontal_step[i], mu_vertical_step[i]), 
                                                             sigma = cov_matrix[[i]],
                                                             lower = c(0, -Inf),
                                                             upper = c(Inf, Inf))[2]
    }
  }
  pseudo_residuals_horizontal_step[1] <- qnorm(delta %*% foo_horizontal_step[1, ])
  pseudo_residuals_vertical_step[1] <- qnorm(delta %*% foo_vertical_step[1, ])
  for(t in 2:T) {
    c <- max(log_alpha[t - 1,])
    a <- exp(log_alpha[t - 1,] - c)
    pseudo_residuals_horizontal_step[t] <- qnorm(t(a) %*% (Gamma / sum(a)) %*% foo_horizontal_step[t, ])
    pseudo_residuals_vertical_step[t] <- qnorm(t(a) %*% (Gamma / sum(a)) %*% foo_vertical_step[t, ])
    
  }
  return(list(pseudo_residuals_vertical_step = pseudo_residuals_vertical_step,
              pseudo_residuals_horizontal_step = pseudo_residuals_horizontal_step)) 
}

pr <- pseudo_residuals(Horizontal_step = y[, 1], 
                       Vertical_step = y[, 2],
                       par = mod$par)

par(mfrow = c(2, 2))
hist(pr$pseudo_residuals_horizontal_step, main = "Horizontal steps")
hist(pr$pseudo_residuals_vertical_step, main = "Vertical steps")
qqnorm(pr$pseudo_residuals_horizontal_ste, main = "Horizontal steps")
abline(a = 0, b = 1, col = "red")
qqnorm(pr$pseudo_residuals_vertical_step, main = "Vertical steps")
abline(a = 0, b = 1, col = "red")

