# now with 5 states. 
# libraries
library(momentuHMM)
library(tidyverse)
library(tseries)
library(MASS)

# Install and load the tictoc package
#install.packages("tictoc")
library(tictoc)

# Disable the prompt for user input after each plot
devAskNewPage(FALSE)

# load data
EagleData <- read.csv('/Users/victorbaekgaard/Desktop/iCloud Drive (Archive)/Documents/Matematik/Kandidat/Project in Stat 23/Eagle_data_with_GPS_positions.csv')

# Rename for prepData. 
colnames(EagleData)[colnames(EagleData) == "Segment_ID"] <- "ID" 

# calculate the vertical steps
EagleData$vertical_step <- unlist(tapply(EagleData$Altitude, EagleData$ID, function(x) c(NA, diff(x))))

# load as prepData object
#Load using prepData
prep <- prepData(data = EagleData,
                 type = "LL",
                 coordNames = c("Longitude", "Latitude")
)

# plot histogram of the step lengths


hist(prep$step, xlab = "step length", main = "") # between 0 and 4

set.seed(123)

# Starting values for the step length parameters
stepMean0 <- c(0.017, 0.25, 0.5, 0.66, 1.11) # initial means (one for each state)
stepSD0 <- c(0.015, 0.168, 0.301, 0.193, 0.237) # initial standard deviations (one for each state)
stepPar0 <- c(stepMean0, stepSD0)



# Starting values for the step length parameters
#stepMean0 <- c(0.25, 0.017, 1.11, 0.5, 0.66) # initial means (one for each state)
#stepSD0 <- c(0.168, 0.015, 0.237, 0.301, 0.193) # initial standard deviations (one for each state)
#stepPar0 <- c(stepMean0, stepSD0)


#prep$vertical_step <- prep$vertical_step/1000 #ONLY RUN ONCE
#hist(prep$vertical_step[abs(prep$vertical_step) < 0.8], xlab="Vertical Step length", main = "")

vertMean0 <- c(-0.064, -0.026, -0.0017, 0.024, 0.071)
vertSD0 <- c(0.087,0.006,0.027,0.055,0.13)
vertPar0 <- c(vertMean0, vertSD0)

#vertMean0 <- c(0.024, -0.0017, -0.064, 0.071, -0.026)
#vertSD0 <- c(0.055,0.027,0.087,0.13,0.06)
#vertPar0 <- c(vertMean0, vertSD0)


Par0 = list(step = stepPar0,
            vertical_step = vertPar0)

dist = list(step = 'gamma',
            vertical_step = 'norm')


m_5_states <- fitHMM(data = prep, 
                     nbStates = 5,
                     dist = dist,
                     Par0 = Par0,
                     estAngleMean = list(angle = TRUE)
)
# 17 min
 
m_5_states
plot(m_5_states)

save(m_5_states, file = "m_5_states.RData")

load("m_5_states.RData")

plot(m_5_states)

#Not a good fit. Try a combination of different starting values

# Number of tries with different starting values
niter <- 5
# Save list of fitted models
allm_new <- list()

for(i in 1:niter) {
  # Step length mean 
  stepMean0 <- runif(5, 
                     c(0.1, 0.5, 1, 2), 
                     c(0.5, 1.0, 1.5, 3))
  # Step length standard deviation
  stepSD0 <- runif(5,
                   min = c(0.01, 0.5),
                   max = c(0.5, 2)) 

  vertMean <- runif(5,
                     c(-0.01, -0.5),
                     c(0.01, 0.5))
  vertSD0 <- runif(5,
                   c(0.01, 0.1),
                   c(0.1, 0.3)) 
  
  # Fit model
  stepPar0 <- c(stepMean0, stepSD0)
  vertPar0 <- c(vertMean0, vertSD0)
  
  allm_new[[i]] <- fitHMM(data = prep, 
                      nbStates = 5,
                      Par0 = list(step = stepPar0, vertical_step = vertPar0),
                      dist = list(step = "gamma", vertical_step = "norm"), 
                      estAngleMean = list(angle = TRUE))
}

# Extract likelihoods of fitted models

allnllk_new <- unlist(lapply(allm_new, function(m) m$mod$minimum))
allnllk_new

# The higher MLE the better. 
# Index of best fitting model (smallest negative log-likelihood)
whichbest_new <- which.min(allnllk_new)
# Best fitting model
mbest <- allm[[whichbest]]
mbest$mod$minimum

#m_new <- allm_new[[2]]
plot(mbest)
#plot(m_new)


mbest <- m_5_states
save(mbest, file = "mbest.RData")
load("/Users/victorbaekgaard/Desktop/iCloud Drive (Archive)/Documents/Matematik/Kandidat/Project in Stat 23/mbest.RData")

###########################################################################################
#                     SIMULATION OF 600 EAGLES WITH NO CORRELATION                        #
###########################################################################################

N <- 5 # Number of states

n_eagles <- 600

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

df_sim <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())

for (n in 1:n_eagles) {
  T <- sample(50:100, size  = 1)
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
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  
  for (t in 1:T) {
    
    y[t, 1] <- y[t, 1] <- rgamma(n = 1, 
                                 shape = (step_mu[s[t]] / step_sigma[s[t]]) ^ 2, 
                                 rate = step_mu[s[t]] / (step_sigma[s[t]] ^ 2))
    y[t, 2] <- rnorm(n = 1,
                     mean = vert_mu[s[t]], 
                     sd = vert_sigma[s[t]])
    
    rows[[t]] <- data.frame("ID" = n, "step" = y[t, 1], "vertical_step" = y[t, 2], "state" = s[t])
    df_sim <- rbind(df_sim, rows[[t]])
}}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim$step, xlab="step length (simulated)", xlim = c(0, 4), main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

###########################################################################################
#                             Actually looks OK                                           #
###########################################################################################

# fit

Prep_data_NoCor <- prepData(data = df_sim,
                           coordNames = NULL)


mod_NoCor <- fitHMM(data = Prep_data_NoCor, 
                   nbStates = 5,
                   dist = list(step = "gamma", 
                               vertical_step = "norm"),
                   Par0 = list(step = c(step_mean, step_var),
                               vertical_step = c(vert_mu, vert_sigma)))

mod_NoCor

save(mod_NoCor, file = "mod_NoCor.RData")

plot(mod_NoCor)


## Classification accuracy

Decoded_states <- viterbi(m = mod_NoCor) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == df_sim$state) / nrow(df_sim) # Classification accuracy
sum(Decoded_states != df_sim$state) / nrow(df_sim) # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod_NoCor)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod_NoCor)
jarque.bera.test(pseudo_res)



###########################################################################################
#                     SIMULATION OF 600 EAGLES WITH HIGH CORRELATION                      #
###########################################################################################

df_sim_cor1 <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())


N <- 5 # Number of states
n_eagles <- 600

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

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
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  # Define the correlation matrix
  rho <- -0.9  # Desired correlation coefficient
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  # Generate correlated random numbers
  correlated_random_numbers <- MASS::mvrnorm(n = T, mu = c(0, 0), Sigma = Sigma, empirical = TRUE)
  
  for (t in 1:T) {
    # Transform the random numbers to have the desired marginal distributions
    y[, 1] <- qgamma(pnorm(correlated_random_numbers[, 1]), shape = (step_mean[s[t]] / step_var[s[t]]) ^ 2, rate = step_mean[s[t]] / (step_var[s[t]] ^ 2))
    y[, 2] <- qnorm(pnorm(correlated_random_numbers[, 2]), mean = vert_mu[s[t]], sd = vert_sigma[s[t]])
    
    rows[[t]] <- data.frame("ID" = n, "step" = y[t, 1], "vertical_step" = y[t, 2], "state" = s[t])
    df_sim_cor1 <- rbind(df_sim_cor1, rows[[t]])
  }
}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim_cor1$step, xlab="step length (simulated)", xlim = c(0, 4), main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim_cor1$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

##############################################################################
# fit

Prep_data_cor1 <- prepData(data = df_sim_cor1,
                      coordNames = NULL)


mod_cor1 <- fitHMM(data = Prep_data_cor1, 
              nbStates = 5,
              dist = list(step = "gamma", 
                          vertical_step = "norm"),
              Par0 = list(step = c(step_mean, step_var),
                          vertical_step = c(vert_mu, vert_sigma)))

mod_cor1

save(mod_cor1, file = "mod_cor1.RData")

#plot(mod_cor1)


## Classification accuracy

Decoded_states <- viterbi(m = mod_cor1) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == df_sim_cor1$state) / nrow(df_sim_cor1) # Classification accuracy
sum(Decoded_states != df_sim_cor1$state) / nrow(df_sim_cor1) # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod_cor1)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod_cor1)
jarque.bera.test(pseudo_res)

##############################################################################

###########################################################################################
#                     SIMULATION OF 600 EAGLES WITH MEDIUM CORRELATION                    #
###########################################################################################

df_sim_cor2 <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())


N <- 5 # Number of states
n_eagles <- 600

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

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
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  # Define the correlation matrix
  rho <- -0.5  # Desired correlation coefficient
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  # Generate correlated random numbers
  correlated_random_numbers <- MASS::mvrnorm(n = T, mu = c(0, 0), Sigma = Sigma, empirical = TRUE)
  
  for (t in 1:T) {
    # Transform the random numbers to have the desired marginal distributions
    y[, 1] <- qgamma(pnorm(correlated_random_numbers[, 1]), shape = (step_mean[s[t]] / step_var[s[t]]) ^ 2, rate = step_mean[s[t]] / (step_var[s[t]] ^ 2))
    y[, 2] <- qnorm(pnorm(correlated_random_numbers[, 2]), mean = vert_mu[s[t]], sd = vert_sigma[s[t]])
    
    rows[[t]] <- data.frame("ID" = n, "step" = y[t, 1], "vertical_step" = y[t, 2], "state" = s[t])
    df_sim_cor2 <- rbind(df_sim_cor2, rows[[t]])
  }
}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim_cor2$step, xlab="step length (simulated)", breaks = 10, xlim = c(0, 4), main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim_cor2$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

##############################################################################

Prep_data_cor2 <- prepData(data = df_sim_cor2,
                           coordNames = NULL)


mod_cor2 <- fitHMM(data = Prep_data_cor2, 
                   nbStates = 5,
                   dist = list(step = "gamma", 
                               vertical_step = "norm"),
                   Par0 = list(step = c(step_mean, step_var),
                               vertical_step = c(vert_mu, vert_sigma)))

mod_cor2

save(mod_cor2, file = "mod_cor2.RData")

#plot(mod_cor2)
  

## Classification accuracy

Decoded_states <- viterbi(m = mod_cor2) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == df_sim_cor2$state) / nrow(df_sim_cor2) # Classification accuracy
sum(Decoded_states != df_sim_cor2$state) / nrow(df_sim_cor2) # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod_cor2)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod_cor2)
jarque.bera.test(pseudo_res)

##############################################################################

###########################################################################################
#                     SIMULATION OF 600 EAGLES WITH LOW CORRELATION   .                   #
###########################################################################################

df_sim_cor3 <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())


N <- 5 # Number of states
n_eagles <- 600

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

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
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  # Define the correlation matrix
  rho <- -0.2  # Desired correlation coefficient
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  # Generate correlated random numbers
  correlated_random_numbers <- MASS::mvrnorm(n = T, mu = c(0, 0), Sigma = Sigma, empirical = TRUE)
  
  for (t in 1:T) {
    # Transform the random numbers to have the desired marginal distributions
    y[, 1] <- qgamma(pnorm(correlated_random_numbers[, 1]), shape = (step_mean[s[t]] / step_var[s[t]]) ^ 2, rate = step_mean[s[t]] / (step_var[s[t]] ^ 2))
    y[, 2] <- qnorm(pnorm(correlated_random_numbers[, 2]), mean = vert_mu[s[t]], sd = vert_sigma[s[t]])
    
    rows[[t]] <- data.frame("ID" = n, "step" = y[t, 1], "vertical_step" = y[t, 2], "state" = s[t])
    df_sim_cor3 <- rbind(df_sim_cor3, rows[[t]])
  }
}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim_cor3$step, xlab="step length (simulated)", xlim = c(0, 4), main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim_cor3$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

##############################################################################
Prep_data_cor3 <- prepData(data = df_sim_cor3,
                           coordNames = NULL)


mod_cor3 <- fitHMM(data = Prep_data_cor3, 
                   nbStates = 5,
                   dist = list(step = "gamma", 
                               vertical_step = "norm"),
                   Par0 = list(step = c(step_mean, step_var),
                               vertical_step = c(vert_mu, vert_sigma)))

mod_cor3

save(mod_cor3, file = "mod_cor3.RData")

#plot(mod_cor3)


## Classification accuracy

Decoded_states <- viterbi(m = mod_cor3) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == df_sim_cor3$state) / nrow(df_sim_cor3) # Classification accuracy
sum(Decoded_states != df_sim_cor3$state) / nrow(df_sim_cor3) # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod_cor3)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod_cor3)
jarque.bera.test(pseudo_res)





































##############################################################################

hist(prep$step/1000, xlab = "step length (data)", main = "") # between 0 and 4
hist(prep$vertical_step/1000, xlab = "Vertical Step length (data)", main = "")

plot(mbest)



###########################################################################################
#                     SIMULATION OF 600 EAGLES WITH LOW CORRELATION   .                   #
###########################################################################################

df_sim_cor4 <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())


N <- 5 # Number of states
n_eagles <- 600

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

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
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  # Define the correlation matrix
  rho <- 0.15  # Desired correlation coefficient
  Sigma <- matrix(c(1, rho, rho, 1), ncol = 2)
  
  # Generate correlated random numbers
  correlated_random_numbers <- MASS::mvrnorm(n = T, mu = c(0, 0), Sigma = Sigma, empirical = TRUE)
  
  for (t in 1:T) {
    # Transform the random numbers to have the desired marginal distributions
    y[, 1] <- qgamma(pnorm(correlated_random_numbers[, 1]), shape = (step_mean[s[t]] / step_var[s[t]]) ^ 2, rate = step_mean[s[t]] / (step_var[s[t]] ^ 2))
    y[, 2] <- qnorm(pnorm(correlated_random_numbers[, 2]), mean = vert_mu[s[t]], sd = vert_sigma[s[t]])
    
    rows[[t]] <- data.frame("ID" = n, "step" = y[t, 1], "vertical_step" = y[t, 2])
    df_sim_cor4 <- rbind(df_sim_cor4, rows[[t]])
  }
}



par(mfrow = c(1,2), ask = FALSE)

hist(df_sim_cor4$step, xlab="step length (simulated)", xlim = c(0, 4), main = "")
#hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim_cor4$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
#hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

##############################################################################

#HUSK AT TJEK RESULTATERNE!!!

##############################################################################










##############################################################################
library(copula)

df_sim <- data.frame(ID = integer(), step = numeric(), vertical_step = numeric())

N <- 5 # Number of states
n_eagles <- 300

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # transition probabilities matrix

for (n in 1:n_eagles) {
  T <- sample(50:300, size = 1)  # Randomly select a value for T between 50 and 500
  
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
  
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2, ]) # Standard deviations of the state-dependent distributions
  
  vert_mu <- unname(mbest$mle$vertical_step[1, ])
  vert_sigma <- unname(mbest$mle$vertical_step[2, ])
  
  # Specify the marginal distributions
  marginal_dist <- c("gamma", "norm")
  marginal_params <- list(
    list(shape = (step_mean[s] / step_var[s]) ^ 2, rate = step_mean[s] / (step_var[s] ^ 2)),
    list(mean = vert_mu[s], sd = vert_sigma[s])
  )
  
  # Specify the desired correlation coefficient
  rho <- -0.95  # Desired correlation coefficient
  
  # Generate uncorrelated data
  sim_data <- MASS::mvrnorm(T, mu = c(0, 0), Sigma = diag(2))
  
  # Apply inverse marginal CDFs for each state
  for (t in 1:T) {
    sim_data[t, 1] <- qgamma(pnorm(sim_data[t, 1]), shape = marginal_params[[1]]$shape[t], rate = marginal_params[[1]]$rate[t])
    sim_data[t, 2] <- qnorm(pnorm(sim_data[t, 2]), mean = marginal_params[[2]]$mean[t], sd = marginal_params[[2]]$sd[t])
    
    rows[[t]] <- data.frame("ID" = n, "step" = sim_data[t, 1], "vertical_step" = sim_data[t, 2])
    df_sim <- rbind(df_sim, rows[[t]])
  }
}


par(mfrow = c(2,2), ask = FALSE)

hist(df_sim$step, xlab="step length (simulated)", xlim = c(0, 4), main = "")
hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(df_sim$vertical_step, xlab="Vertical Step length (simulated)", xlim = c(-2.5, 1), breaks =4 , main = "")
hist(prep$vertical_step, xlab = "Vertical Step length (data)", main = "")

print(cor(df_sim$step, df_sim$vertical_step))
###########################################################################################
#                             CONCLUSION?                                                 #
###########################################################################################

# try with bivariate truncated normal

#compute the covariance matrix

X <- y[, 1]; Y <- y[, 2]

# Construct covariance matrix
cov_matrix <- matrix(0, nrow = 2, ncol = 2)
cov_matrix[1, 1] <- var(log(X))
cov_matrix[1, 2] <- cov(log(X), Y)
cov_matrix[2, 1] <- cov_matrix[1, 2]
cov_matrix[2, 2] <- var(Y)



test <- matrix(0, nrow = 5, ncol = 2)
for (t in 1:5){
  state <- s[t]
  sigma <- matrix(0, 2, 2)
  sigma[1, 1] <- log(1 + (step_var[state]^2 / step_mean[state]^2))
  sigma[1, 2] <- cov(log(1 + (step_var[state]^2 / step_mean[state]^2)), vert_sigma[state])
  sigma[2, 1] <- sigma[1,2]
  sigma[2, 2] <- vert_sigma[state]
  
  test[t, ] <- TruncatedNormal::rtmvnorm(n = 1, 
                                         mu = c(step_mean[state], vert_mu[state]), 
                                         sigma = cov_matrices[[s[t]]],
                                         lb = c(0, -Inf),
                                         ub = c(Inf, Inf))
  gc()
}



mu = list()
# Calculate covariance matrix
sigma <- list()
for (state in 1:5) {
  mu[[state]] <- c(step_mean[state], vert_mu[state])
  cov_matrix <- matrix(0, nrow = 2, ncol = 2)
  cov_matrix[1, 1] <- step_alpha[state] * vert_sigma[state]^2  # Covariance of gamma-distributed variable with itself
  cov_matrix[2, 2] <- vert_sigma[state]^2  # Variance of normally distributed variable
  cov_matrix[1, 2] <- step_beta[state] * vert_sigma[state]  # Covariance between gamma-distributed and normally distributed variables
  cov_matrix[2, 1] <- cov_matrix[1, 2]  # Covariance matrix is symmetric
  
  sigma[[state]] <- cov_matrix
  
}

##############################
# Calculate correlation matrices
cor_matrices <- list()
for (state in 1:5) {
  cor_matrix <- matrix(0, nrow = 2, ncol = 2)
  cor_matrix[1, 1] <- 1  # Correlation of gamma-distributed variable with itself is always 1
  cor_matrix[2, 2] <- 1  # Correlation of normally distributed variable with itself is always 1
  cor_matrix[1, 2] <- step_beta[state] * vert_sigma[state] / (step_alpha[state] * vert_sigma[state]^2)  # Correlation between gamma-distributed and normally distributed variables
  cor_matrix[2, 1] <- cor_matrix[1, 2]  # Correlation matrix is symmetric
  
  cor_matrices[[state]] <- cor_matrix
}

# Calculate covariance matrices from correlation matrices
cov_matrices <- list()
for (state in 1:5) {
  cov_matrix <- matrix(0, nrow = 2, ncol = 2)
  cov_matrix[1, 1] <- step_alpha[state] * vert_sigma[state]^2  # Variance of gamma-distributed variable
  cov_matrix[2, 2] <- vert_sigma[state]^2  # Variance of normally distributed variable
  cov_matrix[1, 2] <- cor_matrices[[state]][1, 2] * sqrt(cov_matrix[1, 1] * cov_matrix[2, 2])  # Covariance between gamma-distributed and normally distributed variables
  cov_matrix[2, 1] <- cov_matrix[1, 2]  # Covariance matrix is symmetric
  
  cov_matrices[[state]] <- cov_matrix
}

# Print the covariance matrices
for (state in 1:num_states) {
  print(cov_matrices[[state]])
}

##############################

# Print the covariance matrices
for (state in 1:5) {
  print(sigma[[state]])
}

ys <- matrix(NA, nrow = T, ncol = 2) # Observation matrix
for(t in 1:T) {
  ys[t, ] <- TruncatedNormal::rtmvnorm(n = 1, 
                                      mu = mu[[s[t]]], 
                                      sigma = cov_matrices[[s[t]]],
                                      lb = c(0, -Inf),
                                      ub = c(Inf, Inf))
}


# fit

Prep_data <- prepData(data = data.frame(Step = y[1:45, 1],
                                        Vertical_step = y[1:45, 2]),
                      coordNames = NULL)


mod <- fitHMM(data = Prep_data, 
              nbStates = 5,
              dist = list(Step = "gamma", 
                          Vertical_step = "norm"),
              Par0 = list(Step = c(step_mean, step_var),
                          Vertical_step = c(vert_mu, vert_sigma)))

mod

plot(mod)

## Classification accuracy

Decoded_states <- viterbi(m = mod) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == s[1:45]) / 45 # Classification accuracy
sum(Decoded_states != s[1:45]) / 45 # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod)
jarque.bera.test(pseudo_res)



# multiple simulations
n_simulations <- 100

simulations_df <- data.frame() # Data frame to store the simulated data

for (sim in 1:n_simulations) {
  T <- sample(50:300, size = 1) # Random number of data points between 50 and 300
  
  s <- rep(NA, times = T) # State vector
  s[1] <- sample(1:N, size = 1, prob = delta)
  
  for (t in 2:T) {
    s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1], ])
  }
  
  cols <- rep(NA, times = T)
  for (i in 1:N) {
    cols[s == i] <- pal[i]
  }
  
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  
  step_alpha <- (step_mean^2) / (step_var^2)
  step_beta <- step_mean / (step_var^2)
  
  for (t in 1:T) {
    y[t, 1] <- rgamma(n = 1, shape = step_alpha[s[t]], rate = step_beta[s[t]])
    y[t, 2] <- rnorm(n = 1, mean = vert_mu[s[t]], sd = vert_sigma[s[t]])
  }
  
  simulation_df <- data.frame(ID = sim, Step_Length = y[, 1], Vertical_Step = y[, 2], states = s)
  
  simulations_df <- rbind(simulations_df, simulation_df)
}

# Print the resulting data frame
print(simulations_df)

Prep_data_sim <- prepData(data = simulations_df,
                      coordNames = NULL)


 
mod_100 <- fitHMM(data = Prep_data_sim, 
              nbStates = 5,
              dist = list(Step_Length = "gamma", 
                          Vertical_Step = "norm"),
              Par0 = list(Step_Length = c(step_mean, step_var),
                          Vertical_Step = c(vert_mu, vert_sigma)))

mod_100
plot(mod_100)

## Classification accuracy

Decoded_states <- viterbi(m = mod_100) # Viterbi algorithm to compute the most likely state sequence
sum(Decoded_states == simulations_df$states) / nrow(simulations_df) # Classification accuracy
sum(Decoded_states != simulations_df$states) / nrow(simulations_df) # Misclassification rate


## Pseudo-residual analysis

pseudo_res <- pseudoRes(mod_100)$Vertical_stepRes
hist(pseudo_res, freq = FALSE, main = "Pseudo-residuals")
lines(x = seq(from = -3, to = 3, length = 100), y = dnorm (x = seq(from = -3, to = 3, length = 100)), type = "l")
plotPR(mod)
jarque.bera.test(pseudo_res)


# correlation 
library(MASS)


N <- 5 # Number of states
T <- 50000 # Number of realizations

delta <- unname(mbest[["mle"]][["delta"]][1,]) # Initial distribution
Gamma <- unname(mbest[["mle"]][["gamma"]]) # Transition probabilities matrix

round(delta * 100)
round(Gamma * 100)

s <- rep(NA, times = T) # State vector
s[1] <- sample(1:N, size = 1, prob = delta)

for (t in 2:T) {
  s[t] <- sample(1:N, size = 1, prob = Gamma[s[t - 1], ])
}

pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cols <- rep(NA, times = T)

for (i in 1:N) {
  cols[s == i] <- pal[i]
}

plot(s, col = cols, xlab = "t", ylab = expression(s[t]))

# Define correlation levels
correlation_levels <- c(0, 0.3, 0.7)

# Simulate data with different correlation levels
for (cor_level in correlation_levels[2]) {
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  
  # Simulate the first variable with state-dependent distribution
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2,]) # Standard deviations of the state-dependent distributions
  step_alpha <- (step_mean^2) / step_var^2
  step_beta <- step_mean / step_var^2
  
  for (t in 1:T) {
    y[t, 1] <- rgamma(n = 1, shape = step_alpha[s[t]], rate = step_beta[s[t]])
  }
  
  # Simulate the second variable with normal distribution and correlation
  vert_mu <- unname(mbest$mle$vertical_step[1,])
  vert_sigma <- unname(mbest$mle$vertical_step[2,])
  
  # Generate correlated samples
  correlated_samples <- mvrnorm(T, mu = c(0, 0), Sigma = matrix(c(1, cor_level, cor_level, 1), nrow = 2))
  
  for (t in 1:T) {
    y[t, 2] <- rnorm(n = 1, mean = vert_mu[s[t]], sd = vert_sigma[s[t]]) + correlated_samples[t, 2]
  }
  
  # Perform analysis with the simulated data (e.g., calculate correlations, plot, etc.)
  # ...
}


par(mfrow = c(2,2), ask = FALSE)

hist(y[, 1][abs(y[, 1]) < 4], xlab="step length (simulated)", xlim = c(0, 4), main = "")
hist(prep$step, xlab = "step length (data)", main = "") # between 0 and 4

hist(y[, 2], xlab="Vertical Step length (simulated)", main = "")
hist(prep$vertical_step[abs(prep$vertical_step) < 0.5], xlab = "Vertical Step length (data)", main = "")


cor(y[, 1], y[, 2])

for (cor_level in correlation_levels[3]) {
  y <- matrix(NA, nrow = T, ncol = 2) # Observation vector
  
  # Simulate the first variable with state-dependent distribution
  step_mean <- unname(mbest$mle$step[1,])  # Means of the state-dependent distributions
  step_var <- unname(mbest$mle$step[2,]) # Standard deviations of the state-dependent distributions
  step_alpha <- (step_mean^2) / step_var^2
  step_beta <- step_mean / step_var^2
  
  for (t in 1:T) {
    y[t, 1] <- rgamma(n = 1, shape = step_alpha[s[t]], rate = step_beta[s[t]])
  }
  
  # Simulate the second variable with normal distribution and correlation
  vert_mu <- unname(mbest$mle$vertical_step[1,])
  vert_sigma <- unname(mbest$mle$vertical_step[2,])
  
  # Generate uncorrelated samples for the second variable
  uncorrelated_samples <- rnorm(T, mean = 0, sd = vert_sigma[s])
  
  # Calculate the Cholesky decomposition of the correlation matrix
  cor_matrix <- matrix(c(1, cor_level, cor_level, 1), nrow = 2)
  chol_decomp <- chol(cor_matrix)
  
  # Generate correlated samples for the second variable using the Cholesky decomposition
  correlated_samples <- uncorrelated_samples %*% chol_decomp
  
  for (t in 1:T) {
    # Add correlated samples to the second variable
    y[t, 2] <- vert_mu[s[t]] + correlated_samples[t, 2]
  }
  
  # Calculate the correlation between the two variables
  correlation <- cor(y[, 1], y[, 2])
  
  # Perform analysis with the simulated data (e.g., plot, calculate correlations, etc.)
  # ...
  
  print(paste("Correlation (", cor_level, "):", correlation))
}
