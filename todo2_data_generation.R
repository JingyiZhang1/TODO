
# data generation for simulation trials
todo2_data_generation <- function (rseed, iscena, ntrial, n1, nsample, peff.m, theta0, delta1, sigma1, tau) {
  
  set.seed(rseed)
  library(rjags)
  
  data1 <- array(0, c(2, 6, ntrial)) # data matrix
  ia.post1 <- fa.post1 <- fa.post2 <- fa.post3 <- matrix(0, nrow = ntrial, ncol = 2) # posterior probabilities
  ia.mean1 <- fa.mean1 <- fa.mean2 <- fa.mean3 <- matrix(0, nrow = ntrial, ncol = 2) # posterior mean estimates
  post1 <- rep(0, ntrial) # posterior probability for between-dose comparison
  peff <- peff.m[iscena, ] # scenario
  
  ## dynamic linear model
  modelstring.ndlm <- " model {
      for(i in 1:narm){
          y[i] ~ dbin(pe[i],n[i])
          probit(pe[i]) <- mu[i]
      }
      mu[1] ~ dnorm(theta,1/sigma1)
      for(i in 2:narm){
          mu[i] ~ dnorm(mu[i-1],1/inv.sig[i-1]/inv.sig[i-1])
          inv.sig[i-1] ~ dt(0,1/tau/tau,1)T(0,)
      }
    }"
  
  ## simulations for data generation and model fitting
  for (trial in 1:ntrial) {
    
    ## data generation
    data1[1:2, 1, trial] <- c(sum(rbinom(n1, 1, peff[1])), sum(rbinom(n1, 1, peff[2])))                      ## number of responders at the interim analysis per arm
    data1[1:2, 3, trial] <- c(sum(rbinom(nsample - n1, 1, peff[1])), sum(rbinom(nsample - n1, 1, peff[2])))  ## number of responders after the interim analysis per arm
    data1[1:2, 5, trial] <- data1[1:2, 1, trial] + data1[1:2, 3, trial]                                      ## total number of responders per arm
    data1[1:2, 2, trial] <- n1             ## number of patients enrolled before interim analysis per arm
    data1[1:2, 4, trial] <- nsample - n1   ## number of patients enrolled after interim analysis per arm
    data1[1:2, 6, trial] <- nsample        ## total number of patients per arm
    
    ## model fitting: interim analysis
    ia.data <- list(y = data1[, 1, trial], n = data1[, 2, trial], narm = dim(data1)[1], theta = qnorm(theta0), sigma1 = sigma1, tau = tau)
    ia <- jags.model(textConnection(modelstring.ndlm), data = ia.data, n.chains = 1, n.adapt = 5000, quiet = TRUE)
    ia.sample <- coda.samples(ia, c("pe"), n.iter = 10000, progress.bar = "none", thin = 5)
    ia1 <- as.matrix(ia.sample)
    ia.post1[trial, ] <- colMeans(ia1 > theta0) ## posterior probability of larger than the null response rate
    ia.mean1[trial, ] <- colMeans(ia1)          ## posterior mean response rate
    
    ## model fitting: final analysis, if both arms pass the interim futility monitoring
    fa.data <- list(y = data1[, 5, trial], n = data1[, 6, trial], narm = dim(data1)[1], theta = qnorm(theta0), sigma1 = sigma1, tau = tau)
    fa <- jags.model(textConnection(modelstring.ndlm), data = fa.data, n.chains = 1, n.adapt = 5000, quiet = TRUE)
    fa.sample <- coda.samples(fa, c("pe"), n.iter = 10000, progress.bar = "none", thin = 5)
    fa1 <- as.matrix(fa.sample)
    fa.post1[trial, ] <- colMeans(fa1 > theta0) ## posterior probability of larger than the null response rate
    fa.mean1[trial, ] <- colMeans(fa1)          ## posterior mean response rate
    post1[trial] <- mean(fa1[, 2] - fa1[, 1] < delta1) ## posterior probability that dose 1 is non-inferior to dose 2
    
    ## model fitting: final analysis, if only dose 2 passes the interim futility monitoring
    fa.data2 <- list(y = c(data1[1, 1, trial], data1[2, 5, trial]), n = c(data1[1, 2, trial], data1[2, 6, trial]), narm = dim(data1)[1], 
                     theta = qnorm(theta0), sigma1 = sigma1, tau = tau)
    fa2 <- jags.model(textConnection(modelstring.ndlm), data = fa.data2, n.chains = 1, n.adapt = 5000, quiet = TRUE)
    fa.sample2 <- coda.samples(fa2, c("pe"), n.iter = 10000, progress.bar = "none", thin = 5)
    fa2 <- as.matrix(fa.sample2)
    fa.post2[trial, ] <- colMeans(fa2 > theta0) ## posterior probability of larger than the null response rate
    fa.mean2[trial, ] <- colMeans(fa2)          ## posterior mean response rate
    
    
    ## model fitting: final analysis, if only dose 1 passes the interim futility monitoring
    fa.data3 <- list(y = c(data1[1, 5, trial], data1[2, 1, trial]), n = c(data1[1, 6, trial], data1[2, 2, trial]), narm = dim(data1)[1], 
                     theta = qnorm(theta0), sigma1 = sigma1, tau = tau)
    fa3 <- jags.model(textConnection(modelstring.ndlm), data = fa.data3, n.chains = 1, n.adapt = 5000, quiet = TRUE)
    fa.sample3 <- coda.samples(fa3, c("pe"), n.iter = 10000, progress.bar = "none", thin = 5)
    fa3 <- as.matrix(fa.sample3)
    fa.post3[trial, ] <- colMeans(fa3 > theta0) ## posterior probability of larger than the null response rate
    fa.mean3[trial, ] <- colMeans(fa3)          ## posterior mean response rate
  }
  re.list <- list(data = data1, 
                  ia.post = ia.post1, fa.post1 = fa.post1, fa.post2 = fa.post2, fa.post3=fa.post3, post = post1,
                  ia.mean = ia.mean1, fa.mean1 = fa.mean1, fa.mean2 = fa.mean2, fa.mean3=fa.mean3)
  return(re.list)
}