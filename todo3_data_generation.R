
# data generation for simulation trials
todo3_data_generation <- function (rseed, iscena, ntrial, m1, nsample, peff.m, theta0, delta1,sigma1,tau){
  
  set.seed(rseed)
  
  NN <- 1000000 # no need to re-generate for each different (y1, y2) combination
  mu1 <- rnorm(NN,qnorm(theta0,0,1),sigma1)
  sigma <- abs(rcauchy(NN, location = 0, scale = tau))
  sigma1 <- abs(rcauchy(NN, location = 0, scale = tau))
  mu2 <- rnorm(NN,mu1,sigma)
  mu3 <- rnorm(NN,mu2,sigma1)
  p1 <- pnorm(mu1)
  p2 <- pnorm(mu2)
  p3 <- pnorm(mu3)
  
  peff <- peff.m[iscena, ]
  narm <- length(peff)
  data1 <- array(0, c(narm, 6, ntrial))
  ia.post1 <- fa.post1 <- fa.post2 <- fa.post3 <- matrix(0, nrow = ntrial, ncol = narm)
  ia.mean1 <- fa.mean1 <- fa.mean2 <- fa.mean3 <- matrix(0, nrow = ntrial, ncol = narm)
  post1 <- post2 <- matrix(0, nrow = ntrial, ncol = 2)
  
  post1f <- function (yy, nn, p1, p2, p3, theta0) {
    
    ll1 <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1)
    ll2 <- yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2)
    ll3 <- yy[3]*log(p3) + (nn[3]-yy[3])*log(1-p3)
    
    if(yy[1]==0){ ll1 <- (nn[1]-yy[1])*log(1-p1) } 
    if(yy[1]>0 & yy[1]<nn[1]){ ll1 <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) } 
    if(yy[1]==nn[1]){ ll1 <- yy[1]*log(p1) } 
    
    if(yy[2]==0){ ll2 <- (nn[2]-yy[2])*log(1-p2) } 
    if(yy[2]>0 & yy[2]<nn[2]){ ll2 <- yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) } 
    if(yy[2]==nn[2]){ ll2 <- yy[2]*log(p2) } 
    
    if(yy[3]==0){ ll3 <- (nn[3]-yy[3])*log(1-p3) } 
    if(yy[3]>0 & yy[3]<nn[3]){ ll3 <- yy[3]*log(p3) + (nn[3]-yy[3])*log(1-p3) } 
    if(yy[3]==nn[3]){ ll3 <- yy[3]*log(p3) }
    
    ll <- ll1 + ll2 + ll3
    
    expll <- exp(ll)
    mean_expll <- mean(expll)
    # ia.post1 <- c(mean(l*(p1>0.2))/mean(l), mean(l*(p2>0.2))/mean(l))
    ia.post1 <- c(mean(expll*(p1>theta0))/mean_expll, mean(expll*(p2>theta0))/mean_expll, mean(expll*(p3>theta0))/mean_expll)
    ia.mean1 <- c(mean(expll*(p1))/mean_expll, mean(expll*(p2))/mean_expll, mean(expll*(p3))/mean_expll)
    
    re.list <- c(post_d1 = ia.post1[1], post_d2 = ia.post1[2], post_d3 = ia.post1[3],
                 mean_d1 = ia.mean1[1], mean_d2 = ia.mean1[2], mean_d3 = ia.mean1[3] )
    return(re.list)
  }
  
  post2f <- function (yy, nn, p1, p2, p3, delta1, theta0) {
    
    ll1 <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1)
    ll2 <- yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2)
    ll3 <- yy[3]*log(p3) + (nn[3]-yy[3])*log(1-p3)
    
    if(yy[1]==0){ ll1 <- (nn[1]-yy[1])*log(1-p1) } 
    if(yy[1]>0 & yy[1]<nn[1]){ ll1 <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) } 
    if(yy[1]==nn[1]){ ll1 <- yy[1]*log(p1) } 
    
    if(yy[2]==0){ ll2 <- (nn[2]-yy[2])*log(1-p2) } 
    if(yy[2]>0 & yy[2]<nn[2]){ ll2 <- yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) } 
    if(yy[2]==nn[2]){ ll2 <- yy[2]*log(p2) } 
    
    if(yy[3]==0){ ll3 <- (nn[3]-yy[3])*log(1-p3) } 
    if(yy[3]>0 & yy[3]<nn[3]){ ll3 <- yy[3]*log(p3) + (nn[3]-yy[3])*log(1-p3) } 
    if(yy[3]==nn[3]){ ll3 <- yy[3]*log(p3) }
    
    ll <- ll1 + ll2 + ll3
    
    expll <- exp(ll)
    mean_expll <- mean(expll)
    ia.post1 <- c(mean(expll*(p1>theta0))/mean_expll, mean(expll*(p2>theta0))/mean_expll, mean(expll*(p3>theta0))/mean_expll)
    ia.mean1 <- c(mean(expll*(p1))/mean_expll, mean(expll*(p2))/mean_expll, mean(expll*(p3))/mean_expll)
    
    post13 <- mean(expll*((p3-p1)<delta1))/mean_expll
    post23 <- mean(expll*((p3-p2)<delta1))/mean_expll
    
    re.list <- c(post_d1 = ia.post1[1], post_d2 = ia.post1[2], post_d3 = ia.post1[3],
                 mean_d1 = ia.mean1[1], mean_d2 = ia.mean1[2], mean_d3 = ia.mean1[3],
                 post1 = post13, post2 = post23)
    return(re.list)
  }
  
  for (trial in 1:ntrial) {
    data1[1:narm, 1, trial] <- c(sum(rbinom(m1, 1, peff[1])), sum(rbinom(m1, 1, peff[2])), sum(rbinom(m1, 1, peff[3])))
    data1[1:narm, 3, trial] <- c(sum(rbinom(nsample - m1, 1, peff[1])), sum(rbinom(nsample - m1, 1, peff[2])), sum(rbinom(nsample - m1, 1, peff[3])))
    data1[1:narm, 5, trial] <- data1[1:narm, 1, trial] + data1[1:narm, 3, trial]
    data1[1:narm, 2, trial] <- m1
    data1[1:narm, 4, trial] <- nsample - m1
    data1[1:narm, 6, trial] <- nsample
    
    
    
    temp_post1f <- post1f(yy = data1[,1,trial], nn = data1[,2,trial], p1, p2, p3, theta0)
    ia.post1[trial, ] <- temp_post1f[1:3]
    ia.mean1[trial, ] <- temp_post1f[4:6]
    
    temp_post2f <- post2f(yy = data1[,5,trial], nn = data1[,6,trial], p1, p2, p3, delta1, theta0)
    fa.post1[trial, ] <- temp_post2f[1:3]
    fa.mean1[trial, ] <- temp_post2f[4:6]
    post1[trial,1] <- temp_post2f[7] ## prob(p3 - p1 < delta1)
    post1[trial,2] <- temp_post2f[8] ## prob(p3 - p2 < delta1)
    
    
    temp_post1f <- post1f(yy = c(data1[1, 1, trial], data1[2, 1, trial], data1[3, 5, trial]), 
                          nn = c(data1[1, 2, trial], data1[2, 2, trial], data1[3, 6, trial]), p1, p2, p3, theta0)
    fa.post2[trial, ] <- temp_post1f[1:3]
    fa.mean2[trial, ] <- temp_post1f[4:6]
    
    
    temp_post2f <- post2f(yy = c(data1[1, 1, trial], data1[2, 5, trial], data1[3, 5, trial]), 
                          nn = c(data1[1, 2, trial], data1[2, 6, trial], data1[3, 6, trial]), p1, p2, p3, delta1, theta0)
    fa.post3[trial, ] <- temp_post2f[1:3]
    fa.mean3[trial, ] <- temp_post2f[4:6]
    post2[trial,1] <- temp_post2f[7] ## prob(p3 - p1 < delta1)
    post2[trial,2] <- temp_post2f[8] ## prob(p3 - p2 < delta1)
  }
  re.list <- list(data = data1, ia.post = ia.post1, fa.post1 = fa.post1, fa.post2 = fa.post2, fa.post3 = fa.post3, 
                  ia.mean = ia.mean1, fa.mean1 = fa.mean1, fa.mean2 = fa.mean2, fa.mean3 = fa.mean3, 
                  post1 = post1, post2 = post2)
  return(re.list)
}