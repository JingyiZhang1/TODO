
# data generation for simulation trials
todo2_data_generation <- function (rseed, iscena, ntrial, m1, nsample, peff.m, theta0, delta1, sigma1, tau) {
  
  
  set.seed(rseed)
  
  NN <- 1000000 # no need to re-generate for each different (y1, y2) combination
  mu1 <- rnorm(NN,qnorm(theta0,0,1),sigma1)
  sigma <- abs(rcauchy(NN, location = 0, scale = tau))
  mu2 <- rnorm(NN,mu1,sigma)
  p1 <- pnorm(mu1)
  p2 <- pnorm(mu2)
  
  post1f <- function (yy, nn, p1, p2, theta0) {
    
    ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
      yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2)
    
    if(yy[1]==0 & yy[2]==0){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        (nn[2]-yy[2])*log(1-p2) } 
    
    if(yy[1]==0 & yy[2]>0 & yy[2]<nn[2]){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) } 
    
    if (yy[1]==0 & yy[2]==nn[2]){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) }
    
    if (yy[1]>0 & yy[1]<nn[1] & yy[2]==0){
      ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
        (nn[2]-yy[2])*log(1-p2) } 
    
    if (yy[1]>0 & yy[1]<nn[1] & yy[2]==nn[2]){
      ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) }
    
    if (yy[1]==nn[1] & yy[2]==0){
      ll <- yy[1]*log(p1) + 
        (nn[2]-yy[2])*log(1-p2) }
    
    if (yy[1]==nn[1] & yy[2]>0 & yy[2]<nn[2]){
      ll <- yy[1]*log(p1) + 
        yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) }
    
    if (yy[1]==nn[1]  & yy[2]==nn[2]){
      ll <- yy[1]*log(p1) + 
        yy[2]*log(p2) } 
    
    expll <- exp(ll)
    # ia.post1 <- c(mean(l*(p1>0.2))/mean(l), mean(l*(p2>0.2))/mean(l))
    ia.post1 <- c(mean(expll*(p1>theta0))/mean(expll), mean(expll*(p2>theta0))/mean(expll))
    ia.mean1 <- c(mean(expll*(p1))/mean(expll), mean(expll*(p2))/mean(expll))
    
    re.list <- c(post_d1 = ia.post1[1], post_d2 = ia.post1[2], mean_d1 = ia.mean1[1], mean_d2 = ia.mean1[2] )
    return(re.list)
  }
  
  post2f <- function (yy, nn, p1, p2, delta1, theta0) {
    
    ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
      yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2)
    
    if(yy[1]==0 & yy[2]==0){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        (nn[2]-yy[2])*log(1-p2) } 
    
    if(yy[1]==0 & yy[2]>0 & yy[2]<nn[2]){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) } 
    
    if (yy[1]==0 & yy[2]==nn[2]){
      ll <- (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) }
    
    if (yy[1]>0 & yy[1]<nn[1] & yy[2]==0){
      ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
        (nn[2]-yy[2])*log(1-p2) } 
    
    if (yy[1]>0 & yy[1]<nn[1] & yy[2]==nn[2]){
      ll <- yy[1]*log(p1) + (nn[1]-yy[1])*log(1-p1) + 
        yy[2]*log(p2) }
    
    if (yy[1]==nn[1] & yy[2]==0){
      ll <- yy[1]*log(p1) + 
        (nn[2]-yy[2])*log(1-p2) }
    
    if (yy[1]==nn[1] & yy[2]>0 & yy[2]<nn[2]){
      ll <- yy[1]*log(p1) + 
        yy[2]*log(p2) + (nn[2]-yy[2])*log(1-p2) }
    
    if (yy[1]==nn[1]  & yy[2]==nn[2]){
      ll <- yy[1]*log(p1) + 
        yy[2]*log(p2) } 
    
    expll <- exp(ll)
    # ia.post1 <- c(mean(l*(p1>0.2))/mean(l), mean(l*(p2>0.2))/mean(l))
    ia.post1 <- c(mean(expll*(p1>theta0))/mean(expll), mean(expll*(p2>theta0))/mean(expll))
    ia.mean1 <- c(mean(expll*(p1))/mean(expll), mean(expll*(p2))/mean(expll))
    
    post1 <- mean(expll*((p2-p1)<delta1))/mean(expll)
    
    re.list <- c(post_d1 = ia.post1[1], post_d2 = ia.post1[2], mean_d1 = ia.mean1[1], mean_d2 = ia.mean1[2], post1 = post1 )
    return(re.list)
  }
  
  data1 <- array(0, c(2, 6, ntrial))
  ia.post1 <- fa.post1 <- fa.post2 <- fa.post3 <- matrix(0, nrow = ntrial, ncol = 2)
  ia.mean1 <- fa.mean1 <- fa.mean2 <- fa.mean3 <- matrix(0, nrow = ntrial, ncol = 2)
  post1 <- rep(0, ntrial)
  peff <- peff.m[iscena, ]
  
  
  for (trial in 1:ntrial) {
    data1[1:2, 1, trial] <- c(sum(rbinom(m1, 1, peff[1])), sum(rbinom(m1, 1, peff[2])))
    data1[1:2, 3, trial] <- c(sum(rbinom(nsample - m1, 1, peff[1])), sum(rbinom(nsample - m1, 1, peff[2])))
    data1[1:2, 5, trial] <- data1[1:2, 1, trial] + data1[1:2, 3, trial]
    data1[1:2, 2, trial] <- m1
    data1[1:2, 4, trial] <- nsample - m1
    data1[1:2, 6, trial] <- nsample
    
    temp_post1f <- post1f(yy = data1[,1,trial], nn = data1[,2,trial], p1, p2, theta0)
    ia.post1[trial, ] <- temp_post1f[1:2]
    ia.mean1[trial, ] <- temp_post1f[3:4]
    
    temp_post2f <- post2f(yy = data1[,5,trial], nn = data1[,6,trial], p1, p2, delta1, theta0)
    fa.post1[trial, ] <- temp_post2f[1:2]
    fa.mean1[trial, ] <- temp_post2f[3:4]
    post1[trial] <- temp_post2f[5]
    
    # dose 1: futile; dose 2: go
    temp_post1f <- post1f(yy = c(data1[1, 1, trial], data1[2, 5, trial]), nn = c(data1[1, 2, trial], data1[2, 6, trial]), p1, p2, theta0)
    fa.post2[trial, ] <- temp_post1f[1:2]
    fa.mean2[trial, ] <- temp_post1f[3:4]
    
    # dose 1: go; dose 2: futile
    temp_post1f <- post1f(yy = c(data1[1, 5, trial], data1[2, 1, trial]), nn = c(data1[1, 6, trial], data1[2, 2, trial]), p1, p2, theta0)
    fa.post3[trial, ] <- temp_post1f[1:2]
    fa.mean3[trial, ] <- temp_post1f[3:4]
  }
  re.list <- list(data = data1, 
                  ia.post = ia.post1, fa.post1 = fa.post1, fa.post2 = fa.post2, fa.post3=fa.post3, post = post1, 
                  ia.mean = ia.mean1, fa.mean1 = fa.mean1, fa.mean2 = fa.mean2, fa.mean3=fa.mean3)
  return(re.list)
}