source('E:/projects/20221213-project-TODO/review/20241119/code/todo3_data_generation.R')
source('E:/projects/20221213-project-TODO/review/20241119/code/todo3_futility.R')

## simulation scenarios
peff.m <- rbind(c(0.20, 0.20, 0.20), c(0.40, 0.40, 0.40) )

## number of scenarios
sim.nscenarios <- dim(peff.m)[1]

## null response rate
theta0 <- 0.2

## non-inferior margin
delta1 <- 0.05

## number of simulation trials
ntrial <- 10000

## number of cutoff combinations for futility monitoring
ntriala <- 4200

## family-wise type I error rate (FWER)
target.alpha <- 10
## target overall power
target.op <- 95

## size of the inconclusive region size (SIR)
sir <- 15

## maximum rate of selecting an inadequate dose (MRID)
mrid <- 25

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.6

for(nsample in 21:22){ ## nsample: sample size per arm
  
  m1_c <- round((nsample/3):(nsample*2/3))
  ass_m1 <- cbind(m1_c, rep(0, length(m1_c)))
  
  ffr2 <- list()
  ffr3 <- ffr4 <- ffr5 <- NULL
  ia12_m1 <- list()
  
  for(im1 in m1_c){
    
    print(paste('n=', nsample, ', m1=', im1,sep=''))
    
    todo3_data <- temp <- list()
    for(iscena in 1:2){
      print(paste('Scenario', iscena))
      todo3_data[[iscena]] <- todo3_data_generation(rseed=1, iscena, ntrial, m1 = im1, nsample, peff.m, theta0, delta1, sigma1, tau)
    }
    
    todo3_a <- tempa <- list()
    for(iscena in 1:2){
      for(iia in 1:ntriala){
        tempa[[iia]] <- todo3_futility(ia12 = iia, fdata=todo3_data[[iscena]]$data, 
                                       ia.post=todo3_data[[iscena]]$ia.post,
                                       fa.post1=todo3_data[[iscena]]$fa.post1, 
                                       fa.post2=todo3_data[[iscena]]$fa.post2, 
                                       fa.post3=todo3_data[[iscena]]$fa.post3, 
                                       post1=todo3_data[[iscena]]$post1, 
                                       post2=todo3_data[[iscena]]$post2, 
                                       ntrial, peff=peff.m[iscena,], theta0, m1 = im1, nsample)
      }
      todo3_a[[iscena]] <- list()
      todo3_a[[iscena]]$pts <- array(0, c(dim(peff.m)[2],6, ntriala))
      todo3_a[[iscena]]$parameter <- matrix(0, nrow = ntriala, ncol = 4)
      todo3_a[[iscena]]$overallpower <- rep(0, ntriala)
      
      for(iia in 1:ntriala){
        todo3_a[[iscena]]$pts[,,iia] <- tempa[[iia]]$pts
        todo3_a[[iscena]]$parameter[iia,] <- tempa[[iia]]$parameter
        todo3_a[[iscena]]$overallpower[iia] <- tempa[[iia]]$overallpower
      }
    }
    
    a_matrix <- cbind(ia = 1:ntriala,
                      a1 = todo3_a[[1]]$parameter[,1], # a1
                      a2 = todo3_a[[1]]$parameter[,2], # a2
                      fwer = todo3_a[[1]]$overallpower, 
                      power = todo3_a[[2]]$overallpower, 
                      ass_scenario1 = colSums(todo3_a[[1]]$pts[,6,]), 
                      ass_scenario2 = colSums(todo3_a[[2]]$pts[,6,]))
    
    a_matrix1 <- a_matrix[which(a_matrix[,4]<target.alpha),] # control fwer in scenario 1
    
    a_matrix2 <- a_matrix1[which( a_matrix1[,5] == max(a_matrix1[,5]) ),] # maximize OP under the alternative scenario
    
    ffr2[[im1]] <- a_matrix2
    
    if(is.vector(ffr2[[im1]])){
      ffr3 <- rbind(ffr3, c(im1, nsample, ffr2[[im1]]))
    }else{
      ffr3 <- rbind(ffr3, cbind(im1, nsample, ffr2[[im1]]))
    }
    
  }
  
  if(is.vector(ffr3)){
    ffr5 <- ffr3
  }else{
    op_star <- max(ffr3[,7]) # maximum overall power
    ffr4 <- ffr3[which(ffr3[,7] >= (op_star-1)),] # allowable overall power
    if(is.vector(ffr4)){
      ffr5 <- ffr4
    }else{
      ffr5 <- ffr4[which(ffr4[,8] == min(ffr4[,8])),] # minimize ASS
    }
  }
  
  if(ffr5[7]>target.op){print(ffr5);break}
}

gamma1 <- log(seq(1,0.5, by = -0.025))/log(0.5)
lambda1 <- seq(0.5,1, by = 0.0025)
a12 <- NULL
for(ia1 in 1:length(lambda1)){ for(ia2 in 1:length(gamma1)){a12 <- rbind(a12,c(lambda1[ia1],gamma1[ia2]))} }

m1_opt <- ffr5[1]
n_opt <- ffr5[2]
a12_opt <- ffr5[3]
gamma_opt <- a12[a12_opt,2]
lambda_opt <- a12[a12_opt,1]
a1_opt <- lambda_opt *(m1_opt/n_opt)^gamma_opt
a2_opt <- lambda_opt

print(paste('optimal (m1, n, gamma, lambda, a1, a2) = (', 
            m1_opt, ', ', n_opt, ', ', gamma_opt, ', ', lambda_opt, ')', sep =''))

save.image("TODO3-n-m1-a1-a2-optimization.RData")
