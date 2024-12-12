# using the functions: todo2_data_generation(), todo2_futility(), todo2()
load('E:/projects/20221213-project-TODO/review/20241205/code-github/TODO2-n-m1-lambda-gamma-optimization.RData')
source('E:/projects/20221213-project-TODO/review/20241205/code-github/todo2.R')

## simulation scenarios
peff.m <- rbind(c(0.20,0.20), c(0.40,0.40), 
                c(0.40,0.45), c(0.40,0.60))

## number of scenarios
sim.nscenarios <- dim(peff.m)[1]

## null response rate
theta0 <- 0.2

## non-inferior margin
delta1 <- 0.05

## incorrect decisions
idr.m <- matrix(0,nrow=sim.nscenarios,ncol=6)
idr.m[1,1:3] <- 1
idr.m[c(2,3),c(2,4)] <- 1
idr.m[4,c(1,4)] <- 1

## number of simulation trials
ntrial <- 10000

## number of cutoff combinations for non-inferior comparison
ntrialc <- 4208

## size of the inconclusive region size (SIR)
sir <- 15

## maximum rate of selecting an inadequate dose (MRID)
mrid <- 20

## sample size per arm
nsample <- n_opt
## conducting interim analysis after enrolling m1 patients per arm
m1 <- m1_opt

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.4

# ntrial <- 10
# todo2_data <- temp <- list()
for(iscena in 3:4){
  print(paste('Scenario', iscena))
  todo2_data[[iscena]] <- todo2_data_generation(rseed=1, iscena, ntrial, m1, nsample, peff.m, theta0, delta1, sigma1, tau)
}


todo2_c <- tempc <- list()
for(iscena in 3:4){
  for(iic in 1:ntrialc){
    tempc[[iic]] <- todo2(ic12 =iic, fdata = todo2_data[[iscena]]$data, 
                            ia.post = todo2_data[[iscena]]$ia.post,
                            fa.post1 = todo2_data[[iscena]]$fa.post1, 
                            fa.post2 = todo2_data[[iscena]]$fa.post2, 
                            fa.post3 = todo2_data[[iscena]]$fa.post3, 
                            post = todo2_data[[iscena]]$post,
                            ia.mean = todo2_data[[iscena]]$ia.mean,
                            fa.mean1 = todo2_data[[iscena]]$fa.mean1,
                            fa.mean2 = todo2_data[[iscena]]$fa.mean2,
                            fa.mean3 = todo2_data[[iscena]]$fa.mean3,
                            ntrial, peff=peff.m[iscena,], 
                            theta0, delta1, ia12=a12_opt, idr.m = idr.m[iscena,], wl, m1, nsample)
  }
  todo2_c[[iscena]] <- list()
  todo2_c[[iscena]]$parameter <- matrix(0, nrow = ntrialc, ncol = 6)
  todo2_c[[iscena]]$present <- matrix(0, nrow = ntrialc, ncol = 10)
  
  for(iic in 1:ntrialc){
    todo2_c[[iscena]]$parameter[iic,] <- tempc[[iic]]$parameter
    todo2_c[[iscena]]$present[iic,] <- tempc[[iic]]$present
  }
}


c_matrix <- cbind(ic = 1:ntrialc,
                  c1 = todo2_c[[3]]$parameter[,5], 
                  c2 = todo2_c[[3]]$parameter[,6], 
                  inconclusive_scenario3 = todo2_c[[3]]$present[,3], 
                  inconclusive_scenario4 = todo2_c[[4]]$present[,3], 
                  weighted_loss_scenario3 = todo2_c[[3]]$present[,5], 
                  weighted_loss_scenario4 = todo2_c[[4]]$present[,5], 
                  select_d1_scenario4 = todo2_c[[4]]$present[,1])

c_matrix1 <- c_matrix[which(c_matrix[,4]<sir & c_matrix[,5]<sir & c_matrix[,8]< mrid),] # control sir and mrid
c12_opt <- c_matrix1[which(rowSums(c_matrix1[,6:7]) == min(rowSums(c_matrix1[,6:7]))),1] # maximize weighted loss
c12 <- NULL
for(ic1 in 1:98){
  for(ic2 in ic1:min(ic1+60,99)){c12 <- rbind(c12,c(ic1,ic2)/100)} }
c1_opt <- c12[c12_opt,1];c2_opt <- c12[c12_opt,2]

print(paste('optimal (m1, n, gamma, lambda, c1, c2) = (', 
            m1_opt, ', ', n_opt, ', ', round(gamma_opt,4), ', ', lambda_opt, ', ', c1_opt, ', ', c2_opt,')', sep =''))

save.image("TODO2-c1-c2-optimization.RData")


