# using the functions: todo3_data_generation(), todo3_futility(), todo3()
load('E:/projects/20221213-project-TODO/review/20241119/code/TODO3-n-m1-a1-a2-optimization.RData')
source('E:/projects/20221213-project-TODO/review/20241119/code/todo3.R')

## simulation scenarios
peff.m <- rbind(c(0.20, 0.20, 0.20), c(0.40, 0.40, 0.40), c(0.40, 0.45, 0.45), 
                c(0.40, 0.40, 0.45), c(0.40, 0.40, 0.60), c(0.40, 0.60, 0.60) )

## number of scenarios
sim.nscenarios <- dim(peff.m)[1]

## null response rate
theta0 <- 0.2

## non-inferior margin
delta1 <- 0.05

## incorrect decisions
idr.m <- matrix(0,nrow=sim.nscenarios,ncol=8)
idr.m[1,1:4] <- 1
idr.m[2:sim.nscenarios,c(1:3,5)] <- 1
opt_eff <- function(peff,theta0,delta1){
  opt <- which(peff>theta0 & peff>=(max(peff)-delta1))
  if(length(opt)>0){return(opt[1])}else{return(0)}
}

for(iscena in 2:sim.nscenarios){
  opt_dose_eff <- opt_eff(peff.m[iscena,],theta0,delta1)
  idr.m[iscena,opt_dose_eff] <- 0
}

## number of simulation trials
ntrial <- 10000

## number of cutoff combinations for non-inferior comparison
ntrialc <- 4208

## size of the inconclusive region size (SIR)
sir <- 15

## maximum rate of selecting an inadequate dose (MRID)
mrid <- 25

## sample size per arm
nsample <- n_opt
## conducting interim analysis after enrolling m1 patients per arm
m1 <- m1_opt

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.6

# ## identify optimal cutoffs for futility monitoring
# a12_opt <- ffr5[3]
# gamma1 <- log(seq(1,0.5, by = -0.025))/log(0.5)
# lambda1 <- seq(0.5,1, by = 0.0025)
# a12 <- NULL
# for(ia1 in 1:length(lambda1)){ for(ia2 in 1:length(gamma1)){a12 <- rbind(a12,c(lambda1[ia1],gamma1[ia2]))} }
# a1_opt <- a12[a12_opt,1]*(im1/nsample)^a12[a12_opt,2]
# a2_opt <- a12[a12_opt,1]

## data generation
# todo3_data <- temp <- list()
for(iscena in 3:6){
  print(paste('Scenario', iscena))
  todo3_data[[iscena]] <- todo3_data_generation(rseed=1, iscena, ntrial, m1, nsample, peff.m, theta0, delta1, sigma1, tau)
}



## identify optimal cutoffs for dose comparison
todo3_c <- tempc <- list()
for(iscena in 3:6){
  for(iic in 1:ntrialc){
    tempc[[iic]] <- todo3(ic12 = iic, fdata = todo3_data[[iscena]]$data, 
                          ia.post = todo3_data[[iscena]]$ia.post,
                          fa.post1 = todo3_data[[iscena]]$fa.post1, 
                          fa.post2 = todo3_data[[iscena]]$fa.post2, 
                          fa.post3 = todo3_data[[iscena]]$fa.post3, 
                          ia.mean = todo3_data[[iscena]]$ia.mean,
                          fa.mean1 = todo3_data[[iscena]]$fa.mean1,
                          fa.mean2 = todo3_data[[iscena]]$fa.mean2,
                          fa.mean3 = todo3_data[[iscena]]$fa.mean3,
                          post1 = todo3_data[[iscena]]$post1,
                          post2 = todo3_data[[iscena]]$post2,
                          ntrial, peff=peff.m[iscena,],
                          theta0, delta1, ia12=a12_opt, idr.m = idr.m[iscena,], wl, m1, nsample)
  }
  todo3_c[[iscena]] <- list()
  todo3_c[[iscena]]$parameter <- matrix(0, nrow = ntrialc, ncol = 4)
  todo3_c[[iscena]]$present <- matrix(0, nrow = ntrialc, ncol = 10)
  
  for(iic in 1:ntrialc){
    todo3_c[[iscena]]$parameter[iic,] <- tempc[[iic]]$parameter
    todo3_c[[iscena]]$present[iic,] <- tempc[[iic]]$present
  }
}


c_matrix <- cbind(ic = 1:ntrialc,
                  c1 = todo3_c[[3]]$parameter[,3], 
                  c2 = todo3_c[[3]]$parameter[,4], 
                  inconclusive_scenario3 = todo3_c[[3]]$present[,4], 
                  inconclusive_scenario4 = todo3_c[[4]]$present[,4], 
                  inconclusive_scenario5 = todo3_c[[5]]$present[,4], 
                  inconclusive_scenario6 = todo3_c[[6]]$present[,4], 
                  weighted_loss_scenario3 = todo3_c[[3]]$present[,6], 
                  weighted_loss_scenario4 = todo3_c[[4]]$present[,6], 
                  weighted_loss_scenario5 = todo3_c[[5]]$present[,6], 
                  weighted_loss_scenario6 = todo3_c[[6]]$present[,6], 
                  select_d1_scenario5 = todo3_c[[5]]$present[,1] + todo3_c[[5]]$present[,2],
                  select_d1_scenario6 = todo3_c[[6]]$present[,1])

c_matrix1 <- c_matrix[which(c_matrix[,4]<sir & c_matrix[,5]<sir  & c_matrix[,6]<sir  & c_matrix[,7]<sir & 
                              c_matrix[,12]< mrid & c_matrix[,13]< mrid),] # control sir and mrid
c12_opt <- c_matrix1[which(rowSums(c_matrix1[,8:11])==min(rowSums(c_matrix1[,8:11]))),1] # maximize weighted loss
c12 <- NULL
for(ic1 in 1:98){
  for(ic2 in ic1:min(ic1+60,99)){c12 <- rbind(c12,c(ic1,ic2)/100)} }
c1_opt <- c12[c12_opt,1];c2_opt <- c12[c12_opt,2]

print(paste('optimal (m1, n, gamma, lambda, c1, c2) = (', 
            m1, ', ', nsample, ', ', a12[a12_opt,2], ', ', a12[a12_opt,1], ', ', c1_opt, ', ', c2_opt,')', sep =''))


save.image("TODO3-c1-c2-optimization.RData")

