# using the functions: todo2_data_generation(), todo2_futility()

## simulation scenarios
peff.m <- rbind(c(0.20,0.20), c(0.40,0.40), c(0.40,0.45),
                c(0.40,0.60), c(0.20,0.40), c(0.40,0.70),
                c(0.40,0.43), c(0.45,0.40), c(0.60,0.40))

## number of scenarios
sim.nscenarios <- dim(peff.m)[1]

## null response rate
theta0 <- 0.2

## non-inferior margin
delta1 <- 0.05

## incorrect decisions
idr.m <- matrix(0,nrow=sim.nscenarios,ncol=6)
idr.m[1,1:3] <- 1
idr.m[c(2,3,7:9),c(2,4)] <- 1
idr.m[4:6,c(1,4)] <- 1

## number of simulation trials
ntrial <- 10000

## number of cutoff combinations for futility monitoring
ntriala <- 3674

## number of cutoff combinations for non-inferior comparison
ntrialc <- 4208

## family-wise type I error rate (FWER)
target.alpha <- 10

## size of the inconclusive region size (SIR)
sir <- 15

## maximum rate of selecting an inadequate dose (MRID)
mrid <- 20

## sample size per arm
nsample <- 24
## conducting interim analysis after enrolling n1 patients per arm
n1 <- 12

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.45

ass_n1 <- cbind(9:15, rep(0, length(9:15)))
for(in1 in 9:15){
  todo2_data <- temp <- list()
  for(iscena in 1:4){
    print(paste('Scenario', iscena))
    todo2_data[[iscena]] <- todo2_data_generation(rseed=1, iscena, ntrial, n1, nsample, peff.m, theta0, delta1, sigma1, tau)
  }
  
  todo2_a <- tempa <- list()
  for(iscena in 1:2){
    for(iia in 1:ntriala){
      tempa[[iia]] <- todo2_futility(ia12 = iia, fdata=todo2_data[[iscena]]$data, 
                                     ia.post=todo2_data[[iscena]]$ia.post,
                                     fa.post1=todo2_data[[iscena]]$fa.post1, 
                                     fa.post2=todo2_data[[iscena]]$fa.post2, 
                                     fa.post3=todo2_data[[iscena]]$fa.post3, 
                                     post=todo2_data[[iscena]]$post,
                                     ntrial, peff=peff.m[iscena,], nsample, theta0)
    }
    todo2_a[[iscena]] <- list()
    todo2_a[[iscena]]$pts <- array(0, c(dim(peff.m)[2],6, ntriala))
    todo2_a[[iscena]]$parameter <- matrix(0, nrow = ntriala, ncol = 2)
    todo2_a[[iscena]]$gpower <- rep(0, ntriala)
    
    for(iia in 1:ntriala){
      todo2_a[[iscena]]$pts[,,iia] <- tempa[[iia]]$pts
      todo2_a[[iscena]]$parameter[iia,] <- tempa[[iia]]$parameter
      todo2_a[[iscena]]$gpower[iia] <- tempa[[iia]]$gpower
    }
  }
  
  
  a_matrix <- cbind(ia = 1:ntriala,
                    a1 = todo2_a[[1]]$parameter[,1], 
                    a2 = todo2_a[[1]]$parameter[,2], 
                    fwer = todo2_a[[1]]$gpower, 
                    power = todo2_a[[2]]$gpower, 
                    ass_scenario1 = colSums(todo2_a[[1]]$pts[,6,]), 
                    ass_scenario2 = colSums(todo2_a[[2]]$pts[,6,]))
  
  
  a_matrix1 <- a_matrix[which(a_matrix[,4]<target.alpha),] # control fwer in scenario 1
  a_matrix2 <- a_matrix1[which(a_matrix1[,5]==min(a_matrix1[,5])),] # maximize power in scenario 2
  if(is.vector(a_matrix2)){ia12 <- a_matrix2[1]}else{ # minimize ess in scenario 1
    a12_opt <- a_matrix2[which(a_matrix2[,6]==max(a_matrix2[,6])),1]
  }
  
  a12 <- NULL
  for(ia1 in 1:98){
    for(ia2 in max(51,ia1):99){a12 <- rbind(a12,c(ia1,ia2)/100)} }
  a1_opt <- a12[a12_opt,1];a2_opt <- a12[a12_opt,2]
  
  
  ass_n1[in1/2-8,2] <- colSums(todo2_a[[1]]$pts[,6,a12_opt])
}

ia_opt <- ass_n1[which(ass_n1[,2] == min(ass_n1[,2])),1]


save.image("TODO2_ia_timing_optimization.RData")

