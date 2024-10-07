setwd('~/ZJY/rce')

## simulation scenarios
peff.m <- rbind(c(0.20, 0.20, 0.20), c(0.40, 0.40, 0.40), c(0.40, 0.45, 0.45), 
                c(0.40, 0.40, 0.45), c(0.40, 0.40, 0.60), c(0.40, 0.60, 0.60), 
                c(0.20, 0.45, 0.45), c(0.15, 0.45, 0.45), c(0.40, 0.40, 0.70) )

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

## family-wise type I error rate (FWER)
target.alpha <- 10

## size of the inconclusive region size (SIR)
sir <- 15

## maximum rate of selecting an inadequate dose (MRID)
mrid <- 25

## sample size per arm
nsample <- 24
## conducting interim analysis after enrolling n1 patients per arm
n1 <- 12

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.7

## optimized design parameters
a12_opt <- 2047
c12_opt <- 2784

## data generation
todo3_data <- temp <- list()
for(iscena in 1:sim.nscenarios){
  print(paste('Scenario', iscena))
  todo3_data[[iscena]] <- todo3_data_generation(rseed=1, iscena, ntrial, n1, nsample, peff.m, theta0, delta1, sigma1, tau)
}

re <- list()
main_result <- bias_mse <- NULL
for(iscena in 1:sim.nscenarios){
  re[[iscena]] <- todo3(ic12 = c12_opt,
                        fdata = todo3_data[[iscena]]$data, 
                        ia.post = todo3_data[[iscena]]$ia.post,
                        fa.post1 = todo3_data[[iscena]]$fa.post1, 
                        fa.post2 = todo3_data[[iscena]]$fa.post2, 
                        fa.post3 = todo3_data[[iscena]]$fa.post3,
                        post1 = todo3_data[[iscena]]$post1,
                        post2 = todo3_data[[iscena]]$post2,
                        ia.mean = todo3_data[[iscena]]$ia.mean, 
                        fa.mean1 = todo3_data[[iscena]]$fa.mean1, 
                        fa.mean2 = todo3_data[[iscena]]$fa.mean2, 
                        fa.mean3 = todo3_data[[iscena]]$fa.mean3, 
                        ntrial, peff = peff.m[iscena,], theta0, delta1, ia12 = a12_opt, idr.m = idr.m[iscena,], wl)
  
  main_result <- rbind(main_result, c(re[[iscena]]$present[1:9], 
                                      go1 = 100-re[[iscena]]$futile[2,1],
                                      go2 = 100-re[[iscena]]$futile[2,2], 
                                      go2 = 100-re[[iscena]]$futile[2,3], 
                                      re[[iscena]]$present[10]))
  
  bias_mse <- rbind(bias_mse, c(bias = mean(re[[iscena]]$bias), mse = mean(re[[iscena]]$mse)))
}
colnames(main_result) <- c('sel1','sel2','sel3','SIR','IDR','WL','ASS1','ASS2','ASS3','go1','go2','go3','none')
round(main_result, 2)
round(bias_mse, 2)

save.image("TODO3.RData")



