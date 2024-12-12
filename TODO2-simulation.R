# using the functions: todo2_data_generation(), todo2()
load('E:/projects/20221213-project-TODO/review/20241205/code-github/TODO2-c1-c2-optimization.RData')

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

## sample size per arm
nsample <- n_opt
## conducting interim analysis after enrolling n1 patients per arm
n1 <- m1_opt

## hyperparameters
sigma1 <- 3
tau <- 1

## discount factor of an inconclusive decision
wl <- 0.4

todo2_data <- temp <- list()
for(iscena in 1:sim.nscenarios){
  print(paste('Scenario', iscena))
  todo2_data[[iscena]] <- todo2_data_generation(rseed=1, iscena, ntrial, m1, nsample, peff.m, theta0, delta1, sigma1, tau)
}

re <- list()
main_result <- bias_mse <- NULL
for(iscena in 1:sim.nscenarios){
  re[[iscena]] <- todo2(ic12 = c12_opt,
                        fdata = todo2_data[[iscena]]$data, 
                        ia.post = todo2_data[[iscena]]$ia.post,
                        fa.post1 = todo2_data[[iscena]]$fa.post1, 
                        fa.post2 = todo2_data[[iscena]]$fa.post2, 
                        fa.post3 = todo2_data[[iscena]]$fa.post3,
                        post = todo2_data[[iscena]]$post,
                        ia.mean = todo2_data[[iscena]]$ia.mean, 
                        fa.mean1 = todo2_data[[iscena]]$fa.mean1, 
                        fa.mean2 = todo2_data[[iscena]]$fa.mean2, 
                        fa.mean3 = todo2_data[[iscena]]$fa.mean3, 
                        ntrial, peff = peff.m[iscena,], theta0, delta1, ia12 = a12_opt, idr.m = idr.m[iscena,], wl, m1, nsample)
  
  main_result <- rbind(main_result, re[[iscena]]$present)
  
  bias_mse <- rbind(bias_mse, c(bias = mean(re[[iscena]]$bias[1], re[[iscena]]$bias[2]),
                                mse = mean(re[[iscena]]$mse[1], re[[iscena]]$mse[2])))
}
round(main_result, 2)
round(bias_mse, 2)

save.image("TODO2.RData")

