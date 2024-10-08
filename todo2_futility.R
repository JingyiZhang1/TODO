

# futility monitoring for interim and final analysis
todo2_futility <- function(ia12, fdata, ia.post, fa.post1, fa.post2, fa.post3, post, ntrial, peff, theta0){
  
  SEL <- matrix(0,ncol=length(peff),nrow=ntrial)
  FUTILE <- array(0,c(2,length(peff),ntrial))
  inconclusive <- matrix(1,ncol=length(peff),nrow=ntrial)
  
  ## all posible cutoff combinations for futility monitoring
  a12 <- NULL
  for(ia1 in 1:98){ for(ia2 in max(51,ia1):99){a12 <- rbind(a12,c(ia1,ia2)/100)} }
  a1 <- a12[ia12,1];a2 <- a12[ia12,2]
  
  ## making desicions for the simulated trials generated by function 'todo2_data_generation()'
  for(trial in 1:ntrial){
    
    ## interim analysis
    setI <- 1:2     ## doses passed the interim futility monitoring
    exclude <- NULL ## doses failed to pass the interim futility monitoring
    if(ia.post[trial,2]<a1){ setI <- NULL;exclude <- 1:2
    }else if(ia.post[trial,1]<a1){ setI <- 2;exclude <- 1  }
    
    if(length(exclude)>0){
      inconclusive[trial,exclude] <- -2
      FUTILE[1,exclude,trial] <- FUTILE[2,exclude,trial] <- 1
      fdata[exclude,3:4,trial] <- 0
      fdata[exclude,5:6,trial] <- fdata[exclude,1:2,trial]
    }
    if(length(exclude)==length(peff)){next}
    
    ## final analysis
    setII <- setI    ## doses passed the final futility monitoring
    exclude <- NULL  ## doses failed to pass the final futility monitoring
    
    if(length(setI)==2){
      setII <- which(ia.post[trial,]>=a1 & fa.post1[trial,]>=a2)
      exclude <- which(ia.post[trial,]<a1 | (ia.post[trial,]>=a1 & fa.post1[trial,]<a2))
    }
    
    if(length(setI)==1){
      if(setI==1){
        if(fa.post3[trial,1]>=a2){
          setII <- 1; exclude <- 2
        }else{ setII <- NULL; exclude <- c(1,2) }
      }
      if(setI==2){
        if(fa.post2[trial,2]>=a2){
          setII <- 2; exclude <- 1
        }else{  setII <- NULL;exclude <- c(1,2) }
      }
    }
    
    if(length(exclude)>0){
      inconclusive[trial,exclude] <- -2
      FUTILE[2,exclude,trial] <- 1
    }
    
    if(length(exclude)==length(peff)){next}
  }
  
  pts <- apply(fdata,c(1,2),mean)
  
  futile <- apply(FUTILE,c(1,2), mean)*100
  colnames(futile) <- c('arm1','arm2');rownames(futile) <- c('IA','FA')
  
  none <- mean(colSums(FUTILE[2,,])==2)*100
  gpower <- 0
  for(trial in 1:ntrial){
    if( sum(FUTILE[2,,trial]==(peff<theta0)) ){gpower <- gpower + 100/ntrial}
  }
  
  result.list <- list(peff=peff, pts=round(pts,2), parameter=parameter, gpower=gpower)
  return(result.list)
  
}