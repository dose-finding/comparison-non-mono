library(nnet)
library(dfcrm)
library(pocrm)

bcrm<-function(pI,p.skel,target,cohortsize,ncohort,start.comb,conf.level,mu,scale){
  conjDist <- function(a,p,y,n){
    dist=exp(-((a-mu)^2/(2*scale))) # normal prior distribution
    for(j in 1:length(p)){ # for each probability of toxicity 
      pj=p[j]**exp(a) # model
      dist=dist*pj^y[j]*(1-pj)^(n[j]-y[j])# conjoint distribution
    }
    return(dist)
  }
  
  bMean <- function(a,p,y,n){
    dist=a*exp(-((a-mu)^2/(2*scale))) # a*normal prior distribution
    for(j in 1:length(p)){ # for each probability of toxicity 
      pj=p[j]**exp(a) # model
      dist=dist*pj^y[j]*(1-pj)^(n[j]-y[j]) # marginal distribution
    }
    return(dist)
  }
  
  bVar <- function(a,p,y,n){
    dist=a^2*exp(-((a-mu)^2/(2*scale))) # a?*normal prior distribution
    for(j in 1:length(p)){ # for each probability of toxicity 
      pj=p[j]**exp(a) # model
      dist=dist*pj^y[j]*(1-pj)^(n[j]-y[j]) # marginal distribution
    }
    return(dist)
  }
  
  ### run a trial 	
  # generate data for a new cohort of patients

  ncomb = length(p.skel) ; 
  y = n = rep(0,ncomb)
  comb.select = rep(0, ncomb)
  ptox.hat = numeric(ncomb) ; 
  comb.curr = start.comb
  
  stop = 0 ;
  k=1
  while (k <= ncohort){
    # Generate a new set of data : n patients and y indicate DLT or not    
    y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,pI[comb.curr]);
    n[comb.curr] = n[comb.curr] + cohortsize;
    
    marginal.tox = integrate(conjDist, lower = -Inf, upper = Inf, p = p.skel, n = n, y = y, abs.tol = 0)$value
    est.tox = integrate(bMean, lower = -10, upper = 10, p = p.skel, n = n, y = y, abs.tol = 0)$value/marginal.tox
    e2.tox = integrate(bVar, lower = -10, upper = 10, p = p.skel, n = n, y = y, abs.tol = 0)$value/marginal.tox
    
    ptox.hat = p.skel**exp(est.tox) # Prob of toxicity for selected order
    post.var = e2.tox-(est.tox)^2
    crit.tox=qnorm(0.5+conf.level/2)
    lb.tox=est.tox-crit.tox*sqrt(post.var)
    ub.tox=est.tox+crit.tox*sqrt(post.var)
    ptoxL=p.skel**exp(ub.tox)
    ptoxU=p.skel**exp(lb.tox)
    
    
    if (ptoxL[1] > target){
      stop = 1
      break
    }
    
    distance=abs(ptox.hat-target)
    comb.best<-which.is.max(-distance)
    comb.curr<-min(comb.best,comb.curr+1)
                   
    k = k+1
  }
  
  if(stop==0){
    comb.select[comb.curr]=comb.select[comb.curr]+1
  }
  
  return(list(comb.select=comb.select,tox.data=y,pt.allocation=n,var_post = post.var
              , Stop=stop))
}


bcrm.sim<-function(pI,p.skel,target,cohortsize,ncohort,start.comb,conf.level,mu,scale,ntrial){
#  pb <- winProgressBar(title = "progress bar", min = 0, max = ntrial, width = 300)
  
  ncomb = length(pI) ; 
  stop <- sel_tox <- n_tox <- matrix(nrow=ntrial, ncol=1) ;
  comb.select <- y <- n <- matrix(nrow=ntrial, ncol=ncomb)
  nstop = 0 ;
  tox_dl <- which(pI > target)
  
  for(j in 1:ntrial){
#    setWinProgressBar(pb, j, title=paste( round(j/ntrial*100, 0),"% done"))
    result<-bcrm(pI,p.skel,target,cohortsize,ncohort,start.comb,conf.level,mu,scale)
    comb.select[j,]=result$comb.select
    y[j, ] = result$tox.data
    n[j, ] = result$pt.allocation
    stop[j, ] <- result$Stop
    sel_tox[j, ] <- sum(result$comb.select[tox_dl])
    n_tox[j, ] <- sum(result$pt.allocation[tox_dl])
    nstop = nstop+result$Stop
  }
  
#close(pb)
trial = ntrial-nstop
resultats <- list()
resultats[[1]] <- round(pI,3);
resultats[[2]] <- formatC( 100 * colMeans(comb.select), digits=1, format="f");
resultats[[3]] <- formatC( colMeans(n), digits=1, format="f");
resultats[[4]] <- formatC( 100 * colMeans(stop), digits=1, format="f");
resultats[[5]] <- formatC( 100 * colMeans(sel_tox), digits=1, format="f");
resultats[[6]] <- formatC( 100 * colMeans(n_tox) / sum(colMeans(n)), digits=1, format="f")

names(resultats) <- c("True tox probability",
                        "Selection percentage",
                        "Average number of patients assigned",
                        "Proportion of early stopping for safety",
                        "Probability of selecting toxic doses",
                        "Percentage of patients assigned to toxic doses")
  return (resultats)
}
