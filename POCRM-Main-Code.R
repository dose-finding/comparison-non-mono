# This code can be used to upload the function for the PO-CRM design
# to reproduce the results of the paper
#
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2019)
#
#
#    Please use Code POCRM_simulations.R to produce the operating characteristics


library(nnet)
library(dfcrm)
library(pocrm)


###Load the function 'bpocrm' 
bpocrm<-function(pI,p.skel,prior.o,target,cohortsize,ncohort,start.comb,conf.level,mu,scale){
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

  ncomb = ncol(p.skel) ; 
  y = n = rep(0,ncomb)
  comb.select = rep(0, ncomb)
  ptox.hat = numeric(ncomb) ; 
  marginal.tox = est.tox = e2.tox = c(rep(NA, nrow(p.skel)))
  comb.curr = start.comb
  
  stop = 0 ; 
  k=1
  while (k <= ncohort){
    # Generate a new set of data : n patients and y indicate DLT or not    
    y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,pI[comb.curr]);
    n[comb.curr] = n[comb.curr] + cohortsize;
    
    for (j in 1:nrow(p.skel)){
      marginal.tox[j] = integrate(conjDist, lower = -Inf, upper = Inf, p = p.skel[j,], n = n, y = y, abs.tol = 0)$value
      est.tox[j] = integrate(bMean, lower = -10, upper = 10, p = p.skel[j,], n = n, y = y, abs.tol = 0)$value/marginal.tox[j]
      e2.tox[j] = integrate(bVar, lower = -10, upper = 10, p = p.skel[j,], n = n, y = y, abs.tol = 0)$value/marginal.tox[j]
    }
    
    ptox.post = (marginal.tox*prior.o) / sum(marginal.tox*prior.o) # Posterior probabilities
    mtox.sel = which.is.max(ptox.post) # Select order which has higher posterior toxicity prob
    
    ptox.hat = p.skel[mtox.sel,]**exp(est.tox[mtox.sel]) # Prob of toxicity for selected order
    post.var = e2.tox[mtox.sel]-(est.tox[mtox.sel])^2
    crit.tox=qnorm(0.5+conf.level/2)
    lb.tox=est.tox[mtox.sel]-crit.tox*sqrt(post.var)
    ub.tox=est.tox[mtox.sel]+crit.tox*sqrt(post.var)
    ptoxL=p.skel[mtox.sel,]**exp(ub.tox)
    ptoxU=p.skel[mtox.sel,]**exp(lb.tox)
    
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
  return(list(ord_sel=mtox.sel, comb.select=comb.select,tox.data=y,pt.allocation=n,var_post = post.var
              , Stop=stop))
}
#'bpocrm' end here

###Load the function 'bpocrm.sim' 
bpocrm.sim<-function(pI,dose,p.skel,prior.o,target,cohortsize,ncohort,start.comb,conf.level,mu,scale,ntrial){
#  pb <- winProgressBar(title = "progress bar", min = 0, max = ntrial, width = 300)
  
  ncomb = length(pI) ; order =  matrix(nrow=ntrial, ncol=1) 
  stop <- sel_tox <- n_tox <- matrix(nrow=ntrial, ncol=1) 
  comb.select <- y <- n <- matrix(nrow=ntrial, ncol=ncomb)
  tox_dl <- which(pI > target)
  for(j in 1:ntrial){
#    setWinProgressBar(pb, j, title=paste( round(j/ntrial*100, 0),"% done"))
    result<-bpocrm(pI,p.skel,prior.o,target,cohortsize,ncohort,start.comb,conf.level,mu,scale)
    comb.select[j,]=result$comb.select
    y[j,]=result$tox.data
    n[j,]=result$pt.allocation
    order[j,] <- result$ord_sel
    stop[j,] <- result$Stop
    sel_tox[j, ] <- sum(result$comb.select[tox_dl])
    n_tox[j, ] <- sum(result$pt.allocation[tox_dl])
  }
#  close(pb)
  
  resultats <- list()
  resultats[[1]] <- round(pI,3)
  resultats[[2]] <- formatC( 100 * colMeans(comb.select), digits=1, format="f")
  resultats[[3]] <- formatC( colMeans(n), digits=1, format="f")
  resultats[[4]] <- formatC( 100 * colMeans(stop), digits=1, format="f")
  resultats[[5]] <- formatC( 100 * colMeans(sel_tox), digits=1, format="f")
  resultats[[6]] <- formatC( 100 * colMeans(n_tox) / sum(colMeans(n)), digits=1, format="f")
  
  names(resultats) <- c("True tox probability",
                        "Selection percentage",
                        "Average number of patients assigned",
                        "Proportion of early stopping for safety",
                        "Probability of selecting toxic doses",
                        "Percentage of patients assigned to toxic doses")

  return (resultats)
}
report <- function(sim_list,lab=" ") {
  anpa <- paste0("(", rbind(sim_list[[3]]), ")")
  res <- cbind(c("scenario",lab,""),rbind(do.call(rbind,sim_list[c(1,2)]),anpa), c("",sim_list[4],""),c("",sim_list[5],""),c("", sim_list[6],""))
  colnames(res) <- c("", paste("d",1:6,sep=""), "Stop" , "SelTox" , "%Tox")
  return(res)}