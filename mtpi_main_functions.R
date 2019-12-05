#### phase1 designs comparison: simulation study 06May2019
#### This is the R code for performing simulation studies based on the mTPI method (Ji et al, 2009). 
#### The users only need to specify an equivalence interval (pT-epsi1, pT+epsi2)
#### where pT is the target tox probility for the MTD and epsi1 and epsi2 are two small values so that any dose in the equivalence interval will be considered potential candidates for the true MTD. 
# The design is then upon specification of 
# p         | vector of true toxicity probability
# D         | the number of dose levels
# pT        | target tox probility is pT
# startdose | the starting dose level 
# eps1;eps2 | (pT-epsi1, pT+epsi2) is the equivalence interval 
# xi        | the cutoff probability for excessive toxicity.
# sampsize  | the max sample size is sampsize
# csize     | cohorts size
# simN      | the number of simulated trials #<-1000  ## 1000 simulations

############################################### Functions needed ####################################
## pava is the pool adjacent violator algorithm to perform isotonic transformation for the posterior means later
pava <- function (x, wt = rep(1, length(x))) {
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

## betavar computes variances of beta distributions 

betavar<-function(a,b){
  resp <- a*b/((a+b)^2*(a+b+1))
  return(resp)
}

####################################### Start simulations ###############
mtpi.sim <- function (p,D,pT,startdose,eps1,eps2,xi,sampsize,csize,simN){
a<-1; b<-1; ## the prior is beta(1,1). Keep a, b fixed at 1 -- essential for this design to work.
datan<-matrix(rep(0,simN*D),ncol=D)
datax<-matrix(rep(0,simN*D),ncol=D)

rez<-rep(0,simN)

for(sim in 1:simN){
  
  x<-rep(0,D); n<-rep(0,D)
  pa<-rep(0, D); pb<-rep(0,D)
  q <- rep(0,3)
  d<-startdose; st<-0; nodose<-0; maxdose<-1; toxdose<-D+1; seldose<-0
  
  
  while(st==0){  ## st = 1 indicates the trial must be terminated
    maxdose<-max(maxdose, d)
    
    ### generate random toxicity response
    xx <- 0
    for(i in 1:csize){
      ttt <- runif(1)
      if(ttt < p[d]) xx <- xx + 1
    }
    
    x[d] <- x[d] + xx; n[d] <- n[d] + csize
    
    #### Update posterior beta distribution
    pa[d]<-x[d]+a; pb[d]<-n[d]-x[d]+b
    pa[d+1]<-x[d+1]+a; pb[d+1]<-n[d+1]-x[d+1]+b
    
    ###Compute the indicator T_{i+1} to see if the next dose is too toxic
    if(d<D){
      temp<-1-pbeta(pT, pa[d+1], pb[d+1])
      if(temp>xi) {tt<-1; toxdose<-d+1} else {tt<-0}
    }

    ##Compute the UPM for three intervals defined by the equivalence interval      
    q[1] <- (1-pbeta(eps2+pT, pa[d], pb[d]))/(1-eps2-pT)
    q[2] <- (pbeta(eps2+pT, pa[d], pb[d]) - pbeta(pT-eps1, pa[d], pb[d]))/(eps2+eps1)
    q[3] <- (pbeta(pT-eps1, pa[d], pb[d])/(pT-eps1))*(1-tt)

    ### implement the dose-assignment rules based on the UPM
    if(d==1){
      ## if the first dose is too toxic, the trial will be terminated
      if((1-pbeta(pT, pa[d], pb[d]))>xi){st=1; nodose<-1;} 
      else{
        if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
        if((q[3]>q[1])&&(q[3]>q[2])){d<-d+1}
      }
    }
    else{
      if(d==D){
        if((q[1]>q[2])&&(q[1]>q[3])){d<-d-1}
        if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
      }
      else{
        if((d>1)&&(d<D)){
          if((q[1]>q[2])&&(q[1]>q[3])){d<-d-1}
          if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
          if((q[3]>q[1])&&(q[3]>q[2])){d<-d+1}
        }
      }
    }
    total<-sum(n)
    if(total >= sampsize){st<-1}
  }

  ### compute the posterior mean from the beta distribution
  if(nodose==0){
    tdose<-min(maxdose, toxdose-1)
    pp<-rep(-100,tdose)
    pp.var<-rep(0, tdose)
    
    for(i in 1:tdose){
      pp[i] <- (x[i]+.005)/(n[i]+.01); pp.var[i] <- betavar(x[i]+0.005, n[i]-x[i]+0.005) ### here adding 0.005 is to use beta(0.005, 0.005) for estimating the MTD, which is different from the dose-finding.
    }
    
    pp<-pava(pp, wt=1/pp.var)  ## perform the isotonic transformation using PAVA with weights being posterior variances
    
    for(i in 2:tdose){
      pp[i] <- pp[i] + i*1E-10 ## by adding an increasingly small number to tox prob at higher doses, it will break the ties and make the lower dose level the MTD if the ties were larger than pT or make the higher dose level the MTD if the ties are smaller than pT
    }
    seldose<-sort(abs(pp-pT), index.return=T)$ix[1]
    ##seldose is the final MTD that is selected based on the order-transofromed posterior means
  }
  rez[sim] <- seldose;
  for(i in 1:D){
    datan[sim,i] <- n[i]
    datax[sim,i] <- x[i]
  }
}
##rez[is.na(rez)]<-0
dtox <- which(p > pT)
mdatan<-aaa<-rep(0,D)
rep(0,length(which(p > pT)))
################## output results ################################
for(i in 1:D){
  aaa[i] <- sum(rez==i)/simN ### aaa is the propotion of be selected as the MTD
}
for(i in which(p > pT) ){
  mdatan[i] <- mean(datan[,i]) # count mean n patients assigned to toxic
}

resultats <- list()
resultats[[1]] <- round(p,3);
resultats[[2]] <- formatC( aaa*100, digits=1, format="f");
resultats[[3]] <- formatC( colMeans(datan), digits=1, format="f");
resultats[[4]] <- formatC( 100*(1-sum( aaa )), digits=1, format="f");
resultats[[5]] <- formatC( 100 * sum(aaa[dtox] ), digits=1, format="f");
resultats[[6]] <- formatC( 100 * sum( mdatan)/sum( colMeans(datan) ), digits=1, format="f")
names(resultats) <- c("True tox probability",
                      "Selection percentage",
                      "Average number of patients assigned",
                      "Proportion of early stopping for safety",
                      "Probability of selecting toxic doses",
                      "Percentage of patients assigned to toxic doses")
return (resultats)
}