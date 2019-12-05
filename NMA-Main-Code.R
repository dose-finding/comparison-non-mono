# This code can be used to upload the function for the NMA design
# to reproduce the results of the paper
#
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2019)
#
#     for the NMA design.
#
#    Please use Code NMA-Result.R to produce the operating characteristics



median.beta<-function(x,n){y<-x/n
  return(y)}

loss.uni<-function(p.tox,target.tox,n=1,kappa=0.5,form,pavel=1){
  if (form=="wde"){
    y1<-(p.tox-target.tox)^2
    y2<-(p.tox^pavel)*((1-p.tox)^(2-pavel))
    y<-(y1/y2)*n^(2*kappa-1)
    return(y)
  }}

which.second.max<-function(x){n<-length(x)
  y<-which(x==sort(x,partial=n-1)[n-1])
  return(y)}

which.third.max<-function(x){n<-length(x)
  y<-which(x==sort(x,partial=n-2)[n-2])
  return(y)} 
 
prob.loss<-function(x){h<-length(x)
  y1<-1/x
  if (any(y1==Inf)){
    Q<-which(y1==Inf)
    y1[Q]<-100}
  y2<-mat.or.vec(h,1)
  y2[which.max(y1)]<-max(y1)
  y2[which.second.max(y1)]<-y1[which.second.max(y1)]
  y<-y2/sum(y2)
  return(y)}

nma.design<-function(P,target,n,cohort,start.dose=1,assigment="cohort-by-cohort",only.first=TRUE,Pr=NULL,initial.prior.value=target,step=NULL,b=1,kappa=0.5,final="plugin",safety.constraint=T,rate=NULL,threshold=target,final.prob=NULL,minimum.sample=T,minimum=1,nsims,trajectories=FALSE,form="wde",pavel=1,futility=FALSE,threshold.low=NULL,rate2=100,final.prob.low=NULL,lambda=0.25){
  final.probability2<-final.prob.low
  down<-1/lambda
  final.probability<-final.prob
  c.star<-threshold
  c.star2<-threshold.low
  
  if (is.null(Pr) & is.null(step) & is.null(initial.prior.value)){
    print("Error! Vector of prior probabilities or the initial value and step between doses should be specified")
    stop()
  }
  
  if (is.null(P)){
    print ("Vector of true probabilities is not specified.")
    stop()
  }
  
  if (is.null(target)){
    print ("Target value is not specified.")
    stop()
  }
  
  if (is.null(cohort)){
    print ("Cohort size is not specified.")
    stop()
  }
  
  if (is.null(assigment)){
    print ("Assigment rule is not specified. Please choose 'cohort-by-cohort` or `randomization`")
    stop()
  }
  
  
  if (is.null(b)){
    print ("Parameters of Beta prior are not specified.")
    stop()
  }
  
  if (is.null(kappa)){
    kappa<-0.5
  }
  
  if (is.null(final)){
    final<-"plugin"
  }
  
  if (is.null(safety.constraint)){
    rate<-NULL
    c.star<-NULL
  } else{
    if (safety.constraint & is.null(rate) & is.null(c.star)){
      print("Safety constaint is not fully specified.")
      stop()
    }
  }
  
  
  if (is.null(minimum.sample)){
    minimum<-NULL
  } else{
    if (minimum.sample & is.null(minimum)){
      print("Minimum sample size is not specified.")
      stop()
    }
  }
  
  
  if (is.null(nsims)){
    print ("Number of simulations is not specified.")
    stop()
  }
  
  # begin of the procedure
  
  sim<-nsims
  
  if (assigment=="cohort-by-cohort"){
    
    
    
    N<-round(n/cohort)                                        
    M<-length(P)
    
    if (is.null(Pr)){    
      Pr<-c(1:M)*(step/(M-1))+(initial.prior.value-step/(M-1))
    } else {}        
    
    if(length(Pr)!=length(P)){
      print("Length of prior vector and true probabilities vector are different")
      cat("Length(True Probabilities)=", length(P))
      cat(" and ")
      cat("Length(Prior)=", length(Pr))
      stop()
    }
    
    
    
    theta<-1
    prob<-mat.or.vec(N+1,M)
    probability2<-mat.or.vec(N+1,M)
    losses<-mat.or.vec(N+1,M)
    result.exp<-mat.or.vec(sim,M)
    rec.all<-mat.or.vec(sim,M)
    result.suc<-mat.or.vec(sim,1)
    bias.final<-mat.or.vec(sim,1)
    several2<-mat.or.vec(sim,1)
    several3<-mat.or.vec(sim,1)
    exp<-mat.or.vec(N+1,M)
    exp.cum<-mat.or.vec(N+1,M)
    suc<-mat.or.vec(N+1,M)
    suc.cum<-mat.or.vec(N+1,M)
    experiment<-array(0,dim=c(N+1,M,sim))
    counter<-0
    
    
    for (z in 1:sim){
      prob<-mat.or.vec(N+1,M)
      losses<-mat.or.vec(N+1,M)
      
      theta<-1
      f.probability<-0.875
      exp[1,]<- exp.cum[1,]<-b
      suc[1,]<-suc.cum[1,]<-Pr*b
      
      
      if(safety.constraint){
        for (v in 1:M){
          prob[1,v]<-pbeta(c.star,suc.cum[1,v]+1,exp.cum[1,v]-suc.cum[1,v]+1, lower.tail = FALSE)
          check<-max(1-rate*exp.cum[1,v],final.probability)
          if (prob[1,v]>check){
            losses[1,v]<-Inf
          }
        }
        if(losses[1,1]==Inf){
          losses[1,]<-Inf
        }
        if(losses[1,2]==Inf){
          losses[1,2:M]<-Inf
        }
        if(losses[1,3]==Inf){
          losses[1,5:6]<-Inf
        }
        if(losses[1,4]==Inf){
          losses[1,6]<-Inf
        }
        if(losses[1,5]==Inf){
          losses[1,6]<-Inf
        }
      }
      
      
      p.est<-median.beta(suc.cum[1,],exp.cum[1,])
      losses[1,]<-loss.uni(p.est,target.tox=target,n=1,kappa,form,pavel)
      
      
      if (all(losses[1,]==Inf))    {print("terminated initially")
        result.exp[z,]<-NA
        result.suc[z]<-NA
        counter<-counter+1
        break} 
      
      else {
        nextdose<-start.dose
        if (length(nextdose)!=1){nextdose<-min(nextdose)}
        exp[2,]<-0
        exp[2,nextdose]<-cohort
        exp.cum[2,]<-exp.cum[1,]
        exp.cum[2,nextdose]<-exp.cum[1,nextdose]+exp[2,nextdose]
        suc[2,]<-0
        response<-rbinom(cohort,1,P[nextdose])
        suc[2,nextdose]<-sum(response)
        suc.cum[2,]<-suc.cum[1,] 
        suc.cum[2,nextdose]<-suc.cum[1,nextdose]+suc[2,nextdose]
      }
      
      # cat("first dose is",nextdose,"\n")
      
      ###########
      j<-2
      ###########
      while (j<N+1){
        losses[j,]<-losses[j-1,]
        p.est<-median.beta(suc.cum[j,]-suc.cum[1,]+suc.cum[1,]/(sum(exp.cum[j,]))^(1/down),exp.cum[j,]-exp.cum[1,]+exp.cum[1,]/(sum(exp.cum[j,]))^(1/down))
        losses[j,]<-loss.uni(p.est,target.tox=target,n=exp.cum[j,],kappa,form,pavel)
        
        if(safety.constraint){
          for (v in 1:M){
            prob[j,v]<-pbeta(c.star,suc.cum[j,v]+1,exp.cum[j,v]-suc.cum[j,v]+1, lower.tail = FALSE)
            check<-max(1-rate*exp.cum[j,v],final.probability)
            if(v==1){
              check<-max(1-rate*exp.cum[j,v],f.probability)
            }
            if (prob[j,v]>check){
              losses[j,v]<-Inf
            }
          }
          if(losses[j,1]==Inf){
            losses[j,]<-Inf
          }
          if(losses[j,2]==Inf){
            losses[j,2:M]<-Inf
          }
          if(losses[j,3]==Inf){
            losses[j,5:6]<-Inf
          }
          if(losses[j,4]==Inf){
            losses[j,6]<-Inf
          }
          if(losses[j,5]==Inf){
            losses[j,6]<-Inf
          }
        }
        
        
        if(futility){
          for (v in 2:M){
            probability2[j,v]<-pbeta(c.star2, suc.cum[j,v]+1,exp.cum[j,v]-suc.cum[j,v]+1, lower.tail = FALSE)
            check<-min(rate2*exp.cum[j,v],final.probability2)
            if (probability2[j,v]<check){
              losses[j,v]<-99
            }
          }
          if(losses[j,1]==99){
            losses[j,1]<-Inf
          }
          if(losses[j,2]==99){
            losses[j,1]<-Inf
          }
          if(losses[j,3]==99){
            losses[j,1:2]<-Inf
          }
          if(losses[j,4]==99){
            losses[j,1:2]<-Inf
          }
          if(losses[j,5]==99){
            losses[j,1:3]<-Inf
          }
          if(losses[j,6]==99){
            losses[j,1:5]<-Inf
          }
        }
        
        prevdose<-nextdose
        loss.cand<-losses[j,]
        
        if (sum(response)>0){   
          if (prevdose==1){
            loss.cand[2:6]<-Inf         
          }
          if(prevdose==2){
            loss.cand[3:6]<-Inf 
          }
          if(prevdose==3){
            loss.cand[5:6]<-Inf
          }
          if(prevdose==4 | prevdose==5){
            loss.cand[6]<-Inf
          }
        }
        
        if (sum(response)==0){
          if(prevdose==1){
            loss.cand[3:6]<-Inf  
          }
          if(prevdose==2){
            loss.cand[1]<-Inf
            loss.cand[5:6]<-Inf  
          }
          if(prevdose==3){
            loss.cand[1:2]<-Inf
            loss.cand[6]<-Inf  
          }
          if(prevdose==4){
            loss.cand[1:2]<-Inf
          }
          if(prevdose==5){
            loss.cand[1:3]<-Inf
          }
          if(prevdose==6){
            loss.cand[1:5]<-Inf
          }
        }
        
        if (all(losses[j,]==Inf)){
          break 
        }  
        else {    
          nextdose<-which.min(loss.cand)
          
          exp[j+1,]<-suc[j+1,]<-0
          exp[j+1,nextdose]<-cohort
          exp.cum[j+1,]<-exp.cum[j,]
          exp.cum[j+1,nextdose]<-exp.cum[j,nextdose]+exp[j+1,nextdose]
          response<-rbinom(cohort,1,P[nextdose])
          suc[j+1,nextdose]<-sum(response)
          suc.cum[j+1,]<-suc.cum[j,] 
          suc.cum[j+1,nextdose]<-suc.cum[j,nextdose]+suc[j+1,nextdose]
        }
        
        j<-j+1
      }
      result.exp[z,]<-(exp.cum[N+1,]-b)
      result.suc[z]<-sum(suc.cum[N+1,]-suc.cum[1,])
      
      
      if (all(losses[j,]==Inf)){
        rec.all[z,]<-0
        counter<-counter+1
        result.exp[z,]<-(exp.cum[j,]-b)
        result.suc[z]<-sum(suc.cum[j,]-suc.cum[1,])
      }else{

        prob[N+1,]<-prob[N,]
        
        j<-N+1

        p.est<-median.beta(suc.cum[j,]-suc.cum[1,]+suc.cum[1,]/(sum(exp.cum[j,]))^(1/down),exp.cum[j,]-exp.cum[1,]+exp.cum[1,]/(sum(exp.cum[j,]))^(1/down))
        losses[j,]<-loss.uni(p.est,target.tox=target,n=exp.cum[j,],kappa,form,pavel)
        
        
        for (u in 1:M){
          if(is.na(losses[j,u])){
            losses[j,u]<-Inf
          }
        }
        
        nextdose<-which.min(losses[j,])
        
        if(safety.constraint){
          for (v in 1:M){
            prob[j,v]<-pbeta(c.star,suc.cum[j,v]+1,exp.cum[j,v]-suc.cum[j,v]+1, lower.tail = FALSE)
            check<-final.probability
            if(v==1){check<-f.probability}
            if (prob[j,v]>check){
              losses[j,v]<-Inf
            }
          }
          if(losses[j,1]==Inf){
            losses[j,]<-Inf
          }
          if(losses[j,2]==Inf){
            losses[j,2:M]<-Inf
          }
          if(losses[j,3]==Inf){
            losses[j,5:6]<-Inf
          }
          if(losses[j,4]==Inf){
            losses[j,6]<-Inf
          }
          if(losses[j,5]==Inf){
            losses[j,6]<-Inf
          }
        }
        
        
        if(futility){
          for (v in 2:M){
            probability2[j,v]<-pbeta(c.star2, suc.cum[j,v]+1,exp.cum[j,v]-suc.cum[j,v]+1, lower.tail = FALSE)
            check<-final.probability2
            if (probability2[j,v]<check){
              losses[j,v]<-99
            }
          }
          if(losses[j,1]==99){
            losses[j,1]<-Inf
          }
          if(losses[j,2]==99){
            losses[j,1]<-Inf
          }
          if(losses[j,3]==99){
            losses[j,1:2]<-Inf
          }
          if(losses[j,4]==99){
            losses[j,1:2]<-Inf
          }
          if(losses[j,5]==99){
            losses[j,1:3]<-Inf
          }
          if(losses[j,6]==99){
            losses[j,1:5]<-Inf
          }
        }
        
        if(minimum.sample){
          enough<-exp.cum[j,]-exp.cum[1,]>=minimum
          for (i in 1:M){
            if (!enough[i]){
              losses[j,i]<-Inf
            }
            else{}
          }
        }
        
        
        
        if (all(losses[j,]==Inf)){
          counter<-counter+1
          rec.all[z,]<-0
          result.exp[z,]<-(exp.cum[j,]-b)
          result.suc[z]<-sum(suc.cum[j,]-suc.cum[1,])
        }
        
        else {
          first<-which.min(losses[j,])
          second<-which.second.max(-losses[j,])[1]
          third<-which.third.max(-losses[j,])[1]
          set2<-cbind(first,second)
          set3<-cbind(first,second,third)
          
          X<-abs(P-target)
          target.doses<-which(X == min(X))
          if(target.doses %in% set2){
            several2[z]<-1
          }else{
            several2[z]<-0
          }
          
          if(target.doses %in% set3){
            several3[z]<-1
          }else{
            several3[z]<-0
          }
          
          nextdose<-which.min(losses[j,])
          if (length(nextdose)==1){}
          else {nextdose<-min(nextdose)}
          rec.all[z,nextdose]<-1
          bias.final[z]<-(median.beta(suc.cum[j,1]-suc.cum[1,1],exp.cum[j,1]-exp.cum[1,1])-median.beta(suc.cum[j,2]-suc.cum[1,2],exp.cum[j,2]-exp.cum[1,2]))-(P[1]-P[2])
          experiment[,,z]<-exp
        }
      }
    }
    X<-abs(P-target)
    target.doses<-which(X == min(X))
    y<-colSums(rec.all)/sim
    correct<-sum(y[target.doses])
    
    output<-list(True.Probabilities=P,number.of.simulation=sim,experimentation=colSums(result.exp)/(n*sim),recommendations=y,correct.recommendation=correct,percentage.of.termination=counter/sim,mean.number.of.patients=mean(rowSums(result.exp)),mean.number.of.events=mean(result.suc))
    
    
    if (trajectories==FALSE){
      return(output)
    } else{
      return(output)
    }
    
    
  } else{
    if (assigment=="randomization"){
    } 
  }
}





