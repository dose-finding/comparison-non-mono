## Simulation study
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2020)

# the following code allow to reproduce results of the simulation study for the PO-CRM design
# as showed in tables 1 -4 of the paper.

######### PARAMETERS #########
#
# pI : vector of true toxiciy probabilities
# dose : vector of dose level number (ex : c(1,2,3,4))
# p.skel : matrix of skeleton values corresponding to the possible orderings
# prior.o : vector of selection probabilities associated to the possible ordering
#           (length = number of dose level)
# target : the target DLT rate
# cohortsize : size of each cohort inclusion
# ncohort : number of cohorts in a sinlge trial ##  cohortsize*ncohort = total number of patients in a trial
# start.comb : first dose level to start trial
# conf.level : confidence level for the confidence interval
# ntrial : number of simulated trials 
# mu : mean of normal prior distribution for A parameter
# scale : sd of normal prior distribution for A parameter
#
##############################

# source the functions
# setwd()
# source("pocrm_main_functions.R")

order <- matrix(nrow = 3, ncol = 6)
order[1,] = c(1,2,3,4,5,6)
order[2,] = c(1,2,4,3,5,6)
order[3,] = c(1,2,3,5,4,6)

p1 <- c(0.30, 0.40, 0.50, 0.60, 0.90, 0.90) # Scenario 1
p2 <- c(0.14, 0.30, 0.40, 0.50, 0.70, 0.75) # Scenario 2
p31 <- c(0.05, 0.10, 0.20, 0.30, 0.45, 0.70) # Scenario 3.1
p32 <- c(0.05, 0.10, 0.30, 0.20, 0.45, 0.70) # Scenario 3.2
p33 <- c(0.05, 0.10, 0.20, 0.45, 0.30, 0.70) # Scenario 3.3
p4 <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30) # Scenario 4
p5 <- c(0.001,0.01, 0.05, 0.07, 0.10, 0.15) # Scenario 5
p6 <- c(0.50, 0.70, 0.80, 0.90, 0.95, 0.95) # Scenario 6
p71 <- c(0.01, 0.05, 0.20, 0.50, 0.60, 0.70) # Scenario 7.1
p72 <- c(0.01, 0.05, 0.50, 0.20, 0.60, 0.70) # Scenario 7.2
p81 <- c(0.01, 0.05, 0.10, 0.25, 0.45, 0.60) # Scenario 8.1
p82 <- c(0.01, 0.05, 0.25, 0.10, 0.45, 0.60) # Scenario 8.2
p83 <- c(0.01, 0.05, 0.10, 0.45, 0.25, 0.60) # Scenario 8.3
p91 <- c(0.01, 0.05, 0.08, 0.10, 0.25, 0.45) # Scenario 9.1
p92 <- c(0.01, 0.05, 0.08, 0.25, 0.10, 0.45) # Scenario 9.2


# Results for the paper
set.seed(2306) ; pocrm.s1  <- bpocrm.sim(pI=p1 ,  dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s2  <- bpocrm.sim(pI=p2 ,  dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s31 <- bpocrm.sim(pI=p31 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s32 <- bpocrm.sim(pI=p32 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s33 <- bpocrm.sim(pI=p33 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s4  <- bpocrm.sim(pI=p4 ,  dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s5  <- bpocrm.sim(pI=p5 ,  dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s6  <- bpocrm.sim(pI=p6 ,  dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s71 <- bpocrm.sim(pI=p71 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s72 <- bpocrm.sim(pI=p72 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s81 <- bpocrm.sim(pI=p81 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s82 <- bpocrm.sim(pI=p82 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s83 <- bpocrm.sim(pI=p83 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s91 <- bpocrm.sim(pI=p91 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; pocrm.s92 <- bpocrm.sim(pI=p92 , dose=c(1:ncol(order)), p.skel=getwm(order, getprior(0.05,0.30,2,6)), prior.o = rep(1/nrow(order),nrow(order)), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)


pocrm.s33
# Printing the Results as in Table 1
kable(report(pocrm.s1,"PO-CRM"),row.names = F)
kable(report(pocrm.s1,"PO-CRM"),row.names = F, format="latex")
