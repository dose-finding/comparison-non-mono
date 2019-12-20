## Simulation study
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2020)

# the following code allow to reproduce results of the simulation study for the CRM design
# as showed in tables 1 -4 of the paper.

# source the functions
# setwd()
source("crm_main_functions.R")
 

set.seed(2306)

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
set.seed(2306) ; crm.s1  <- bcrm.sim(pI=p1  , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s2  <- bcrm.sim(pI=p2  , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s31 <- bcrm.sim(pI=p31 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s32 <- bcrm.sim(pI=p32 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s33 <- bcrm.sim(pI=p33 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s4  <- bcrm.sim(pI=p4  , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s5  <- bcrm.sim(pI=p5  , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s6  <- bcrm.sim(pI=p6  , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s71 <- bcrm.sim(pI=p71 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s72 <- bcrm.sim(pI=p72 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s81 <- bcrm.sim(pI=p81 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s82 <- bcrm.sim(pI=p82 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s83 <- bcrm.sim(pI=p83 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s91 <- bcrm.sim(pI=p91 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)
set.seed(2306) ; crm.s92 <- bcrm.sim(pI=p92 , p.skel=getprior(0.05,0.30,2,6), target=0.30, cohortsize=3,ncohort=10, start.comb=2, conf.level=0.80, mu=0, scale=0.75, ntrial=10000)

crm.s1
