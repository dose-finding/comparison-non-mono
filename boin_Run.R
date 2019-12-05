## Simulation study

# the following code allow to reproduce results of the simulation study for the BOIN design
# as showed in tables 1 -4 of the paper.

# source the functions
# install.packages("BOIN")
library(BOIN)


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

boin.design <- get.oc(target=0.3, p.true=p1, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
boin.design


# Results for the paper
set.seed(2306) ; boin.s1 <- get.oc(target=0.3, p.true=p1, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s2 <- get.oc(target=0.3, p.true=p2, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s31 <- get.oc(target=0.3, p.true=p31, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s32 <- get.oc(target=0.3, p.true=p32, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s33 <- get.oc(target=0.3, p.true=p33, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s4 <- get.oc(target=0.3, p.true=p4, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s5 <- get.oc(target=0.3, p.true=p5, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s6 <- get.oc(target=0.3, p.true=p6, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s71 <- get.oc(target=0.3, p.true=p71, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s72 <- get.oc(target=0.3, p.true=p72, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s81 <- get.oc(target=0.3, p.true=p81, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s82 <- get.oc(target=0.3, p.true=p82, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s83 <- get.oc(target=0.3, p.true=p83, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s91 <- get.oc(target=0.3, p.true=p91, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)
set.seed(2306) ; boin.s92 <- get.oc(target=0.3, p.true=p92, startdose=2,ncohort=10, cohortsize=3, extrasafe=F, ntrial=10000)

