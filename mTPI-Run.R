## Simulation study
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2020)

# the following code allow to reproduce results of the simulation study for the mTPI design
# as showed in tables 1 -4 of the paper.


######### PARAMETERS #########
#
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
#
##############################

# source the functions
# setwd()
# source("mtpi_main_functions.R")

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
set.seed(2306) ; mtpi.s1	  <- mtpi.sim(p=p1	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s2	  <- mtpi.sim(p=p2	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s31	  <- mtpi.sim(p=p31	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s32	  <- mtpi.sim(p=p32	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s33	  <- mtpi.sim(p=p33	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s4	  <- mtpi.sim(p=p4	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s5	  <- mtpi.sim(p=p5	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s6	  <- mtpi.sim(p=p6	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s71	  <- mtpi.sim(p=p71	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s72	  <- mtpi.sim(p=p72	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s81	  <- mtpi.sim(p=p81	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s82	  <- mtpi.sim(p=p82	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s83	  <- mtpi.sim(p=p83	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s91	  <- mtpi.sim(p=p91	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)
set.seed(2306) ; mtpi.s92	  <- mtpi.sim(p=p92	,D=6, pT=.3, startdose = 2, eps1=.10, eps2=.10, xi=.9,sampsize = 30, csize =3, simN=10000)

mtpi.s1

