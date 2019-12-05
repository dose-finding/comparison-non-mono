# This code can be used to reproduce the results of the paper
#
#     A comparison of Phase I dose-escalation designs in clinical
#     trials with monotonicity assumption violation
#
#     by

#     Abbas, Rossoni, Jaki, Paoletti, and Mozgunov (2019)
#
#     for the NMA design.
#
#    Code NMA-Main-Code.R should be run first.

target<-0.30                                         # Target Toxicity Level
prior<-c(0.20, 0.23, 0.26, 0.29, 0.32, 0.35)         # Prior parameter $\nu$

# Specify the scenario
#
# true<-c(0.3,0.4,0.5,0.6,0.9,0.9)        # Scenario 1
# true<-c(0.14,0.3,0.4,0.5,0.7,0.75)      # Scenario 2
# true<-c(0.05,0.1,0.2,0.3,0.45,0.7)      # Scenario 3.1
# true<-c(0.05,0.1,0.3,0.2,0.45,0.7)      # Scenario 3.2
true<-c(0.05,0.1,0.2,0.45,0.3,0.7)      # Scenario 3.3
# true<-c(0.01,0.05,0.1,0.15,0.2,0.3)     # Scenario 4
# true<-c(0.00,0.01,0.05,0.07,0.1,0.15)   # Scenario 5
# true<-c(0.5,0.7,0.8,0.9,0.95,0.95)      # Scenario 6
# true<-c(0.01,0.05,0.20,0.50,0.60,0.70)  # Scenario 7.1
# true<-c(0.01,0.05,0.50,0.20,0.60,0.70)  # Scenario 7.2
# true<-c(0.01,0.05,0.10,0.25,0.45,0.60)  # Scenario 8.1
# true<-c(0.01,0.05,0.25,0.10,0.45,0.60)  # Scenario 8.2
# true<-c(0.01,0.05,0.10,0.45,0.25,0.60)  # Scenario 8.3
# true<-c(0.01,0.05,0.08,0.10,0.25,0.45)  # Scenario 9.1
# true<-c(0.01,0.05,0.08,0.25,0.10,0.45)  # Scenario 9.2
 

set.seed(123)
trial<-nma.design(P=true, # scenario
                  target=0.30, # target toxicity
                  n=30, # sample size
                  cohort=3, # cohort size
                  start.dose=2, # starting dose
                  Pr=prior, # prior parameters $\nu$
                  lambda=0.25, # parameter $\lambda$
                  nsims=10000, # number of simulations
                  final.prob=0.90,rate=0.005, # parameters of the safety constraint (2.7)
                  futility=T,threshold.low=0.25,final.prob.low=0.30) # parameters of the futility constraint (2.8)

# Creating matrix with the results as in Table  1
result.table<-mat.or.vec(3,10)
result.table[1,2:7]<-true
result.table[2,2:7]<-round(trial$recommendations*100,1)
result.table[2,8]<-round(trial$percentage.of.termination*100,1)
result.table[2,9]<-round(sum(trial$recommendations[which(true>target)])*100,1)
result.table[3,2:7]<-y2<-round(trial$experimentation*trial$mean.number.of.patients,1)
result.table[2,10]<-round(100*sum(y2[which(true>target)])/sum(y2),1)
result.table[1,8]<-"Stop"
result.table[1,9]<-"SelTox"
result.table[1,10]<-"%Tox"
result.table[1,1]<-"True"
result.table[2,1]<-"Selection"
result.table[3,1]<-"Patients"

# Printing the Results as in Table 1
# library("xtable")
xtable(result.table,digits=0)


