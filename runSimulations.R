#script to generate simulated data and run IV models, for 50 replications
#followed by code to summarize the posterior draws of the fitted model


##If Using R Package Installed From the .tar.zip File
##
##for MAC or linux OS
#####
library(BayesIV)
library(ivpack)


##For windows users uncomment following
#####
#####
# library(Rcpp)
# library(MASS)
# library(parallel)
# library(ivpack)
# library(msm) ##for sampling from truncated normal distribution
# library(DPpackage) ##for sampling from dirichlet process mixture prior
# 
# source('MCMC_IV_DirichletLatentFactor.R')
# source('MCMC_NormalLatentFactor.R')
# source('MCMC_IV_Normal_DPM.R')
# source('getATE_posterior.R')
# source('getATE_posterior_NDPM.R)
# source('getTT_posterior.R')
# source('getTT_posterior_NDPM.R)
# source('Gibbs_Reg.R)

##simulation code along with code to summarize the results from simulation
##n = 100
##replications = 50

source('Simulation_strongIV_n_100_example.R')
source('SummarizeSimulations_n_100.R')

##simulation code along with code to summarize the results from simulation
##n = 500
##replications = 50

source('Simulation_strongIV_n_500_example.R')
source('SummarizeSimulations_n_500.R')

##simulation code along with code to summarize the results from simulation
##n = 2000
##replications = 50

source('Simulation_strongIV_n_2000_example.R')
source('SummarizeSimulations_n_2000.R')


