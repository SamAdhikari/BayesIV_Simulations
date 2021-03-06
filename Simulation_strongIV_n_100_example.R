#If running individually uncomment following to load the required libraries

#library(BayesIV)
#library(ivpack)


##setting up MCMC parameters
reps = 50 #repitions
niter = 20000 #MCMC  iterations
burnin = 5000 #first 5000 samples will be discarded
thin = 20 ## thinning parameter

Chain = seq(burnin, niter, by = thin)

##SETTING UP PARAMETERS FOR DATA GENERATION
#############   #####           #########
n = 100 ##sample size

####Confounders
#heterogeneity is based on X_3
X = as.matrix( data.frame(rep(1, n),
                rnorm(n, 0, 1),
                rnorm(n, -1, 1.5),
                c( rep(1, n/2), rep(0, n/2))
        )     )

#Binary Instrument
Z = rbinom(n, 1, 0.5)

#Coefficients
beta1 = c(100, -0.5, 1.5, 10)
alpha1 = 0.1
beta0 = c(90, -0.5, 1.5, 0)
alpha0 = 0.1
alphaD = 0.2
gammaD = 1.5

pdfname = paste('Gamma_Strong_simulationRep_n', n, sep = '')

Theta = rnorm(n, 0, 0.1)

##SPECIFY THESE PARAMETERS FOR ERRORS WITH MIXTURE DISTRIBUTION INSTEAD OF GAMMA DISTRIBUTION
#components = sample(1:4, prob = c(0.15, 0.4, 0.25, 0.05),
#                    size = n, replace = TRUE)
#mus = c(-0.1, 1, 1, 10)
#sds = sqrt( c(0.01, 0.1, 0.1, 10))

##DATA GENERATION AND MODEL FITTING OVER 50 REPITITIONS
for(rr in 1:reps)
{
    ##specify errors on the potential outcomes
    err1 = rgamma(n, 3, 0.1)
    err0 = rgamma(n, 3, 0.1)
    
    ##generate potential outcomes
    Y1 = X%*%beta1 + alpha1*Theta  + err1
    Y0 = X%*%beta0 + alpha0*Theta + err0
    
    ##generate treatment assignment
    Dstar = -0.5 + X[,3] + Z*gammaD  + alphaD * Theta + rnorm(n,0,1)
    D = 1*(Dstar > 0)
    
    ##Observed outcome based on treatment assignment
    Yobs = D*Y1 + (1-D)*Y0
    
    ###MODEL FITTING

    ##1.  PROPOSED MODEL
    runModel_Normal_DPM = mcmcRun_Normal_DPM(Yobs = Yobs, Tr = D, X = X[,-1],
                                        Z = Z, niter = niter)
    
    save(runModel_Normal_DPM,
    file = paste('NormalDPM', pdfname, 'rep', rr, '.RData', sep = '')
    )
    
    ##2. BAYESIAN IV MODEL WITH NORMAL ERRORS
    runModel_Normal = mcmcRun_Normal(Yobs = Yobs, Tr = D, X = X, Z = Z,
                        niter = niter)
    
    save(runModel_Normal, file = paste('Normal', pdfname, 'rep', rr, '.RData', sep = '')
    )
    
    ##3. Two stage least squares mdodel
    IVfit0 = ivreg(Yobs ~ D + X[, 2:4]|Z + X[, 2:4])
    robust.se(IVfit0)
    
    ##also fit two stage separately for males and females for treatment effect heterogeneity
    femaleind = which(X[,4] == 0)
    XFemale = X[femaleind,]
    
    maleind = which(X[,4] == 1)
    XMale = X[maleind,]
    
    IVfit_Female = ivreg(Yobs[femaleind] ~ D[femaleind] +
                    X[femaleind, 2:3]| Z[femaleind] +
                    X[femaleind, 2:3] )
    
    IVfit_Male = ivreg(Yobs[maleind] ~ D[maleind] +
                X[maleind, 2:3]| Z[maleind] +
                        X[maleind, 2:3] )
    
    save(IVfit0, IVfit_Female, IVfit_Male,
                Y1,Y0, Yobs, D, X, Z,
                file = paste('Data', pdfname, 'rep', rr, '.RData', sep = '')
    )
    
    rm(runModel_Normal_DPM,runModel_Normal)
    
}
