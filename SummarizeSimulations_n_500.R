library(ivpack)
source('helperfunctions.R')

ate_rslt_all = list()
CATE_rslt_all = list()


n = 500
type = 'Strong'
reps =  50
niter = 20000
burnin = 5000
thin = 10
Chain = seq(from = burnin, to = niter, by = thin)

runModels_DP_list = list()
runModel_Normal_DPM_list = list()
runModel_Normal_list = list()
dataList = list()
IVfit = list()
IVfit_Female_list = list()
IVfit_Male_list = list()

for(rr in 1:reps){
    load( paste('DataGamma_', type, '_simulationRep_n', n, 'rep', rr, '.RData', sep = '') )
    load( paste('NormalDPMGamma_', type, '_simulationRep_n', n, 'rep', rr, '.RData', sep = '') )
    load( paste('NormalGamma_', type, '_simulationRep_n', n, 'rep', rr, '.RData', sep = '') )
    
    dataList[[rr]] = list(Yobs = Yobs, D = D, Y1 = Y1, Y0 = Y0)
    
    runModel_Normal_DPM_list[[rr]] = runModel_Normal_DPM
    runModel_Normal_list[[rr]] = runModel_Normal
    
    IVfit[[rr]] = IVfit0
    IVfit_Female_list[[rr]] = IVfit_Female
    IVfit_Male_list[[rr]] = IVfit_Male
}

runModel_Normal_DPM = lapply(1:reps, function(x) runModel_Normal_DPM_list[[x]] )
runModels_Normal = lapply(1:reps, function(x) runModel_Normal_list[[x]] )

IVfit = lapply(1:reps, function(x) IVfit[[x]] )
IVfit_Female = lapply(1:reps, function(x) IVfit_Female_list[[x]] )
IVfit_Male = lapply(1:reps, function(x) IVfit_Male_list[[x]] )

dataList = lapply(1:reps, function(x) dataList[[x]] )

#####   AVERAGE TREATMENT EFFECTS   #######
#######################################

ATE_post_NDPM = lapply(1:reps, function(ww){
    ATEpost_NDPM = getATE_posterior_NDPM(fittedModel = runModel_Normal_DPM[[ww]],
        X = X[,-1],niter = niter,
        burnin = burnin, thin = thin)
    ATE_Chain_N_DPM = ATEpost_NDPM$ATE
    ATEhat_NDPM = mean(ATE_Chain_N_DPM)
    ATEci_NDPM = quantile(ATE_Chain_N_DPM, c(0.025, 0.5, 0.975) )
    
    return( list( ATEhat = ATEhat_NDPM, ATEci = ATEci_NDPM))
})

####
ATE_post_Normal = lapply(1:reps, function(x) {
    ATE_post = getATE_posterior(fittedModel = runModels_Normal[[x]],
        X = X, niter = niter, burnin = burnin, thin = thin)
    
    ATE_Chain_Normal = ATE_post$ATE
    ##Average treatment effect
    ATEhat_Normal = mean(ATE_Chain_Normal)
    ATEci_Normal = quantile(ATE_Chain_Normal, c(0.025, 0.5, 0.975) )
    
    return(list(ATEhat = ATEhat_Normal, ATEci = ATEci_Normal))
})

ATEci_2SLS = lapply(1:reps , function(x) {
    c(robust.se( IVfit[[x]])[2] - 1.96 * robust.se(IVfit[[x]])[2,2],
        robust.se( IVfit[[x]])[2],
        robust.se( IVfit[[x]])[2] + 1.96 * robust.se(IVfit[[x]])[2,2])
})


##compute bias on ATE
TrueATE = sapply(1:reps, function(x) mean( dataList[[x]]$Y1 - dataList[[x]]$Y0) )

#bias
biasNDPM = sapply( 1:reps , function(x) TrueATE[x] - ATE_post_NDPM[[x]]$ATEhat)
biasNormal = sapply(1:reps , function(x) TrueATE[x] -  ATE_post_Normal[[x]]$ATEhat)
bias2SLS =  sapply( 1:reps , function(x) TrueATE[x] - robust.se(IVfit[[x]])[2])

#width of credible interval
CIwidthNDPM = sapply(1:reps ,
                    function(x) ATE_post_NDPM[[x]]$ATEci[3] - ATE_post_NDPM[[x]$ATEci[1] )
CIwidthNormal = sapply(1:reps ,
                    function(x) ATE_post_Normal[[x]]$ATEci[3] - ATE_post_Normal[[x]]$ATEci[1] )
CIwidth2SLS = sapply(1:reps ,
                    function(x) ATEci_2SLS[[x]][3] - ATEci_2SLS[[x]][1] )


#\% coverage of 95% credible interval
coverageNDPM = sum( sapply(1:reps, function(x)
                1 * (TrueATE[x] >= ATE_post_NDPM[[x]]$ATEci[1]  &
                    TrueATE[x] <= ATE_post_NDPM[[x]]$ATEci[3])) ) / reps

coverageNormal = sum( sapply(1:reps, function(x)
                1 * (TrueATE[x] >= ATE_post_Normal[[x]]$ATEci[1]  &
                    TrueATE[x] <= ATE_post_Normal[[x]]$ATEci[3])) ) / reps

coverage2SLS = sum( sapply(1:reps, function(x)
                    1 * (TrueATE[x] >= ATEci_2SLS[[x]][1] &
                    TrueATE[x] <= ATEci_2SLS[[x]][3])) ) / reps

##SUMMARIZE THE RESULTS
#mean overall bias
mean(abs(biasNDPM))
mean(abs(biasNormal))
mean(abs(bias2SLS))

#mean overall CI width
mean(CIwidthNDPM)
mean(CIwidthNormal)
mean(CIwidth2SLS)

#% Coverage of CI
coverageNDPM
coverageNormal
coverage2SLS


ate_rslt = data.frame( array( NA, dim = c(3, 4)))

ate_rslt[, 1] = c('NPDM', 'Normal', '2SLS')

ate_rslt[, 2] = c(mean( abs(biasNDPM)),
                mean( abs(biasNormal)),
                mean( abs(bias2SLS)) )

ate_rslt[, 3] = c( mean(CIwidthNDPM),
                    mean(CIwidthNormal),
                    mean(CIwidth2SLS) )

ate_rslt[, 4] = c(coverageNDPM,
                coverageNormal,
                coverage2SLS )

###HETEROGENOUS TREATMENT EFFECT##########
###########################################
#### ATE X = 1 and X = 0
femaleind = which(X[, 4] == 0)
XFemale = X[femaleind, ]

maleind = which(X[, 4] == 1)
XMale = X[maleind, ]

ATE_female_true =sapply(1:reps,function(x)
            mean( dataList[[x]]$Y1[femaleind] -
                    dataList[[x]]$Y0[femaleind]) )

ATE_male_true = sapply(1:reps,function(x)
                mean(   dataList[[x]]$Y1[maleind]-
                    dataList[[x]]$Y0[maleind])  )

obs_Fem = mean( Yobs[which(D == 1 & X[, 3] == 0)]) -
        mean( Yobs[which(D == 0 & X[, 3] == 0)])

obs_Male = mean( Yobs[which(D == 1 & X[, 3] == 1)]) -
        mean( Yobs[which(D == 0 & X[, 3] == 1)])

HATEci_2SLS_F = lapply(1:reps,function(x){
    c( robust.se(IVfit_Female[[x]])[2] - 1.96 * robust.se(IVfit_Female[[x]])[2, 2],
        robust.se(IVfit_Female[[x]])[2],
        robust.se(IVfit_Female[[x]])[2] + 1.96 * robust.se(IVfit_Female[[x]])[2, 2])
}   )

HATEci_2SLS_M = lapply(1:reps, function(x){
    c(  robust.se(IVfit_Male[[x]])[2] - 1.96 * robust.se(IVfit_Male[[x]])[2, 2],
        robust.se(IVfit_Male[[x]])[2] ,
    robust.se(IVfit_Male[[x]])[2] + 1.96 * robust.se(IVfit_Male[[x]])[2, 2])
}    )

HATE_Normal = lapply(1:reps,function(x){
    ATE_post_Normal = getATE_posterior( fittedModel = runModel_Normal,
                            X = X, niter = niter, burnin = burnin, thin = thin)
    
    Y1hat_Chain_Normal = ATE_post_Normal$Y1hat
    Y0hat_Chain_Normal = ATE_post_Normal$Y0hat
    
    Y1hatmean_Normal = apply(Y1hat_Chain_Normal, 1 ,mean)
    Y0hatmean_Normal = apply(Y0hat_Chain_Normal, 1, mean)
    
    ##ATE X
    ATE_X_Chain_Normal = sapply(1:length(Chain),function(x)
            mean( Y1hat_Chain_Normal[which(X[,4] == 1), x] -
                Y0hat_Chain_Normal[which(X[,4] == 1), x]) )
                
    ATE_X_Normal = mean(ATE_X_Chain_Normal)
    ATE_X_Normal_ci = quantile(ATE_X_Chain_Normal, c(0.025, 0.5, 0.975))
    
    ATE_F_Chain_Normal  = sapply(1:length(Chain), function(x)
            mean(Y1hat_Chain_Normal[which(X[,4] == 0), x] -
                Y0hat_Chain_Normal[which(X[,4] == 0), x]) )
                
    ATE_F_Normal = mean(ATE_F_Chain_Normal)
    ATE_F_Normal_ci = quantile(ATE_F_Chain_Normal, c(0.025, 0.5, 0.975))
    
    return( list( ATE_X_Normal = ATE_X_Normal, ATE_X_Normal_ci = ATE_X_Normal_ci,
                ATE_F_Normal = ATE_F_Normal, ATE_F_Normal_ci = ATE_F_Normal_ci) )
})


####ATE X
HATE_NDPM = lapply(1:reps,function(x)
{
    mean1 = apply(runModel_Normal_DPM[[x]]$muY1, 2 , mean)
    mean0 = apply(runModel_Normal_DPM[[x]]$muY0, 2, mean)
    
    M1 = array(NA, dim = c(n, niter))
    M0 = array(NA, dim = c(n, niter))
    
    M1[which(dataList[[x]]$D == 1), ] = runModel_Normal_DPM[[x]]$muY1
    for(z in 1:niter){
        M1[which(dataList[[x]]$D == 0), z] = mean1[z]
    }
    
    M0[which(dataList[[x]]$D == 0),] = runModel_Normal_DPM[[x]]$muY0
    for(y in 1:niter){
        M0[which(dataList[[x]]$D == 1), y] = mean0[y]
    }
    
    fittedModel = runModel_Normal_DPM[[x]]
    
    Y1hat_Chain_NDPM = sapply(Chain, function(zz){
        X[,-1] %*% fittedModel$Beta1[, zz] +
        fittedModel$Theta[, zz] * fittedModel$Alpha1[zz] +
        M1[,zz]    })
    
    Y0hat_Chain_NDPM = sapply(Chain, function(zz){
        X[,-1] %*% fittedModel$Beta0[, zz] +
        fittedModel$Theta[, zz] * fittedModel$Alpha0[zz] +
        M0[, zz]   })
    
    ATE_X_Chain_NDPM = sapply(1:length(Chain), function(zz)
                mean(Y1hat_Chain_NDPM[which(X[, 4] == 1), zz]) -
                mean(Y0hat_Chain_NDPM[which(X[, 4] == 1), zz]))
    
    ATE_X_NDPM = mean(ATE_X_Chain_NDPM)
    ATE_X_NDPM_ci = quantile(ATE_X_Chain_NDPM, c(0.025, 0.5, 0.975))
    
    ATE_X_Chain_NDPM = sapply(1:length(Chain),function(zz)
                mean(Y1hat_Chain_NDPM[which(X[, 4] == 1), zz]) -
                    mean(Y0hat_Chain_NDPM[which(X[, 4] == 1), zz]))
    
    ATE_X_NDPM = mean(ATE_X_Chain_NDPM)
    ATE_X_NDPM_ci = quantile(ATE_X_Chain_NDPM, c(0.025, 0.5, 0.975))
    
    ATE_F_Chain_NDPM = sapply(1:length(Chain), function(zz)
                mean(Y1hat_Chain_NDPM[which(X[,4] == 0), zz]) -
                    mean(Y0hat_Chain_NDPM[which(X[,4] == 0), zz]))
    
    ATE_F_NDPM = mean(ATE_F_Chain_NDPM)
    ATE_F_NDPM_ci = quantile(ATE_F_Chain_NDPM, c(0.025, 0.5, 0.975))
    
    return(list(ATE_X_NDPM = ATE_X_NDPM, ATE_X_NDPM_ci = ATE_X_NDPM_ci,
                        ATE_F_NDPM = ATE_F_NDPM, ATE_F_NDPM_ci = ATE_F_NDPM_ci))
})

##BIAS HATE female
HATE_F_biasNDPM = sapply(1:reps, function(x) ATE_female_true[x] - HATE_NDPM[[x]]$ATE_F_NDPM )
HATE_F_biasNormal = sapply(1:reps, function(x) ATE_female_true[x] -  HATE_Normal[[x]]$ATE_F_Normal )
HATE_F_bias2SLS =  sapply(1:reps, function(x) ATE_female_true[x] - HATEci_2SLS_F[[x]][2] )

#CI width
CIwidthNDPM_F = sapply(1:reps, function(x) HATE_NDPM[[x]]$ATE_F_NDPM_ci[3] - HATE_NDPM[[x]]$ATE_F_NDPM_ci[1] )
CIwidthNormal_F = sapply(1:reps, function(x) HATE_Normal[[x]]$ATE_F_Normal_ci[3] - HATE_Normal[[x]]$ATE_F_Normal_ci[1] )
CIwidth2SLS_F = sapply(1:reps, function(x) HATEci_2SLS_F[[x]][3] - HATEci_2SLS_F[[x]][1] )

#% Coverage
coverageNDPM_F = sum( sapply(1:reps, function(x)
                        1 * (ATE_female_true[x] >= HATE_NDPM[[x]]$ATE_F_NDPM_ci[1] &
                            ATE_female_true[x] <= HATE_NDPM[[x]]$ATE_F_NDPM_ci[3])) ) / reps

coverageNormal_F = sum( sapply(1:reps, function(x)
                    1 * (ATE_female_true[x] >= HATE_Normal[[x]]$ATE_F_Normal_ci[1] &
                        ATE_female_true[x] <= HATE_Normal[[x]]$ATE_F_Normal_ci[3])) ) / reps

coverage2SLS_F = sum( sapply(1:reps, function (x)
                    1 * (ATE_female_true[x] >= HATEci_2SLS_F[[x]][1] &
                        ATE_female_true[x] <= HATEci_2SLS_F[[x]][3])) )/ reps

######Bias, CI width, Coverage MALE
##BIAS HATE female
HATE_X_biasNDPM = sapply(1:reps,function(x)ATE_male_true[x] - HATE_NDPM[[x]]$ATE_X_NDPM)
HATE_X_biasNormal = sapply(1:reps,function(x)ATE_male_true[x] -  HATE_Normal[[x]]$ATE_X_Normal)
HATE_X_bias2SLS =  sapply(1:reps,function(x)ATE_male_true[x] - robust.se(IVfit_Male[[x]])[2])

#width of CI
CIwidthNDPM_X = sapply(1:reps, function(x) HATE_NDPM[[x]]$ATE_X_NDPM_ci[3] - HATE_NDPM[[x]]$ATE_X_NDPM_ci[1])
CIwidthNormal_X = sapply(1:reps, function(x) HATE_Normal[[x]]$ATE_X_Normal_ci[3] - HATE_Normal[[x]]$ATE_X_Normal_ci[1])
CIwidth2SLS_X = sapply(1:reps, function(x) HATEci_2SLS_M[[x]][3] - HATEci_2SLS_M[[x]][1])

#Coverage
coverageNDPM_X = sum( sapply(1:reps, function(x)
                    1 * (ATE_male_true[x] >= HATE_NDPM[[x]]$ATE_X_NDPM_ci[1] &
                        ATE_male_true[x] <= HATE_NDPM[[x]]$ATE_X_NDPM_ci[3])) ) / reps

coverageNormal_X = sum( sapply(1:reps, function(x)
                    1 * (ATE_male_true[x] >= HATE_Normal[[x]]$ATE_X_Normal_ci[1] &
                        ATE_male_true[x] <= HATE_Normal[[x]]$ATE_X_Normal_ci[3])) ) / reps

coverage2SLS_X = sum(sapply(1:reps, function(x)
                    1 * (ATE_male_true[x] >= HATEci_2SLS_M[[x]][1] &
                        ATE_male_true[x] <= HATEci_2SLS_M[[x]][3])) ) / reps

#mean bias for females
mean( abs(HATE_F_biasNormal) )
mean( abs(HATE_F_biasNDPM) )
mean( abs(HATE_F_bias2SLS) )

#mean bias for males
mean( abs(HATE_X_biasNormal) )
mean( abs(HATE_X_biasNDPM) )

coverageNDPM_F
coverageNDPM_X
coverageNormal_F
coverageNormal_X
coverage2SLS_F


CATE_rslt = data.frame(array(NA,dim = c(6, 4)) )
dimnames(CATE_rslt)[[2]] = c('Method', 'Mean Abs Bias', 'Mean CI Width','% Coverage')
CATE_rslt[, 1] = c('DPMLIV_F',
                'NormalLIV_F',
                '2SLS_F',
                'DPMLIV_M',
                'NormalLIV_M',
                '2SLS_M')

CATE_rslt[, 2] = c(mean(abs(HATE_F_biasNDPM)),
                    mean(abs(HATE_F_biasNormal)),
                    mean(abs(HATE_F_bias2SLS)),
                    mean(abs(HATE_X_biasNDPM)),
                    mean(abs(HATE_X_biasNormal)),
                    mean(abs(HATE_X_bias2SLS)))

CATE_rslt[, 3] = c(mean(CIwidthNDPM_F),
                    mean(CIwidthNormal_F),
                    mean(CIwidth2SLS_F),
                    mean(CIwidthNDPM_X),
                    mean(CIwidthNormal_X),
                    mean(CIwidth2SLS_X))

CATE_rslt[, 4] = c(coverageNDPM_F,
                coverageNormal_F,
                coverage2SLS_F,
                coverageNDPM_X ,
                coverageNormal_X,
                coverage2SLS_X)


rm(runModel_Normal,
runModel_Normal_DPM,
runModel_Normal_DPM_list,
runModel_Normal_list,
runModels_DP_list,
runModels_Normal)



save(ate_rslt, CATE_rslt,
            file=paste('RsltGamma', type, '_n', n, '.RData',sep = ''))


