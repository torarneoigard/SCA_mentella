# SCA_mentella.R in TMB
# Program designed to perform a statistical catch-at-age assessment based upon
# catch data (total catch and pelagic fleet catches) as reported to ICES and
# survey data (0-group, winter, ecosyste and Russian groundfish surveys),
#
# The present R file is fitting the model in TMB (together with SCA_mentella_model.cpp)
# the file SCA_mentella_data.R is preparing the data
# the file SCA_mentella_dataplots.R is plotting the input data
# the file SCA_mentella_plots.R is plotting the results
#
# The version in the directory "~/Documents/Work/Redfish/ICES/AFWG2016/mentella/SCA_mentella_TMB_2016/Version0" reproduce the SCA outputs from the AFWG2014
# Benjamin Planque, March-April 2015
# Benjamin Planque, Tor Arne ??ig??rd and Alf Harbitz, February 2016

load(paste(getwd(),"/SCA_mentella_data.Rdata",sep = ""))
YearSpan=data$minYear:data$maxYear

#data$pa0PropInit
#data$logb1PropInit
#data$logb2PropInit
#data$logb1PropMap

parameters <- list(
  logNY1=rep(17.1,data$nAges),    # log-numbers in year 1 (for all ages) 
  logNA1fe=rep(18.2,(data$nYears-1)), # log-numbers at age one (for all years except the first one) <- only to be used with the fixed effects model
  DemlogFY=rep(-4,data$nYears),   # partial log-Mortality by year for the Demersal Fleet (separable mortality)
  PellogFY=c(rep(-1e6,14),rep(-4,data$nYears-14)),  # partial log-Mortality by year for the Pelagic Fleet (separable mortality)
  pDema50=0,                      # Demersal fleet selectivity coefficient 1 (this should be bounded between 6 and maxAge)
  Demlogw=1,        		          # Demersal fleet selectivity coefficient 2
  pPela50=0,                      # Pelagic fleet selectivity coefficient 1 (this should be bounded between 6 and maxAge)
  Pellogw=1,        		          # Pelagic fleet selectivity coefficient 2
  DemlogVarLogC=0,                # (0) log of the variance of the logCatch, demersal fleet
  PellogVarLogC=0,       	        # (0) log of the variance of the logCatch, pelagic fleet
  logVarLogI=rep(-0.1,3),         # log of the variance of the logSurveyIndex
  logQSurvey = data$logQSurveyInit,
  # optional parameters (can be switched off)
  Demsplus=10,                    # Demersal fleet selectivity coefficient 3 (selectivity for +group), to be switched off
  Pelsplus=10,                    # Pelagic fleet selectivity coefficient 3 (selectivity for +group), to be switched off
  logM2=-3,					              # log of Natural mortality (should be switched off)
  pa0 = data$pa0Init,             # Survey selectivity coefficients a (probit-transformed age at which the derivative of the survey selectivity function is flat)
  pa0Prop = data$pa0PropInit,             # Survey selectivity coefficients a (probit-transformed age at which the derivative of the survey selectivity function is flat)
  logb1 = data$logb1Init,         # Survey selectivity coefficients b1 (log-transformed upward slope of the selectivity function, for age<pa)
  logb2 = data$logb2Init,         # Survey selectivity coefficients b2 (log-ransformed downward slope of the selectivity function, for age>pa)
  logb1Prop = data$logb1PropInit, # Survey selectivity coefficients b1 (log-transformed upward slope of the selectivity function, for age<pa)
  logb2Prop = data$logb2PropInit, # Survey selectivity coefficients b2 (log-ransformed downward slope of the selectivity function, for age>pa)
  palogNA1 = 1,                   # Probit-transformed coffecient a of the autoregressive recruitment model (random effects)
  logSigmalogNA1=0,               # log-transformed standard deviation of log-recruitment in the autoregressive recruitment model (random effects)
  ulogNA1=rep(0,(data$nYears-1))  # random effects proper 
)

compile("SCA_mentella_model.cpp")           # compile SCA model
dyn.load(dynlib("SCA_mentella_model"))      # load model

NoPellFy <- 1992:2005       # Fy for pelagic fleet is not estimated for the first 14y (1992-2005)
ind <- which(NoPellFy %in% YearSpan)
#Objective function depends on wether or not you want logNA1 as random or fixed effect
if(data$REswitch == 0){
  obj <- MakeADFun(data,parameters,DLL="SCA_mentella_model",checkParameterOrder = TRUE,map=list(
    PellogFY=factor(c(rep(NA,length(ind)),1:(length(YearSpan)-length(ind)))),       # Fy for pelagic fleet is not estimated for the first 14y (1992-2005)
    logQSurvey = data$logQSurveyMap,          # Survey scaling factor for the ecosystem survey is fixed
    Demsplus=factor(NA),                      # Demersal fleet selectivity for the +group is fixed
    Pelsplus=factor(NA),                      # Pelagic fleet selectivity for the +group is fixed
    logb1 = data$logb1Map,
    logb2 = data$logb2Map,
    logM2=factor(NA),                         # Natural mortality is fixed
    palogNA1 = factor(NA),                    # parameters for random effect are not estimated 
    logSigmalogNA1=factor(NA),                # parameters for random effect are not estimated 
    ulogNA1=factor(rep(NA,(data$nYears-1))))  # parameters for random effect are not estimated 
  )
} else {
  obj <- MakeADFun(data,parameters,random=c("ulogNA1"),DLL="SCA_mentella_model",checkParameterOrder = FALSE,map=list(
    PellogFY=factor(c(rep(NA,length(ind)),1:(length(YearSpan)-length(ind)))),       # Fy for pelagic fleet is not estimated for the first 14y (1992-2005)
    logQSurvey= data$logQSurveyMap,          # Survey scaling factor for the ecosystem survey is fixed
    #logQSurvey2=factor(NA),                   # Survey scaling factor for the ecosystem survey is fixed
    Demsplus=factor(NA),                      # Demersal fleet selectivity for the +group is fixed
    Pelsplus=factor(NA),                      # Pelagic fleet selectivity for the +group is fixed
    logb1 = data$logb1Map,
    logb2 = data$logb2Map,
    logM2=factor(NA),                         # Natural mortality is fixed
    logNA1fe=factor(rep(NA,(data$nYears-1))))   # Fixed effect is not estimated 
  )
}


# obj$fn()
# obj$gr()
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr,control = list(eval.max = 1e6,maxit = 1e6)))

report <- sdreport(obj)
R.code <- scan('SCA_mentella_model.R',what="",sep="\n")  # reads the current file and store it into the variable 'code'
cpp.code <- scan('SCA_mentella_model.cpp',what="",sep="\n")  # reads the current file and store it into the variable 'code'
date.flag=date()
model=list(date.flag=date.flag,data=data,parameters=parameters,R.code=R.code,cpp.code=cpp.code,obj=obj,opt=opt,rep=report)
save(model,file='SCA_mentella_model.Rdata')
