# SCA_mentella_projections.R in TMB
# Program designed to project the stock using catch-at-age assessment based upon
# catch data (total catch and pelagic fleet catches) as reported to ICES and
# survey data (0-group, winter, ecosyste and Russian groundfish surveys),
#
# Benjamin Planque, April 2017, revised 2018

# WARNING: THIS IS NOT WORKING PROPERLY
# THE CODE IS REVISED TO DEAL WITH VARYING SELECTIVITY
# OF THE DEMERSAL FLEET BUT THE CALCULATED CATCHES ARE NOT IDENTICAL
# TO THOSE FROM SCA OUTPUT

require(tidyverse)

Fmultiplier=1 #2.051282 #1.320588 # Multiplier for the fishing mortality for both fleets (1.67 is the value to reach F0.1=0.039)

load('SCA_mentella_model.Rdata')            # Load model results from TMB
data=model$data                             # extract data
YearSpan=data$minYear:data$maxYear
AgeSpan=data$minAge:data$maxAge
parameters=model$parameters                 # and parameters

rep.matrix <- summary(model$rep)            # get the list of reported values and standard deviations
rep.rnames <- rownames(rep.matrix)          # get the names of variables
indlogNY1 <- which(rep.rnames=="logNY1")    # extract line numbers for numbers in year one
indlogNA1 <- which(rep.rnames=="logNA1")    # extract line numbers for numbers at age one
# if(data$RElogNA == 0){
# indlogNA1 <- which(rep.rnames=="logNA1")    # extract line numbers for numbers at age one
# } else {
#   indlogNA1 <- which(rep.rnames == "logNA1re")
# }
if (REDemFishMort==0){
  indDemlogFY <- which(rep.rnames=="DemlogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
} else{
  indDemlogFY <- which(rep.rnames=="DemlogFYRE") # extract line numbers for demersal fishing mortality in years (Fy's)
}
if (REPelFishMort==0){
  indPellogFY <- which(rep.rnames=="PellogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
}else{
  indPellogFY <- which(rep.rnames=="PellogFYRE") # extract line numbers for demersal fishing mortality in years (Fy's)
}
if (REDemFishSel==0){
  indlogitDemFAfe <- which(rep.rnames=="logitDemFAfe") # extract line numbers for demersal fishing mortality in years (Fy's)
}else{
  indlogitDemFARE <- which(rep.rnames=="logitDemFARE") # extract line numbers for demersal fishing mortality in years (Fy's)
}
indlogitPelFA <- which(rep.rnames=="logitPelFA") # extract line numbers for demersal fishing mortality in years (Fy's)
indSA <- which(rep.rnames=="SA")
indlogSA <- which(rep.rnames=="logSA")
indlogSSB <- which(rep.rnames=="logSSB")    # extract line numbers for SSB
indlogTriN <- which(rep.rnames=="logTriN")  # extract line numbers for triangular population matrix

# output vectors for plotting
logNY1 <- rep.matrix[indlogNY1,1]
logNY1.sd <- rep.matrix[indlogNY1,2]
logNA1 <- c(logNY1[1],rep.matrix[indlogNA1,1])
logNA1.sd <- c(logNY1.sd[1],rep.matrix[indlogNA1,2])
DemlogFY <- rep.matrix[indDemlogFY,1]
DemlogFY.sd <- rep.matrix[indDemlogFY,2]
PellogFY <- c(rep(-1e6,14),rep.matrix[indPellogFY,1])
PellogFY.sd <- c(rep(0,14),rep.matrix[indPellogFY,2])
if (REDemFishSel==0){
  logitDemFA <-  rep.matrix[indlogitDemFA,1]
  logitDemFA.sd <-  rep.matrix[indlogitDemFARE,2]
} else {
  logitDemFARE <-  rep.matrix[indlogitDemFARE,1]
  logitDemFARE.sd <-  rep.matrix[indlogitDemFARE,2]
}
logitPelFA <-  rep.matrix[indlogitPelFA,1]
logitPelFA.sd <-  rep.matrix[indlogitPelFA,2]
SA <- matrix(rep.matrix[indSA,1],nrow = length(surveys),byrow = FALSE)
logSA <- matrix(rep.matrix[indlogSA,1],nrow = length(surveys),byrow = FALSE)
logSA.sd <- matrix(rep.matrix[indlogSA,2],nrow = length(surveys),byrow = FALSE)

logSSB <- rep.matrix[indlogSSB,1]
logSSB.sd <- rep.matrix[indlogSSB,2]

logQSurvey <- array(NA,length(surveys))
logQSurveytemp = rep.matrix[rep.rnames=="logQSurvey",1]
oneless <- 0
for (i in 1:length(surveys)){
  if(!is.na(data$logQSurveyMap[i])){
    logQSurvey[i] <- logQSurveytemp[i-oneless]
  }  else{
    logQSurvey[i] <- parameters$logQSurvey[i]
    oneless <- oneless + 1
  } 
}

# reformat the triangular population matrix
logTriN <- rep.matrix[indlogTriN,1]
logTriN.std <- rep.matrix[indlogTriN,2]
logTriNmatrix <-  (matrix(logTriN,nrow = (data$nYears+1),byrow = FALSE))
logTriNmatrix.std <- (matrix(logTriN.std,nrow = (data$nYears+1),byrow = FALSE))

# Fishing mortalities
FY=data.frame(Year=data$minYear:data$maxYear,
              DemFY=exp(DemlogFY),
              PelFY=exp(PellogFY))
new.FY=tail(FY,1)*Fmultiplier
FY=rbind(FY,tail(FY,1)) # keep the last year F
for (i in 2:5) FY=rbind(FY,new.FY)

FARE=data.frame(year=rep(YearSpan,18),
                age=rep(2:19,each=data$nYears),
                DemFARE=exp(logitDemFARE)/(1+exp(logitDemFARE)),
                PelFA=rep(exp(logitPelFA)/(1+exp(logitPelFA)),each=data$nYears))

# Filling up the fishing mortality matrix
DemFA=as.matrix(dplyr::select(FARE,year,age,DemFARE) %>% tidyr::spread(age,DemFARE) %>% dplyr::select(-year)) # Demersal fishing mortality prepared as a matrix
for(i in 1:5){DemFA=rbind(DemFA,DemFA[26,])} # extra 5 years added with identical selectivity pattern
PelFA=as.matrix(dplyr::select(FARE,year,age,PelFA) %>% tidyr::spread(age,PelFA) %>% dplyr::select(-year)) # Pelagic fishing mortality prepared as a matrix
for(i in 1:5){PelFA=rbind(PelFA,DemFA[26,])} # extra 5 years added with identical selectivity pattern
DemFY=FY$DemFY%*%t(rep(1,length(logNY1))) # DemFy as a matrix
PelFY=FY$PelFY%*%t(rep(1,length(logNY1))) # PelFy as a matrix
F.Matrix=DemFA*DemFY+PelFA*PelFY # Total mortality Matrix

# Projection of the population matrix 5y into the future ------------------
Pop.Matrix=matrix(data=NA,nrow=length(YearSpan)+5,ncol=length(logNY1)) # initialise population matrix
Pop.Matrix[,1]=c(exp(logNA1),rep(exp(mean((logNA1))),5)) # recruits (last 5y are Long Term Geometric Mean)
Pop.Matrix[1,]=exp(logNY1) # number of different ages in year one

# Calculating numbers at age
last.age=dim(Pop.Matrix)[2]
for (year in 2:(dim(Pop.Matrix)[1])){ # loop on years
  for (age in 2:(dim(Pop.Matrix)[2]-1)){ # loop on ages, leaving out the +group
    Z=F.Matrix[year-1,age-1]+exp(parameters$logM2)
    Pop.Matrix[year,age]=Pop.Matrix[year-1,age-1]*exp(-Z)
  }
  Z=F.Matrix[year-1,last.age-1]+exp(parameters$logM2)
  Zplus=F.Matrix[year-1,last.age]+exp(parameters$logM2) # mortality in the plus group
  Pop.Matrix[year,age+1]=Pop.Matrix[year-1,last.age-1]*exp(-Z)+Pop.Matrix[year-1,last.age]*exp(-Zplus) # numbers in the plus group
}

# Calculating SSB & TSB
MaturityAtAge=data$MaturityAtAge
for (i in 1:5) MaturityAtAge=rbind(MaturityAtAge,tail(MaturityAtAge,1))
WeightAtAge=data$WeightAtAge
for (i in 1:5) WeightAtAge=rbind(WeightAtAge,tail(WeightAtAge,1))
TSB.Matrix=Pop.Matrix*WeightAtAge
TSB=rowSums(TSB.Matrix)/1000 # SSB (in tonnes)
SSB.Matrix=Pop.Matrix*MaturityAtAge*WeightAtAge
SSB=rowSums(SSB.Matrix)/1000 # SSB (in tonnes)

# Calculating F12-18 and F19
F.12.18=rowMeans(F.Matrix[,11:17])
F.19=F.Matrix[,18]

# Calculating Catches
Catch.Matrix=(F.Matrix/(F.Matrix+exp(parameters$logM2)))*(1-exp(-(F.Matrix+exp(parameters$logM2))))*Pop.Matrix
Catch.in.Tonnes=Catch.Matrix*WeightAtAge/1000
Total.Catch=rowSums(Catch.in.Tonnes)

# prepare output table for export
Table2export=data.frame(Year=c(data$minYear:(data$maxYear+5)),
                        Recruits.at.age.2=round(Pop.Matrix[,1]),
                        Recruits.at.age.6=round(Pop.Matrix[,5]),
                        Total.Stock.Biomass=round(TSB),
                        Spawning.Stock.Biomass=round(SSB),
                        F12.18=round(F.12.18, digits=3),
                        F19=round(F.19,digits=3),
                        Total.Catch.in.Tonnes=round(Total.Catch),
                        Fmultiplier=c(rep(NA,data$nYears),rep(Fmultiplier,5)))

write.table(Table2export,'SCA_mentella_projections.txt',row.names = FALSE,quote = FALSE,sep='\t')

# quick plots
par(mfrow=c(1,2))
plot(SSB) # the new calculations
lines(exp(logSSB)) # the output from SCAA

plot(Total.Catch)
lines(data$TotalCatches[,2]) # input data


# export table D9 (pop matrix)
nyears=length(YearSpan)
nyears=length(YearSpan)
Table9=data.frame(format(Pop.Matrix[1:nyears,]/1000,digits = 4))
colnames(Table9)=2:19
Table9=cbind(round(FY[1:nyears,c(2,3,1)]*1000)/1000,Table9)

# STOP HERE
# FA2=t(rbind(data.frame(DemFA=rep(NA,3),PelFA=rep(NA,3)),round(FA[,2:3]*1000)/1000))
# colnames(FA2)=colnames(Table9)
# Table9=rbind(FA2,Table9)
# Table9
# write.table(Table9,'Table9_pop.matrix.txt',row.names = FALSE,quote = FALSE,sep='\t')
