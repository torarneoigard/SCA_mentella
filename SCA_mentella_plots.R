# SCA_mentella_plot.R in TMB
# Program designed to plot the results from catch-at-age assessment based upon
# catch data (total catch and pelagic fleet catches) as reported to ICES and
# survey data (0-group, winter, ecosyste and Russian groundfish surveys),
#
# The present R file is fitting the model in TMB (together with SCA_mentella_model.cpp)
# the file SCA_mentella_data.R is preparing the data
# the file SCA_mentella_plots.R is plotting the results
#
# The version in the directory "~/Documents/Work/Redfish/ICES/AFWG2015/Mentella/SCAA_mentella_tmb" reproduce the SCA outputs from the AFWG2014
# Benjamin Planque, February 2016

load('SCA_mentella_model.Rdata')            # Load model results from TMB
data=model$data                             # extract data
parameters=model$parameters                 # and parameters

rep.matrix <- summary(model$rep)            # get the list of reported values and standard deviations
rep.rnames <- rownames(rep.matrix)          # get the names of variables
indlogNY1 <- which(rep.rnames=="logNY1")    # extract line numbers for numbers in year one
if(REswitch == 0){
indlogNA1 <- which(rep.rnames=="logNA1")    # extract line numbers for numbers at age one
} else {
  indlogNA1 <- which(rep.rnames == "logNA1re")
}
indDemlogFY <- which(rep.rnames=="DemlogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
indPellogFY <- which(rep.rnames=="PellogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
indlogitDemFA <- which(rep.rnames=="logitDemFA") # extract line numbers for demersal fishing mortality in years (Fy's)
indlogitPelFA <- which(rep.rnames=="logitPelFA") # extract line numbers for demersal fishing mortality in years (Fy's)
indSAWinter <- which(rep.rnames=="SAWinter") # extract line numbers for winter survey selectivity at age
indSAEco <- which(rep.rnames=="SAEco") # extract line numbers for ecosystem survey selectivity at age
indSARussian <- which(rep.rnames=="SARussian") # extract line numbers for Russian groundfish survey selectivity at age
indlogSAWinter <- which(rep.rnames=="logSAWinter") # extract line numbers for winter survey selectivity at age
indlogSAEco <- which(rep.rnames=="logSAEco") # extract line numbers for ecosystem survey selectivity at age
indlogSARussian <- which(rep.rnames=="logSARussian") # extract line numbers for Russian groundfish survey selectivity at age
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
logitDemFA <-  rep.matrix[indlogitDemFA,1]
logitPelFA <-  rep.matrix[indlogitPelFA,1]
logitDemFA.sd <-  rep.matrix[indlogitDemFA,2]
logitPelFA.sd <-  rep.matrix[indlogitPelFA,2]
SAWinter <- rep.matrix[indSAWinter,1]
SAEco <- rep.matrix[indSAEco,1]
SARussian <- rep.matrix[indSARussian,1]
logSAWinter <- rep.matrix[indlogSAWinter,1]
logSAWinter.sd <- rep.matrix[indlogSAWinter,2]
logSAEco <- rep.matrix[indlogSAEco,1]
logSAEco.sd <- rep.matrix[indlogSAEco,2]
logSARussian <- rep.matrix[indlogSARussian,1]
logSARussian.sd <- rep.matrix[indlogSARussian,2]
logSSB <- rep.matrix[indlogSSB,1]
logSSB.sd <- rep.matrix[indlogSSB,2]


logQSurvey1=rep.matrix[rep.rnames=="logQSurvey1",1]
logQSurvey2=parameters$logQSurvey2
logQSurvey3=rep.matrix[rep.rnames=="logQSurvey3",1]


# reformat the triangular population matrix
logTriN <- rep.matrix[indlogTriN,1]
logTriN.std <- rep.matrix[indlogTriN,2]
logTriNmatrix <-  (matrix(logTriN,nrow = (data$nYears+1),byrow = FALSE))
logTriNmatrix.std <- (matrix(logTriN.std,nrow = (data$nYears+1),byrow = FALSE))

#########
# PLOTS #
#########
pdf('SCA_mentella_plots.pdf')

# Population age structure in the last year+1 of the assessment
PredictedNinLastYear=data.frame(Age=data$minAge:(data$maxAge+data$nYears),
                                Npred=exp(logTriNmatrix[data$nYears+1,])/1e6,
                                N05=exp(logTriNmatrix[data$nYears+1,]-2*logTriNmatrix.std[data$nYears+1,])/1e6,
                                N95=exp(logTriNmatrix[data$nYears+1,]+2*logTriNmatrix.std[data$nYears+1,])/1e6)
#quartz("",7,5)
ggplot(data=PredictedNinLastYear,aes(x=Age,y=Npred))+
  geom_bar(stat="identity",fill="gray75")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=1)+
  labs(x='Age (year)',y='Numbers (millions)',title=paste('Numbers-at-age in',data$maxYear+1))


# Numbers-at-age in year 1
NY1=data.frame(Age=data$minAge:data$maxAge,
               NY1=exp(logNY1)/1e6,
               N05=exp(logNY1-2*logNY1.sd)/1e6,
               N95=exp(logNY1+2*logNY1.sd)/1e6)
#quartz("",7,5)
ggplot(data=NY1,aes(x=Age,y=NY1))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=0.1)+
  labs(x='Age (year)',y='Numbers (millions)',title=paste('Numbers-at-age in',data$minYear))


# Numbers-at-age 1 (2y old)
NA1=data.frame(Year=data$minYear:data$maxYear,
               NA1=exp(logNA1)/1e6,
               N05=exp(logNA1-2*logNA1.sd)/1e6,
               N95=exp(logNA1+2*logNA1.sd)/1e6)
#quartz("",7,5)
ggplot(data=NA1,aes(x=Year,y=NA1))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=0.1)+
  labs(x='Year',y='Numbers (millions)',title='Numbers of 2y old')


# SSB
SSB=data.frame(Year=data$minYear:data$maxYear,
               SSB=exp(logSSB)/1e3,
               SSB05=exp(logSSB-2*logSSB.sd)/1e3,
               SSB95=exp(logSSB+2*logSSB.sd)/1e3)
#quartz("",7,5)
ggplot(data=SSB,aes(x=Year,y=SSB))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = SSB05, ymax = SSB95),colour='black',size=.5,fatten=0.1)+
  labs(x='Year',y='Biomass (tonnes)',title='Spawning Stock Biomass')


# Fishing mortalities
FY=data.frame(Year=data$minYear:data$maxYear,
               DemFY=exp(DemlogFY),
               DemFY05=exp(DemlogFY-2*DemlogFY.sd),
               DemFY95=exp(DemlogFY+2*DemlogFY.sd),
               PelFY=exp(PellogFY),
               PelFY05=exp(PellogFY-2*PellogFY.sd),
               PelFY95=exp(PellogFY+2*PellogFY.sd))

#quartz("",7,5)
ggplot(data=FY,aes(x=Year))+
  geom_line(aes(y=DemFY))+
  geom_ribbon(aes(y = DemFY,ymin = DemFY05, ymax = DemFY95),fill='lightblue',alpha=0.5)+
  geom_line(aes(y=PelFY),colour='red')+
  geom_ribbon(aes(y = PelFY,ymin = PelFY05, ymax = PelFY95),fill='lightpink',alpha=0.5)+
  labs(x='Year',y='Fishing mortality (Fys)',title='Annual fishing mortality')


# Fishing selectivities
# ... to be continued, these need to be reported as 'sdreport' in the cpp code of the model
FA=data.frame(Age=data$minAge:data$maxAge,
              DemFA=exp(logitDemFA)/(1+exp(logitDemFA)),
              DemFA05=exp(logitDemFA-2*logitDemFA.sd)/(1+exp(logitDemFA-2*logitDemFA.sd)),
              DemFA95=exp(logitDemFA+2*logitDemFA.sd)/(1+exp(logitDemFA+2*logitDemFA.sd)),
              PelFA=exp(logitPelFA)/(1+exp(logitPelFA)),
              PelFA05=exp(logitPelFA-2*logitPelFA.sd)/(1+exp(logitPelFA-2*logitPelFA.sd)),
              PelFA95=exp(logitPelFA+2*logitPelFA.sd)/(1+exp(logitPelFA+2*logitPelFA.sd))
              )
#quartz("",7,5)
ggplot(data=FA,aes(x=Age))+
  geom_line(aes(y=DemFA),colour='blue')+
  geom_ribbon(aes(y = DemFA,ymin = DemFA05, ymax = DemFA95),fill='lightblue',alpha=0.5)+
  geom_line(aes(y=PelFA),colour='red')+
  geom_ribbon(aes(y = PelFA,ymin = PelFA05, ymax = PelFA95),fill='lightpink',alpha=0.5)+
  labs(x='Age (year)',y='Fleet selectivity (Sa)',title='Fleet selectivities-at-age')+
  lims(y=c(0,1))




# Survey selectivities
# with confidence intervals
SA=data.frame(Age=data$minAge:data$maxAge,
              SAWinter=exp(logSAWinter),
              SAWinter05=exp(logSAWinter-2*logSAWinter.sd),
              SAWinter95=exp(logSAWinter+2*logSAWinter.sd),
              SAEco=exp(logSAEco),
              SAEco05=exp(logSAEco-2*logSAEco.sd),              
              SAEco95=exp(logSAEco+2*logSAEco.sd),              
              SARussian=exp(logSARussian),
              SARussian05=exp(logSARussian-2*logSARussian.sd),              
              SARussian95=exp(logSARussian+2*logSARussian.sd))
#quartz("",7,5)
ggplot(data=SA,aes(x=Age))+
  geom_ribbon(aes(ymin = SAWinter05, ymax = SAWinter95), fill = "lightblue",alpha=0.5)+
  geom_ribbon(aes(ymin = SAEco05, ymax = SAEco95), fill = "gray70",alpha=0.5)+
  geom_ribbon(aes(ymin = SARussian05, ymax = SARussian95), fill = "lightpink",alpha=0.5)+
  geom_line(aes(y=SAWinter),colour='blue')+
  geom_line(aes(y=SAEco),colour='black')+
  geom_line(aes(y=SARussian),colour='red')+
  labs(x='Age (year)',y='Survey selectivity (Sa)',title='Survey selectivities-at-age')+
  lims(y=c(0,1.2))


###########################
# Diagnostics - Residuals #
###########################

## CATCHES

# preliminary calculations:

# Natural mortalities
M=exp(parameters$logM2)
# Fishing mortalities
DemF=FY$DemFY%*%t(FA$DemFA)
PelF=FY$PelFY%*%t(FA$PelFA)
TotF=DemF+PelF
Z=TotF+M

# Abundance Matrix
Abundance=matrix(data=NA,nrow=data$nYears,ncol=18)
row.names(Abundance)=data$YearSpan
colnames(Abundance)=2:19
# filling up first row
Abundance[1,]=exp(logNY1)
# filling up first column
Abundance[,1]=exp(logNA1)

# filling up the rest of the matrix
for(year in 2:data$nYears){
  for (age in 2:17){
    Abundance[year,age]=Abundance[year-1,age-1]*exp(-Z[year-1,age-1])
  }
  Abundance[year,18]=Abundance[year-1,17]*exp(-Z[year-1,17])+Abundance[year-1,18]*exp(-Z[year-1,18])
}

Abundance2=cbind(exp(logTriNmatrix[1:data$nYears,1:17]),
                 data.frame(plus.group=colSums(exp(logTriNmatrix[1:data$nYears,18:(17+data$nYears)]))))
Abundance2=as.matrix(Abundance2)
colnames(Abundance2)=2:19

# Observed Catch Matrix
DemersalCatchAtAge=cast(data=as.data.frame(data$CatchAtAge),formula=Year~Age,value="Demersal")
DemersalCatchAtAge=as.matrix(DemersalCatchAtAge[,2:15])
PelagicCatchAtAge=cast(data=as.data.frame(data$CatchAtAge),formula=Year~Age,value="Pelagic")
PelagicCatchAtAge=as.matrix(PelagicCatchAtAge[,2:15])

# Predicted Catch Matrix
DemlogC=PellogC=matrix(nrow=data$nYears,ncol=14)
for (year in 1:data$nYears){
  for (age in 5:18){
    DemlogC[year,age-4]=log(DemF[year,age])-log(TotF[year,age]+M)+log(1-exp(-TotF[year,age]-M))+log(Abundance[year,age])
    PellogC[year,age-4]=log(PelF[year,age])-log(TotF[year,age]+M)+log(1-exp(-TotF[year,age]-M))+log(Abundance[year,age])
  }
}
DemC=exp(DemlogC)
PelC=exp(PellogC)
PredlogC=log(DemC+PelC)

###########################
# Diagnostics - Residuals #
###########################

##############################
# Catches in the demersal fleet

#quartz()
par(mfrow=c(2,2))
# predicted vs observed catches
obs=DemersalCatchAtAge
pred=DemC
valid=(obs>0)&(pred>0)
obs.vec=log(as.vector(obs[valid]))
pred.vec=log(as.vector(pred[valid]))
Demersal.obs=obs.vec
Demersal.pred=pred.vec

plot(obs.vec~pred.vec,main='Demersal fleet',xlim=c(6,18),ylim=c(6,18),ylab='log(observed catches)',xlab='log(fitted catches)')
abline(0,1,lty=2,col='red')
# Anomalies in log-catches
delta=as.matrix(log(obs)-log(pred))
delta[is.finite(delta)==FALSE]=NA
# Residuals by Age
boxplot(delta,xlab='Age',ylab='catch residuals (log)',names=c('6','7','8','9','10','11','12','13','14','15','16','17','18','19'),main='residuals by age')
abline(0,0,lty=2,col='red')
# Residuals by Year
boxplot(t(delta),xlab='Year',ylab='catch residuals (log)',names=data$minYear:data$maxYear,main='residuals by year',las=2)
abline(0,0,lty=2,col='red')
# Bubble-plot by Age*Year
Xpos=rep(1,data$nYears)%*%t(6:19)
Ypos=(data$minYear:data$maxYear)%*%t(rep(1,14))
plot(Xpos,Ypos,cex=delta^.5,pch=19,col=hsv(0.6,1,1,alpha=0.75),
     xlab='Age',ylab='Year',main='residuals by age & year',
     xlim=c(1.5,19.5),ylim=c(data$minYear-0.5,data$maxYear+0.5))
points(Xpos,Ypos,cex=(-delta)^.5,pch=19,col=hsv(.95,1,1,alpha=0.75))

##############################
# Catches in the pelagic fleet

#quartz()
par(mfrow=c(2,2))
# predicted vs observed catches
obs=PelagicCatchAtAge
pred=PelC
valid=(obs>0)&(pred>0)
obs.vec=log(as.vector(obs[valid]))
pred.vec=log(as.vector(pred[valid]))
Pelagic.obs=obs.vec
Pelagic.pred=pred.vec

plot(obs.vec~pred.vec,main='Pelagic fleet',xlim=c(6,18),ylim=c(6,18),ylab='log(observed catches)',xlab='log(fitted catches)')
abline(0,1,lty=2,col='red')
# Anomalies in log-catches
delta=as.matrix(log(obs)-log(pred))
delta[is.finite(delta)==FALSE]=NA
# Residuals by Age
boxplot(delta,xlab='Age',ylab='catch residuals (log)',names=c('6','7','8','9','10','11','12','13','14','15','16','17','18','19'),main='residuals by age')
abline(0,0,lty=2,col='red')
# Residuals by Year
boxplot(t(delta),xlab='Year',ylab='catch residuals (log)',names=data$minYear:data$maxYear,main='residuals by year',las=2)
abline(0,0,lty=2,col='red')
# Bubble-plot by Age*Year
Xpos=rep(1,data$nYears)%*%t(6:19)
Ypos=(data$minYear:data$maxYear)%*%t(rep(1,14))
plot(Xpos,Ypos,cex=delta^.5,pch=19,col=hsv(0.6,1,1,alpha=0.75),
     xlab='Age',ylab='Year',main='residuals by age & year',
     xlim=c(1.5,19.5),ylim=c(data$minYear-0.5,data$maxYear+0.5))
points(Xpos,Ypos,cex=(-delta)^.5,pch=19,col=hsv(.95,1,1,alpha=0.75))

##############################
## SURVEYS

## WINTER SURVEY

#quartz()
par(mfrow=c(2,2))
# predicted vs observed catches
obs=as.matrix(Winter[1:length(YearSpan),2:15])
Lq=logQSurvey1
LSs=logSAWinter # log-selectivity
pred=obs
for (year in 1:length(YearSpan)){
  for (age in 2:15){
    pred[year,age-1]=exp(Lq+LSs[age-1]+log(Abundance[year,age-1])-0.5*Z[year,age-1])
  }
}
valid=(obs>0)&(pred>0)
obs.vec=log(as.vector(obs[valid]))
pred.vec=log(as.vector(pred[valid]))
Winter.obs=obs.vec
Winter.pred=pred.vec

plot(obs.vec~pred.vec,main='Winter survey',xlim=c(3,15),ylim=c(3,15),ylab='log(observed indices)',xlab='log(fitted indices)')
abline(0,1,lty=2,col='red')
# Anomalies in log-catches
delta=as.matrix(log(obs)-log(pred))
delta[is.finite(delta)==FALSE]=NA
# Residuals by Age
boxplot(delta,xlab='Age',ylab='indices residuals (log)',names=c('2','3','4','5','6','7','8','9','10','11','12','13','14','15'),main='residuals by age')
abline(0,0,lty=2,col='red')
# Residuals by Year
boxplot(t(delta),xlab='Year',ylab='indices residuals (log)',names=YearSpan,main='residuals by year',las=2)
abline(0,0,lty=2,col='red')
# Bubble-plot by Age*Year
Xpos=rep(1,length(YearSpan))%*%t(2:15)
Ypos=YearSpan%*%t(rep(1,14))
plot(Xpos,Ypos,cex=delta^.5,pch=19,col=hsv(0.6,1,1,alpha=0.75),xlab='Age',ylab='Year',main='residuals by age & year',xlim=c(1.5,19.5),ylim=c(min(YearSpan)-0.5,max(YearSpan)+0.5))
points(Xpos,Ypos,cex=(-delta)^.5,pch=19,col=hsv(.95,1,1,alpha=0.75))

## ECOSYSTEM SURVEY

#quartz()
par(mfrow=c(2,2))
# predicted vs observed catches
obs=as.matrix(Ecosystem[1:length(YearSpan),2:15])
Lq=logQSurvey2
LSs=logSAEco # log-selectivity
pred=obs
for (year in 1:length(YearSpan)){
  for (age in 2:15){
    pred[year,age-1]=exp(Lq+LSs[age-1]+log(Abundance[year,age-1])-0.5*Z[year,age-1])
  }
}
valid=(obs>0)&(pred>0)
obs.vec=log(as.vector(obs[valid]))
pred.vec=log(as.vector(pred[valid]))
Ecosystem.obs=obs.vec
Ecosystem.pred=pred.vec

plot(obs.vec~pred.vec,main='Ecosystem survey',xlim=c(3,15),ylim=c(3,15),ylab='log(observed indices)',xlab='log(fitted indices)')
abline(0,1,lty=2,col='red')
# Anomalies in log-catches
delta=as.matrix(log(obs)-log(pred))
delta[is.finite(delta)==FALSE]=NA
# Residuals by Age
boxplot(delta,xlab='Age',ylab='indices residuals (log)',names=c('2','3','4','5','6','7','8','9','10','11','12','13','14','15'),main='residuals by age')
abline(0,0,lty=2,col='red')
# Residuals by Year
boxplot(t(delta),xlab='Year',ylab='indices residuals (log)',names=YearSpan,main='residuals by year',las=2)
abline(0,0,lty=2,col='red')
# Bubble-plot by Age*Year
Xpos=rep(1,length(YearSpan))%*%t(2:15)
Ypos=YearSpan%*%t(rep(1,14))
plot(Xpos,Ypos,cex=delta^.5,pch=19,col=hsv(0.6,1,1,alpha=0.75),xlab='Age',ylab='Year',main='residuals by age & year',xlim=c(1.5,19.5),ylim=c(min(YearSpan)-0.5,max(YearSpan)+0.5))
points(Xpos,Ypos,cex=(-delta)^.5,pch=19,col=hsv(.95,1,1,alpha=0.75))

## RUSSIAN GROUNDFISH SURVEY

#quartz()
par(mfrow=c(2,2))
# predicted vs observed catches
obs=as.matrix(Russian2[1:length(YearSpan),2:15])
Lq=logQSurvey3
LSs=logSARussian # log-selectivity
pred=obs
for (year in 1:length(YearSpan)){
  for (age in 2:15){
    pred[year,age-1]=exp(Lq+LSs[age-1]+log(Abundance[year,age-1])-0.5*Z[year,age-1])
  }
}
valid=(obs>0)&(pred>0)
obs.vec=log(as.vector(obs[valid]))
pred.vec=log(as.vector(pred[valid]))
Russian.obs=obs.vec
Russian.pred=pred.vec

plot(obs.vec~pred.vec,main='Russian Groundfish survey',xlim=c(-4,4),ylim=c(-4,4),ylab='log(observed indices)',xlab='log(fitted indices)')
abline(0,1,lty=2,col='red')
# Anomalies in log-catches
delta=as.matrix(log(obs)-log(pred))
delta[is.finite(delta)==FALSE]=NA
# Residuals by Age
boxplot(delta,xlab='Age',ylab='indices residuals (log)',names=c('2','3','4','5','6','7','8','9','10','11','12','13','14','15'),main='residuals by age')
abline(0,0,lty=2,col='red')
# Residuals by Year
boxplot(t(delta),xlab='Year',ylab='indices residuals (log)',names=YearSpan,main='residuals by year',las=2)
abline(0,0,lty=2,col='red')
# Bubble-plot by Age*Year
Xpos=rep(1,length(YearSpan))%*%t(2:15)
Ypos=YearSpan%*%t(rep(1,14))
plot(Xpos,Ypos,cex=delta^.5,pch=19,col=hsv(0.6,1,1,alpha=0.75),xlab='Age',ylab='Year',main='residuals by age & year',xlim=c(1.5,19.5),ylim=c(min(YearSpan)-0.5,max(YearSpan)+0.5))
points(Xpos,Ypos,cex=(-delta)^.5,pch=19,col=hsv(.95,1,1,alpha=0.75))

# Statistical distributions of residuals
library(gamlss)
par(mfrow=c(3,2))
histDist((Demersal.obs-Demersal.pred),'NO',nbins=25,main='Demersal Catches')
histDist((Pelagic.obs-Pelagic.pred),'NO',nbins=10,main='Pelagic Catches')
histDist((Winter.obs-Winter.pred),'NO',nbins=25,main='Winter Survey')
histDist((Ecosystem.obs-Ecosystem.pred),'NO',nbins=25,main='Ecosystem Survey')
histDist((Russian.obs-Russian.pred),'NO',nbins=25,main='Russian Groundfish Survey')

par(mfrow=c(3,2))
qqnorm((Demersal.obs-Demersal.pred),main='Demersal Catches')
qqnorm((Pelagic.obs-Pelagic.pred),main='Pelagic Catches')
qqnorm((Winter.obs-Winter.pred),main='Winter Survey')
qqnorm((Ecosystem.obs-Ecosystem.pred),main='Ecosystem Survey')
qqnorm((Russian.obs-Russian.pred),main='Russian Groundfish Survey')

dev.off()
