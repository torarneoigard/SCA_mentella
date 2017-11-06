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
# Benjamin Planque, February 2016, September 2017

# data preparation --------------------------------------------------------
load('SCA_mentella_model.Rdata')            # Load model results from TMB
data=model$data                             # extract data
YearSpan=data$minYear:data$maxYear
parameters=model$parameters                 # and parameters

rep.matrix <- summary(model$rep)            # get the list of reported values and standard deviations
rep.rnames <- rownames(rep.matrix)          # get the names of variables
indlogNY1 <- which(rep.rnames=="logNY1")    # extract line numbers for numbers in year one
indlogNA1 <- which(rep.rnames=="logNA1")    # extract line numbers for numbers at age one
indDemlogFY <- which(rep.rnames=="DemlogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
indPellogFY <- which(rep.rnames=="PellogFY") # extract line numbers for demersal fishing mortality in years (Fy's)
indlogitDemFA <- which(rep.rnames=="logitDemFA") # extract line numbers for demersal fishing mortality in years (Fy's)
indlogitPelFA <- which(rep.rnames=="logitPelFA") # extract line numbers for demersal fishing mortality in years (Fy's)
indSA <- which(rep.rnames=="SA")            # extract line numbers for survey selectivities at age (SA's)
indlogSA <- which(rep.rnames=="logSA")      # extract line numbers for survey log-selectivities at age (SA's)
indSAProp <- which(rep.rnames=="SAProp")    # extract line numbers for proportion-survey selectivities at age (SA's)
indlogSSB <- which(rep.rnames=="logSSB")    # extract line numbers for SSB
indlogTriN <- which(rep.rnames=="logTriN")  # extract line numbers for triangular population matrix
indIndexPropTruncMa <- which(rep.rnames=="IndexPropTruncMa")

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
SA <- matrix(rep.matrix[indSA,1],nrow = length(surveys),byrow = FALSE)
logSA <- matrix(rep.matrix[indlogSA,1],nrow = length(surveys),byrow = FALSE)
logSA.sd <- matrix(rep.matrix[indlogSA,2],nrow = length(surveys),byrow = FALSE)
SAProp <- matrix(rep.matrix[indSAProp,1],nrow = length(surveysProp),byrow = FALSE)

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

#########
# PLOTS #
#########
pdf.filename=paste('SCA_mentella_plots',model$date.flag,'.pdf',sep="")
pdf(pdf.filename)

# Population age structure ------------------------------------------------
# Population age structure in the last year+1 of the assessment
PredictedNinLastYear=data.frame(Age=data$minAge:(data$maxAge+data$nYears),
                                Npred=exp(logTriNmatrix[data$nYears+1,])/1e6,
                                N05=exp(logTriNmatrix[data$nYears+1,]-2*logTriNmatrix.std[data$nYears+1,])/1e6,
                                N95=exp(logTriNmatrix[data$nYears+1,]+2*logTriNmatrix.std[data$nYears+1,])/1e6)
#quartz("",7,5)
print(ggplot(data=PredictedNinLastYear,aes(x=Age,y=Npred))+
  geom_bar(stat="identity",fill="gray75")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=1)+
  labs(x='Age (year)',y='Numbers (millions)',title=paste('Numbers-at-age in',data$maxYear+1)))

# Numbers-at-age in year 1
NY1=data.frame(Age=data$minAge:data$maxAge,
               NY1=exp(logNY1)/1e6,
               N05=exp(logNY1-2*logNY1.sd)/1e6,
               N95=exp(logNY1+2*logNY1.sd)/1e6)
#quartz("",7,5)
print(ggplot(data=NY1,aes(x=Age,y=NY1))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=0.1)+
  labs(x='Age (year)',y='Numbers (millions)',title=paste('Numbers-at-age in',data$minYear)))


# Numbers-at-age 1 (2y old)
NA1=data.frame(Year=data$minYear:data$maxYear,
               NA1=exp(logNA1)/1e6,
               N05=exp(logNA1-2*logNA1.sd)/1e6,
               N95=exp(logNA1+2*logNA1.sd)/1e6)
#quartz("",7,5)
print(ggplot(data=NA1,aes(x=Year,y=NA1))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = N05, ymax = N95),colour='black',size=.5,fatten=0.1)+
  labs(x='Year',y='Numbers (millions)',title='Numbers of 2y old'))


# SSB
SSB=data.frame(Year=data$minYear:data$maxYear,
               SSB=exp(logSSB)/1e3,
               SSB05=exp(logSSB-2*logSSB.sd)/1e3,
               SSB95=exp(logSSB+2*logSSB.sd)/1e3)
#quartz("",7,5)
print(ggplot(data=SSB,aes(x=Year,y=SSB))+
  geom_bar(stat="identity",fill="gray80")+
  geom_pointrange(aes(ymin = SSB05, ymax = SSB95),colour='black',size=.5,fatten=0.1)+
  labs(x='Year',y='Biomass (tonnes)',title='Spawning Stock Biomass'))


# Fishing mortalities -----------------------------------------------------
FY=data.frame(Year=data$minYear:data$maxYear,
               DemFY=exp(DemlogFY),
               DemFY05=exp(DemlogFY-2*DemlogFY.sd),
               DemFY95=exp(DemlogFY+2*DemlogFY.sd),
               PelFY=exp(PellogFY),
               PelFY05=exp(PellogFY-2*PellogFY.sd),
               PelFY95=exp(PellogFY+2*PellogFY.sd))


#quartz("",7,5)
print(ggplot(data=FY,aes(x=Year))+
  geom_line(aes(y=DemFY))+
  geom_ribbon(aes(y = DemFY,ymin = DemFY05, ymax = DemFY95),fill='lightblue',alpha=0.5)+
  geom_line(aes(y=PelFY),colour='red')+
  geom_ribbon(aes(y = PelFY,ymin = PelFY05, ymax = PelFY95),fill='lightpink',alpha=0.5)+
  labs(x='Year',y='Fishing mortality (Fys)',title='Annual fishing mortality'))


# Fishing selectivities
FA=data.frame(Age=data$minAge:data$maxAge,
              DemFA=exp(logitDemFA)/(1+exp(logitDemFA)),
              DemFA05=exp(logitDemFA-2*logitDemFA.sd)/(1+exp(logitDemFA-2*logitDemFA.sd)),
              DemFA95=exp(logitDemFA+2*logitDemFA.sd)/(1+exp(logitDemFA+2*logitDemFA.sd)),
              PelFA=exp(logitPelFA)/(1+exp(logitPelFA)),
              PelFA05=exp(logitPelFA-2*logitPelFA.sd)/(1+exp(logitPelFA-2*logitPelFA.sd)),
              PelFA95=exp(logitPelFA+2*logitPelFA.sd)/(1+exp(logitPelFA+2*logitPelFA.sd))
              )
#quartz("",7,5)
print(ggplot(data=FA,aes(x=Age))+
  geom_line(aes(y=DemFA),colour='blue')+
  geom_ribbon(aes(y = DemFA,ymin = DemFA05, ymax = DemFA95),fill='lightblue',alpha=0.5)+
  geom_line(aes(y=PelFA),colour='red')+
  geom_ribbon(aes(y = PelFA,ymin = PelFA05, ymax = PelFA95),fill='lightpink',alpha=0.5)+
  labs(x='Age (year)',y='Fleet selectivity (Sa)',title='Fleet selectivities-at-age')+
  lims(y=c(0,1)))


# Survey selectivities ----------------------------------------------------
# with confidence intervals
SAdf <- list(Age=data$minAge:data$maxAge)
SAdf.names <- c("Age")
counter = 2
for (i in 1:length(surveys)){
  SAdf[[counter]] <- exp(logSA[i,])
  SAdf.names <- append(SAdf.names,paste("SA",surveys[i],sep = ""))
  SAdf[[counter+1]] <- exp(logSA[i,]-2*logSA.sd[i,])
  SAdf.names <- append(SAdf.names,paste("SA",surveys[i],"05",sep = ""))
  SAdf[[counter+2]] <- exp(logSA[i,]+2*logSA.sd[i,])
  SAdf.names <- append(SAdf.names,paste("SA",surveys[i],"95",sep = ""))
  counter <- counter + 3
}
names(SAdf) <- SAdf.names
SAdf <- as.data.frame(SAdf)

ggplot(data=SAdf,aes(x=Age))+
  geom_ribbon(aes(ymin = SAWinter05, ymax = SAWinter95), fill = "lightblue",alpha=0.5)+
  geom_ribbon(aes(ymin = SAEcosystem05, ymax = SAEcosystem95), fill = "gray70",alpha=0.5)+
  geom_ribbon(aes(ymin = SARussian05, ymax = SARussian95), fill = "lightpink",alpha=0.5)+
  geom_line(aes(y=SAWinter),colour='blue')+
  geom_line(aes(y=SAEcosystem),colour='black')+
  geom_line(aes(y=SARussian),colour='red')+
  labs(x='Age (year)',y='Survey selectivity (Sa)',title='Survey selectivities-at-age')+
  lims(y=c(0,1.2))



## CATCHES

# Diagnostics residuals ---------------------------------------------------
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


# Demersal fleet#############################
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

# Pelagic fleet#############################
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

# Surveys#############################
## SURVEYS
obs.vec <- list()
pred.vec <- list()
for (i in 1:length(surveys)){
  par(mfrow=c(2,2))
  # predicted vs observed catches
  indDF <- which(names(data)== surveys[i])
  obs=as.matrix(data[[indDF]][1:length(YearSpan),2:15])
  Lq=logQSurvey[i]
  LSs=logSA[i,] # log-selectivity
  pred=obs
  for (year in 1:length(YearSpan)){
    for (age in 2:15){
      pred[year,age-1]=exp(Lq+LSs[age-1]+log(Abundance[year,age-1])-0.5*Z[year,age-1])
    }
  }
  valid=(obs>0)&(pred>0)
  obs.vec[[i]]=log(as.vector(obs[valid]))
  pred.vec[[i]]=log(as.vector(pred[valid]))
  #Winter.obs=obs.vec
  #Winter.pred=pred.vec
  
  plot(obs.vec[[i]]~pred.vec[[i]],main=paste(surveys[i],"survey"),
       xlim=c(min(pred.vec[[i]],na.rm=TRUE)*0.8,max(pred.vec[[i]],na.rm=TRUE)*1.2),
       ylim=c(min(obs.vec[[i]],na.rm=TRUE)*0.8,max(obs.vec[[i]],na.rm=TRUE)*1.2),
       ylab='log(observed indices)',xlab='log(fitted indices)')
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
}

# Statistical distributions of residuals
library(gamlss)
par(mfrow=c(3,2))
histDist((Demersal.obs-Demersal.pred),'NO',nbins=25,main='Demersal Catches')
histDist((Pelagic.obs-Pelagic.pred),'NO',nbins=10,main='Pelagic Catches')

for (i in 1:length(surveys)){
  histDist((obs.vec[[i]]-pred.vec[[i]]),'NO',nbins=25,main=paste(surveys[i],"survey"))
}
#histDist((Winter.obs-Winter.pred),'NO',nbins=25,main='Winter Survey')
#histDist((Ecosystem.obs-Ecosystem.pred),'NO',nbins=25,main='Ecosystem Survey')
#histDist((Russian.obs-Russian.pred),'NO',nbins=25,main='Russian Groundfish Survey')

par(mfrow=c(3,2))
qqnorm((Demersal.obs-Demersal.pred),main='Demersal Catches')
qqnorm((Pelagic.obs-Pelagic.pred),main='Pelagic Catches')

for (i in 1:length(surveys)){
  qqnorm((obs.vec[[i]]-pred.vec[[i]]),main=paste(surveys[i],"survey"))
}

# Surveys in proportion ---------------------------------------------------
if(length(surveysProp)>0){
  # Calculate expected numbers at age based on numbers in the population and survey selectivity at age
  nAgesInTriMatrix=dim(logTriNmatrix)[2]
  YearSpanInTriMatrix=c(YearSpan,tail(YearSpan,1)+1)
  if("WGIDEEPS"%in%surveysProp){
    Year=rep(YearSpanInTriMatrix,nAgesInTriMatrix) # recode the data matrix into a 4 column table
    Age=rep(data$minAge:(data$minAge+nAgesInTriMatrix-1),each = length(YearSpanInTriMatrix)) # fill in the 'Age' vector
    Survey=rep(1,nAgesInTriMatrix*length(YearSpanInTriMatrix)) # fill in the 'Survey' vector
    TriNmatrix=exp(logTriNmatrix);TriNmatrix[TriNmatrix==1]=0  # get population numbers from the triangular matrix
    SAProp.df=as.data.frame(SAProp) # format survey selectivity into a data.frame
    SelectivityMatrix=SAProp.df[rep(seq_len(nrow(SAProp.df)), each=length(YearSpanInTriMatrix)),] # Construct a Matrix of survey selectivity (years*ages)
    SelectivityMatrix[,(dim(SAProp)[2]+1):(nAgesInTriMatrix)]=1 # set selectivity for 19+ ages to 1
    TriNmatrix2=TriNmatrix*SelectivityMatrix # compute the Numbers-at-age  corrected by the selectivity
    
    # Recode the triangular matrix in proportions and age-blocks instead of ages
    Prop <- NULL
    AgeBlocks=as.data.frame(data$AgeBlocks) # format age-blocks into a data.frame
    Age2Blocks=function(XProp,AgeBlocks,StartYear,PlusGroupInStartYear){ # recode a dataset by age in years into a dataset by age in age blocks
      # Xprop is the original table with data by age
      # AgeBlocks is the table that contains the definition of the age blocks to use
      # StartYear is the first year to run the model (1992)
      # PlusGroupInStartYear is what is says (19y)
      Years=unique(XProp$Year)
      Surveys=unique(XProp$Survey) 
      XPropBlocks={}
      for (S in Surveys){ # loop on surveys
        for (Y in Years){ # loop on Years
          for (A in AgeBlocks[,1]){ # loop on Age blocks
            if (AgeBlocks$maxAge[A]<PlusGroupInStartYear+Y-StartYear){ # test if the age block does not contains the +group for that year
              Xsub=subset(XProp,(Survey==S)&(Year==Y)&(Age>=AgeBlocks$minAge[A])&(Age<=AgeBlocks$maxAge[A]))
              if(dim(Xsub)[1]>0){
                XPB=data.frame(Year=Y,AgeBlock=A,Survey=S,IndexProp=sum(Xsub$IndexProp)) #
                XPropBlocks=rbind(XPropBlocks,XPB)
              }
            } else { # +group case
              if (AgeBlocks$minAge[A]<=PlusGroupInStartYear+Y-StartYear){
                Xsub=subset(XProp,(Survey==S)&(Year==Y)&(Age>=AgeBlocks$minAge[A]))
                if(dim(Xsub)[1]>0){
                  XPB=data.frame(Year=Y,AgeBlock=A,Survey=S,IndexProp=sum(Xsub$IndexProp)) #
                  XPropBlocks=rbind(XPropBlocks,XPB)
                }
              }
            }
          }
        }
      }
      XPropBlocks
    } # function to recode age into age-blocks
    PropMatrix=TriNmatrix2/rowSums(TriNmatrix2) # calculate proportions
    IndexProp=as.numeric(as.matrix(PropMatrix)) # fill in the 'Index' vector
    Prop=rbind(Prop,data.frame(Year=Year,Age=Age,Survey=Survey,IndexProp=IndexProp)) 
    PropBlocks=Age2Blocks(Prop,AgeBlocks,YearSpan[1],19) # construct the same table but with AgeBlocks instead of Ages in years
  }
}
SurveyProps=as.data.frame(model$data$SurveyProps) # observed proportions-at-age block in the survey(s)
plot.new()
print(ggplot()+
  geom_point(data = PropBlocks,aes(y=Year,x=(AgeBlocks$minAge[PropBlocks$AgeBlock]+AgeBlocks$maxAge[PropBlocks$AgeBlock])/2,size=IndexProp*100),shape = 21, fill = "lightblue",alpha=0.5)+
  geom_point(data = SurveyProps,aes(y=Year,x=(AgeBlocks$minAge[SurveyProps$AgeBlock]+AgeBlocks$maxAge[SurveyProps$AgeBlock])/2,size=IndexProp*100),shape = 21, fill = "lightpink3",alpha=0.5)+
  scale_size_area(name='Proportion (%)')+
  labs(x='Age (year)',y='Year',title='Modelled vs. Observed proportions-at-age'))
  
# ADDITIONAL PLOTS AND TABLES############################
# 
#############################

# STOCK SUMMARY PLOT AND TABLE

# POP MATRIX
N=exp(logTriNmatrix[-dim(logTriNmatrix)[1],]) # numbers-at-age (until last year of data)
N[N==1]=0 # replace ones by zeros (log<->exp transformation issue)

# WEIGHT MATRIX
W=matrix(nrow=dim(N)[1],ncol=dim(N)[2])
W[,1:(data$maxAge-data$minAge+1)]=data$WeightAtAge
W[,(data$maxAge-data$minAge+2):dim(W)[2]]=data$WeightAtAge[,(data$maxAge-data$minAge+1)]%*%t(rep(1,dim(W)[2]-(data$maxAge-data$minAge+1)))

# Maturity matrix
Mat=matrix(nrow=dim(N)[1],ncol=dim(N)[2])
Mat[,1:(data$maxAge-data$minAge+1)]=data$MaturityAtAge
Mat[,(data$maxAge-data$minAge+2):dim(Mat)[2]]=data$MaturityAtAge[,(data$maxAge-data$minAge+1)]%*%t(rep(1,dim(Mat)[2]-(data$maxAge-data$minAge+1)))

# Biomass
TSB=rowSums(N*W)/1000 # Total Stock Biomass in tonnes
SSB=rowSums(N*W*Mat)/1000 # Spawning Stock Biomass in tonnes

# Recruits
Rec2=exp(logTriNmatrix[1:length(YearSpan),(2-data$minAge+1)])/1e6 # in millions
Rec6=exp(logTriNmatrix[1:length(YearSpan),(6-data$minAge+1)])/1e6 # in millions

# Fishing mortality
DemFY=exp(DemlogFY)
DemFA=exp(logitDemFA)/(1+exp(logitDemFA))
PelFY=exp(PellogFY)
PelFA=exp(logitPelFA)/(1+exp(logitPelFA))
Ftot=(DemFY)%*%t(DemFA)+(PelFY)%*%t(PelFA)
F12.18=rowMeans(Ftot[,((12:18)-data$minAge+1)])
F19=Ftot[,(19-data$minAge+1)]

# Stock Summary Table
SS=data.frame(Year=YearSpan,
              Recruits.at.age.2=round(Rec2),
              Recruits.at.age.6=round(Rec6),
              Total.Stock.Biomass=round(TSB),
              Spawning.Stock.Biomass=round(SSB),
              F12.18=round(F12.18,3),
              F19=round(F19,3))
write.table(SS,'SCA_mentella_stock_summary.txt',row.names = FALSE,quote = FALSE,sep='\t')

# Stock Summary table extended
SSE=data.frame(Year=YearSpan,
              Rec.age.2=round(exp(logNA1)),
              Rec.age.2.025=round(exp(logNA1-logNA1.sd*1.96)),
              Rec.age.2.975=round(exp(logNA1+logNA1.sd*1.96)),
              SSB=round(exp(logSSB)),
              SSB.025=round(exp(logSSB-logSSB.sd*1.96)),
              SSB.975=round(exp(logSSB+logSSB.sd*1.96)),
              TSB=round(TSB),
              Catches.in.Tonnes=round(data$TotalCatches[,2]),
              Catches.Pelagic=round(data$TotalCatches[,3]),
              Catches.Other=round(data$TotalCatches[,4]),
              F12.18=round(F12.18,3),
              F19=round(F19,3))
write.table(SSE,'SCA_mentella_stock_summary4graphs.txt',row.names = FALSE,quote = FALSE,sep='\t')

plot.new()
ggplot_overlay=function(p1,p2){
  # extract gtable
  g1 <- ggplot_gtable(ggplot_build(p1))
  g2 <- ggplot_gtable(ggplot_build(p2))
  
  # overlap the panel of 2nd plot on that of 1st plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                       pp$l, pp$b, pp$l)
  
  # axis tweaks
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0, "cm")
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
  g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
  
  # draw it
  grid.draw(g)
}

# plot stock summary
p1=ggplot(data=SS,aes(x=Year))+
  geom_ribbon(aes(ymin = 0, ymax = Total.Stock.Biomass/1000),fill='lightblue',colour='black',alpha=0.7)+
  geom_ribbon(aes(ymin = 0, ymax = Spawning.Stock.Biomass/1000),fill='darkblue',colour='black',alpha=0.7)+
  labs(x='Year',y='Biomass (1000 tonnes)',title='S. mentella in ICES subareas 1 and 2 - summary')+
  xlim(model$data$minYear-0.5,model$data$maxYear+.5)+
  theme_bw()+
  theme(panel.grid=element_blank(),panel.spacing=unit(rep(0,4),rep("pt",4)))

p2=ggplot(data=SS, aes(x=Year, y=Recruits.at.age.2))+ 
  geom_bar(stat="identity",fill="yellow",colour='black',alpha=0.5,width=0.6)+
  labs(y='Recruits at age 2y (millions)')+
  xlim(model$data$minYear-0.5,model$data$maxYear+.5)+
  theme_bw()+
  theme(panel.background = element_rect(fill = NA),panel.grid=element_blank(),panel.border=element_blank())

print(ggplot_overlay(p1,p2))
dev.off()

#######
# plot age-structure in last year of the assessment
#######
LastYear=data.frame(Age=model$data$minAge+(1:dim(N)[2])-1,
                    N=N[dim(N)[1],],
                    TSB=N[dim(N)[1],]*W[dim(N)[1],]/1e6,
                    SSB=N[dim(N)[1],]*W[dim(N)[1],]*Mat[dim(N)[1],]/1e6)
LastYear=LastYear[-dim(LastYear)[1],]
#scaling.coef=1.1*max(LastYear$TSB)/max(LastYear$N)
#LastYear$tN=LastYear$N*scaling.coef
p1=ggplot(data=LastYear,aes(x=Age))+
  geom_ribbon(aes(ymin = 0, ymax = TSB),fill='lightblue',colour='black',alpha=0.7)+
  geom_ribbon(aes(ymin = 0, ymax = SSB),fill='darkblue',colour='black',alpha=0.7)+
#  geom_bar(stat='identity',aes(y=tN),width=.5,fill='yellow',colour='black',alpha=0.5)+
  labs(x='Age',y='Biomass (1000 tonnes)',title=paste('Age structure in ',model$data$maxYear,sep=''))+
  xlim(data$minAge-0.5,data$minAge+dim(N)[2]-1+0.5)+
  theme_bw()+
  theme(panel.grid=element_blank())

p2=ggplot(data=LastYear, aes(x=Age, y=N/1e6))+ 
  geom_bar(stat="identity",fill="yellow",colour='black',alpha=0.5,width=0.6)+
  labs(y='Numbers (millions)')+
  xlim(data$minAge-0.5,data$minAge+dim(N)[2]-1+0.5)+
  theme_bw()+
  theme(panel.background = element_rect(fill = NA),panel.grid=element_blank(),panel.border=element_blank())

pdf.filename2=paste('SCA_mentella_plots',model$date.flag,'2.pdf',sep="")
pdf(pdf.filename2)
print(ggplot_overlay(p1,p2))
dev.off()

