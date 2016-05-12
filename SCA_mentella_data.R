# SCA_mentella.R in TMB
# Program designed to perform a statistical catch-at-age assessment based upon
# catch data (total catch and pelagic fleet catches) as reported to ICES and
# survey data (0-group, winter, ecosyste and Russian groundfish surveys),
#
# The present R file is doing the data prepapration 
# the file SCA_mentella_model.R is runing the SCA in TMB (together with SCA_mentella_model.cpp)
# the file SCA_mentella_plots.R is plotting the results
#
# The version in the directory "~/Documents/Work/Redfish/ICES/AFWG2015/Mentella/SCAA_mentella_tmb" reproduce the SCA outputs from the AFWG2014
# Benjamin Planque, March 2015, updated February 2016


#####################
## CATCH DATA      ##
#####################
## Total Catch in tonnes 
TotalCatches=read.delim('TotalCatches.txt') # total catches in tonnes as reported in ICES AFWG (Table 6.1)
# selection of data within the year span
TotalCatches=subset(TotalCatches,Year%in%YearSpan)

## Loading total catch-at-age data
TotalCatchAtAge=read.table("CatchAtAge.txt",header=TRUE) # note that the last group is a +group
# selection of data within the year span
TotalCatchAtAge=subset(TotalCatchAtAge,Year%in%YearSpan)

# reshaping into a three column vector with Year, Age, Catch in number
# modification by Alf 20.04.2016

TCA = TotalCatchAtAge                           # helping variable
siz = dim(TCA)                                 # size of TCA, = [nrow ncol]
ny = siz[1]                                     # no of ages
na = siz[2]-1                                   # no or ages, first column contains years
age1 = as.numeric(substr(colnames(TCA)[2],2,5)) # first age
AgeSpan = age1:(age1+na-1)                      # parenthesis necessary!
yrvec = as.matrix(rep(YearSpan,na))             # repeats yrspan na time in col.vector
agevec = as.matrix(rep(AgeSpan,each = ny))      # repeats each age ny times in col.vector
catchvec = array(as.matrix(TCA[,2:(na+1)]))     # makes columnvector of catchmatrix for col3 in new table
TotalCatchAtAge2 = cbind(yrvec,agevec,catchvec*1000) # new table
#

# The above modification replaces the following commands
#
#TotalCatchAtAge2=matrix(nrow=prod(dim(TotalCatchAtAge)[1],dim(TotalCatchAtAge)[2]-1),ncol=3)
#StartAge=as.numeric(substr(colnames(TotalCatchAtAge)[2],2,5))
#for (age in 2:dim(TotalCatchAtAge)[2])
#{
#  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),1]=TotalCatchAtAge[,1] # Year
#  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),2]=StartAge+age-2 # Age
#  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),3]=TotalCatchAtAge[1:dim(TotalCatchAtAge)[1],age]*1000 # Numbers
#}
# end of removed section due to modification

colnames(TotalCatchAtAge2)=c("Year","Age","Catch.in.number")
TotalCatchAtAge2=as.data.frame(TotalCatchAtAge2)
# selection of data within the year span
TotalCatchAtAge2=subset(TotalCatchAtAge2,Year%in%YearSpan)

## Loading pelagic fleet catch-at-age data
PelagicCatchAtAge=read.table("PelagicCatchAtAge.txt",header=TRUE) # note that the last group is a +group
# selection of data within the year span
PelagicCatchAtAge=subset(PelagicCatchAtAge,Year%in%YearSpan)

# reshaping into a three column vector with Year, Age, Catch in number
PelagicCatchAtAge2=matrix(nrow=prod(dim(PelagicCatchAtAge)[1],dim(PelagicCatchAtAge)[2]-1),ncol=3)
PelagicStartAge=as.numeric(substr(colnames(PelagicCatchAtAge)[2],2,5))
for (age in 2:dim(PelagicCatchAtAge)[2])
{
  PelagicCatchAtAge2[((age-2)*(dim(PelagicCatchAtAge)[1])+1):((age-1)*(dim(PelagicCatchAtAge)[1])),1]=PelagicCatchAtAge[,1] # Year
  PelagicCatchAtAge2[((age-2)*(dim(PelagicCatchAtAge)[1])+1):((age-1)*(dim(PelagicCatchAtAge)[1])),2]=PelagicStartAge+age-2 # Age
  PelagicCatchAtAge2[((age-2)*(dim(PelagicCatchAtAge)[1])+1):((age-1)*(dim(PelagicCatchAtAge)[1])),3]=PelagicCatchAtAge[1:dim(PelagicCatchAtAge)[1],age]*1000 # Numbers
}
colnames(PelagicCatchAtAge2)=c("Year","Age","Catch.in.number")
PelagicCatchAtAge2=as.data.frame(PelagicCatchAtAge2)
# selection of data within the year span
PelagicCatchAtAge2=subset(PelagicCatchAtAge2,Year%in%YearSpan)

DemersalCatchAtAge=TotalCatchAtAge-PelagicCatchAtAge
DemersalCatchAtAge$Year=TotalCatchAtAge$Year
# selection of data within the year span
DemersalCatchAtAge=subset(DemersalCatchAtAge,Year%in%YearSpan)

# Assembly of the total and pelagic catches into one data frame
CatchAtAge=cbind(TotalCatchAtAge2,PelagicCatchAtAge2$Catch.in.number,TotalCatchAtAge2$Catch.in.number-PelagicCatchAtAge2$Catch.in.number)
colnames(CatchAtAge)=c("Year","Age","Total","Pelagic","Demersal")

# selection of positive total catch only
CatchAtAge2=subset(CatchAtAge,Total>0)

#####################
## WEIGHT DATA     ##
#####################
## Loading weight data
WeightAtAge=read.table("WeightAtAgeInCatch.txt",header=TRUE) # note that the last group is a +group
# selection of data within the year span
WeightAtAge=subset(WeightAtAge,Year%in%YearSpan)
# selection of ages 2 to 19+
WeightAtAge=WeightAtAge[,4:21]

#####################
## MATURITYT DATA  ##
#####################
## Loading maturity data
MaturityAtAge=read.table("MaturityAtAge.txt",header=TRUE) # note that the last group is a +group
# selection of data within the year span
MaturityAtAge=subset(MaturityAtAge,Year%in%YearSpan)
# selection of ages 2 to 19+
MaturityAtAge=MaturityAtAge[,4:21]
#####################
## SURVEY DATA     ##
#####################

logQSurveyInit <- NULL
logQSurveyMap <- NULL

X <- NULL
SurveyTime <- NULL
pa0Init <- NULL
logb1Init <- NULL
logb2Init <- NULL
logb1Map <- NULL
logb2Map <- NULL

lowerAgeBoundary <- NULL
upperAgeBoundary <- NULL

## Reshaping data into a 5 column vector with
## Year, Age, Survey, Gear selectivity, Survey Index
surveyCounter <- 1

if("Winter"%in%surveys){
  ## Loading survey data files
  Winter=read.table("WinterSurvey.txt",header=TRUE) # survey indices 1992-2010 for ages 2-15 (YC: 1977-2008)
  # selection of data within the year span

  Winter=subset(Winter,Year%in%YearSpan)

  ## Reshaping data into a 5 column vector with
  ## Year, Age, Survey, Gear selectivity, Survey Index
  
  Year=rep(Winter$Year,dim(Winter)[2]-1)
  #Age=as.vector(rep(1,dim(Winter)[1])%*%t(2:15))
  Age<-rep(2:15,each = dim(Winter)[1])
  Survey=rep(surveyCounter,dim(Winter)[1]*(dim(Winter)[2]-1))
  Index=as.numeric(as.matrix(Winter[,2:15]))
  
  X=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)
  #Xa=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)
  SurveyTime <- append(SurveyTime,0.12)
  logQSurveyInit <- append(logQSurveyInit,-8)
  logQSurveyMap <- append(logQSurveyMap,1)
  
  pa0Init <- append(pa0Init,log((7.0887-2)/(15-7.0887)))
  
  logb1Init <- append(logb1Init,-10)
  logb1Map <- append(logb1Map,NA)
  
  logb2Init <- append(logb2Init,-2.0311)
  logb2Map <- append(logb2Map,1)
  
  lowerAgeBoundary <- append(lowerAgeBoundary,2)
  upperAgeBoundary <- append(upperAgeBoundary,13)
  surveyCounter <- surveyCounter + 1
}

if("Ecosystem"%in%surveys){
  Ecosystem=read.table("EcosystemSurvey.txt",header=TRUE) # survey indices 1996-2009 for ages 2-15 (YC: 1981-2007)
  # selection of data within the year span
  Ecosystem=subset(Ecosystem,Year%in%YearSpan)
  
  ## Reshaping data into a 5 column vector with
  ## Year, Age, Survey, Gear selectivity, Survey Index
  Year=rep(Ecosystem$Year,dim(Ecosystem)[2]-1)
  #Age=as.vector(rep(1,dim(Ecosystem)[1])%*%t(2:15))
  Age=rep(2:15,each = dim(Ecosystem)[1])
  Survey=rep(surveyCounter,dim(Ecosystem)[1]*(dim(Ecosystem)[2]-1))
  Index=as.numeric(as.matrix(Ecosystem[,2:15]))
  
  #Xb <- data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)
  X <- rbind(X,data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index))
  SurveyTime <- append(SurveyTime,0.75)
  logQSurveyInit <- append(logQSurveyInit,-8.160)
  logQSurveyMap <- append(logQSurveyMap,NA)
  
  pa0Init <- append(pa0Init,log((8.5883-2)/(15-8.5883)))
  
  logb1Init <- append(logb1Init,-10)
  logb1Map <- append(logb1Map,NA)
  
  logb2Init <- append(logb2Init,-0.97676)
  logb2Map <- append(logb2Map,2)
  
  lowerAgeBoundary <- append(lowerAgeBoundary,2)
  upperAgeBoundary <- append(upperAgeBoundary,13)
  surveyCounter <- surveyCounter + 1
  
}

if("Russian"%in%surveys){
  
  Russian=read.table("RussianSurvey.txt",header=TRUE) # Year-class indices 1974-2008
  # preparation of Russian groundfish survey data
  Russian2=matrix(nrow=dim(Russian)[1],ncol=dim(Russian)[2])
  StartAge=2
  Russian2[,1]=Russian[,1]+StartAge
  for (age in 2:dim(Russian)[2])
  {
    Russian2[(age-1):dim(Russian)[1],age]=Russian[1:(dim(Russian)[1]-age+2),age]
  }
  colnames(Russian2)=colnames(Russian)
  colnames(Russian2)[1]="Year"
  Russian2=as.data.frame(Russian2)
  # selection of data within the year span
  Russian2=subset(Russian2,Year%in%YearSpan)
  
  Year=rep(Russian2$Year,dim(Russian2)[2]-1)
  #Age=as.vector(rep(1,dim(Russian2)[1])%*%t(2:15))
  Age=rep(2:15,each = dim(Russian2)[1])
  Survey=rep(surveyCounter,dim(Russian2)[1]*(dim(Russian2)[2]-1))
  Index=as.numeric(as.matrix(Russian2[,2:15]))
  
  #Xc=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)
  
  X=rbind(X,data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index))
  
  SurveyTime <- append(SurveyTime,0.9)
  logQSurveyInit <- append(logQSurveyInit,-16)
  logQSurveyMap <- append(logQSurveyMap,2)
  pa0Init <- append(pa0Init,log((8.0891-2)/(11-8.0891)))
  logb1Init <- append(logb1Init,-0.15861)
  logb1Map <- append(logb1Map,1)
  logb2Init <- append(logb2Init,-10)
  logb2Map <- append(logb2Map,NA)
  
  lowerAgeBoundary <- append(lowerAgeBoundary,2)
  upperAgeBoundary <- append(upperAgeBoundary,9)
  surveyCounter <- surveyCounter + 1
  
}

if("WGIDEEPS"%in%surveys){
  WGIDEEPS=read.table("WGIDEEPS.txt",header=TRUE) # survey indices 2008,2009,2013 for ages 7-75
  # selection of data within the year span
  WGIDEEPS=subset(WGIDEEPS,Year%in%YearSpan)
  
  Year=rep(WGIDEEPS$Year,dim(WGIDEEPS)[2]-1)
  Age=rep(7:75,each = dim(WGIDEEPS)[1])
  Survey=rep(surveyCounter,dim(WGIDEEPS)[1]*(dim(WGIDEEPS)[2]-1))
  Index=as.numeric(as.matrix(WGIDEEPS[,2:70]))
  
  X=rbind(X,data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index))
  
  SurveyTime <- append(SurveyTime,0.67)
  
  #THESE NEED TO BE SET FOR THIS SURVEY
  logQSurveyInit <- append(logQSurveyInit,-16)
  logQSurveyMap <- append(logQSurveyMap,2)
  pa0Init <- append(pa0Init,log((8.0891-2)/(11-8.0891)))
  logb1Init <- append(logb1Init,-0.15861)
  logb1Map <- append(logb1Map,1)
  logb2Init <- append(logb2Init,-10)
  logb2Map <- append(logb2Map,NA)
  
  lowerAgeBoundary <- append(lowerAgeBoundary,2)
  upperAgeBoundary <- append(upperAgeBoundary,9)
  
  surveyCounter <- surveyCounter + 1
  
}
  

#X=rbind(Xa,Xb,Xc,Xd) # combining data from the different surveys

# remove lines with zero data
X=subset(X,Index>0)
# selection of data within the year span
X=subset(X,Year%in%YearSpan)

# Survey Timing (as fraction of the year)
#SurveyTime=c(0.12,0.75,0.9,0.67)

####################################
## EXPORT DATA  TO BE USED IN TMB ##
####################################

data=list()
#data$Date=date()
data$REswitch=REswitch
data$minYear=YearSpan[1]
data$maxYear=YearSpan[length(YearSpan)]
data$minAgeInCatch=min(CatchAtAge2$Age)
data$maxAgeInCatch=max(CatchAtAge2$Age)
data$minAgeInSurvey=min(X$Age)
data$maxAgeInSurvey=max(X$Age)
data$minAge=min(data$minAgeInCatch,data$minAgeInSurvey)
#data$maxAge=max(data$maxAgeInCatch,data$maxAgeInSurvey) # needs revision
data$maxAge=19 # needs revision
data$TotalCatches=as.matrix(TotalCatches)
data$CatchAtAge=as.matrix(CatchAtAge2)
data$CatchNrow=dim(data$CatchAtAge)[1]
data$SurveyIndex=as.matrix(X)
data$SurveyNrow=dim(data$SurveyIndex)[1]
data$WeightAtAge=as.matrix(WeightAtAge)
data$MaturityAtAge=as.matrix(MaturityAtAge)
data$nYears=length(data$minYear:data$maxYear)
data$nYearsPel=length(2006:data$maxYear)
data$nAgesInCatch=length(data$minAgeInCatch:data$maxAgeInCatch)
data$nAges=length(data$minAge:data$maxAge)
data$nSurveys=length(unique(X$Survey))
data$SurveyTime=SurveyTime
data$lowerAgeBoundary = lowerAgeBoundary
data$upperAgeBoundary = upperAgeBoundary
data$pa0Init <- pa0Init
data$logb1Init <- logb1Init
data$logb2Init <- logb2Init
# Additional data needed for the plots - THIS NEED TO BE MADE MORE DYNAMIC
if("Winter"%in%surveys) data$Winter=Winter
if("Ecosystem"%in%surveys) data$Ecosystem=Ecosystem
if("Russian"%in%surveys) data$Russian=Russian2
if("WGIDEEPS"%in%surveys) data$WGIDEEPS=WGIDEEPS

data$logQSurveyInit <- logQSurveyInit      #Initial value for logQSurveys
data$logQSurveyMap <- factor(logQSurveyMap)
data$logb1Map <- factor(logb1Map)
data$logb2Map <- factor(logb2Map)

save(data,file='SCA_mentella_data.Rdata')

