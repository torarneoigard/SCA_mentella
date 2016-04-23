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
TotalCatchAtAge2=matrix(nrow=prod(dim(TotalCatchAtAge)[1],dim(TotalCatchAtAge)[2]-1),ncol=3)
StartAge=as.numeric(substr(colnames(TotalCatchAtAge)[2],2,5))
for (age in 2:dim(TotalCatchAtAge)[2])
{
  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),1]=TotalCatchAtAge[,1] # Year
  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),2]=StartAge+age-2 # Age
  TotalCatchAtAge2[((age-2)*(dim(TotalCatchAtAge)[1])+1):((age-1)*(dim(TotalCatchAtAge)[1])),3]=TotalCatchAtAge[1:dim(TotalCatchAtAge)[1],age]*1000 # Numbers
}
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
## Loading survey data files
Winter=read.table("WinterSurvey.txt",header=TRUE) # survey indices 1992-2010 for ages 2-15 (YC: 1977-2008)
# selection of data within the year span
Winter=subset(Winter,Year%in%YearSpan)

Ecosystem=read.table("EcosystemSurvey.txt",header=TRUE) # survey indices 1996-2009 for ages 2-15 (YC: 1981-2007)
# selection of data within the year span
Ecosystem=subset(Ecosystem,Year%in%YearSpan)

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

## Reshaping data into a 5 column vector with
## Year, Age, Survey, Gear selectivity, Survey Index

# a. Winter
Year=rep(Winter$Year,dim(Winter)[2]-1)
Age=as.vector(rep(1,dim(Winter)[1])%*%t(2:15))
Survey=rep(1,dim(Winter)[1]*(dim(Winter)[2]-1))
Index=as.numeric(as.matrix(Winter[,2:15]))

Xa=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)

# b. Ecosystem
Year=rep(Ecosystem$Year,dim(Ecosystem)[2]-1)
Age=as.vector(rep(1,dim(Ecosystem)[1])%*%t(2:15))
Survey=rep(2,dim(Ecosystem)[1]*(dim(Ecosystem)[2]-1))
Index=as.numeric(as.matrix(Ecosystem[,2:15]))

Xb=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)

# c. Russian groundfish
Year=rep(Russian2$Year,dim(Russian2)[2]-1)
Age=as.vector(rep(1,dim(Russian2)[1])%*%t(2:15))
Survey=rep(3,dim(Russian2)[1]*(dim(Russian2)[2]-1))
Index=as.numeric(as.matrix(Russian2[,2:15]))

Xc=data.frame(Year=Year,Age=Age,Survey=Survey,Index=Index)

X=rbind(Xa,Xb,Xc) # the Ogroup survey is no longer included

# remove lines with zero data
X=subset(X,Index>0)
# selection of data within the year span
X=subset(X,Year%in%YearSpan)

# Survey Timing
SurveyTime=c(0.12,0.75,0.9)

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
data$maxAge=max(data$maxAgeInCatch,data$maxAgeInSurvey)
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
# Additional data needed for the plots
data$Winter=Winter
data$Ecosystem=Ecosystem
data$Russian2=Russian2

save(data,file='SCA_mentella_data.Rdata')

