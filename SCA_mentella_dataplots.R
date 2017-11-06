# SCA_mentella_dataplots.R in TMB
# Program designed to produce plots of the
# data used in the catch-at-age assessment
#
# Benjamin Planque, April 2016, October 2017

graphics.off()
load('SCA_mentella_data.Rdata')            # load model results
YearSpan=data$minYear:data$maxYear         # extract year span

pdf(file='SCA_mentella_dataplots.pdf')     # open file to save all plots

# total catches
TotalCatches=as.data.frame(data$TotalCatches)
print(ggplot(data=TotalCatches,aes(x=Year))+
  geom_ribbon(aes(ymin=0,ymax=Catch),fill='orangered3',colour='black')+
  geom_ribbon(aes(ymin=0,ymax=Other),fill='royalblue4',colour='black')+
  labs(y='Catch in tonnes'))
  
# catches-at-age
X=as.data.frame(data$CatchAtAge)
X1=data.frame(Year=X$Year,Age=X$Age,Fleet='Pelagic',Catch=X$Pelagic/1e6)
X1=subset(X1,Catch>0)
X2=data.frame(Year=X$Year,Age=X$Age,Fleet='Demersal',Catch=X$Demersal/1e6)
X2=subset(X2,Catch>0)
CatchAtAge=rbind(X1,X2)
CatchAtAge$Age=as.factor(CatchAtAge$Age)
CatchAtAge$Fleet=as.factor(CatchAtAge$Fleet)

# catch-at-age total numbers
nages=nlevels(CatchAtAge$Age)
CatchAtAge$Age=factor(CatchAtAge$Age,levels=levels(CatchAtAge$Age)[nages:1])
print(ggplot(data=CatchAtAge,aes(x=Year))+
  facet_grid(Fleet~.)+
  theme(panel.grid=element_blank())+
  geom_bar(stat='identity',aes(y=Catch,fill=Age),width=1,position='stack',colour='black')+
  #scale_fill_discrete()+
  scale_fill_grey(start = 1, end = 0.25)+
  xlim(data$minYear-0.5,data$maxYear+0.5)+
  labs(y='Catch-at-age (millions)'))

# catch-at-age proportions
print(ggplot(data=CatchAtAge,aes(x=Year))+
  facet_grid(Fleet~.,scales='free_y')+
  theme(panel.grid=element_blank())+
  geom_bar(stat='identity',aes(y=Catch,fill=Age),width=1,position='fill',colour='black')+
  #scale_fill_discrete()+
  scale_fill_grey(start = 1, end = 0.25)+
  xlim(data$minYear-0.5,data$maxYear+0.5)+
  labs(y='Catch-at-age (%)'))

SurveyIndex=as.data.frame(data$SurveyIndex)
SurveyIndex$Age=as.factor(SurveyIndex$Age)
SurveyIndex$Survey=as.factor(SurveyIndex$Survey)

#SurveyIndex$Survey=mapvalues(SurveyIndex$Survey, from = c("1", "2", "3","4"), to = c("Winter", "Ecosystem", "Russian","WGIDEEPS"))
SurveyIndex$Survey=mapvalues(SurveyIndex$Survey, from = as.character(1:length(surveys)), to = surveys)

# survey indices by age, total numbers
nages=nlevels(SurveyIndex$Age)
SurveyIndex$Age=factor(SurveyIndex$Age,levels=levels(SurveyIndex$Age)[nages:1])
ggplot(data=SurveyIndex,aes(x=Year))+
  facet_grid(Survey~.,scales='free_y')+
  theme(panel.grid=element_blank())+
  geom_bar(stat='identity',aes(y=Index,fill=Age),width=1,position='stack',colour='black')+
  #scale_fill_discrete()+
  scale_fill_grey(start = 1, end = 0.25)+
  xlim(data$minYear-0.5,data$maxYear+0.5)

# survey indices by age, proportions
ggplot(data=SurveyIndex,aes(x=Year))+
  facet_grid(Survey~.,scales='free_y')+
  theme(panel.grid=element_blank())+
  geom_bar(stat='identity',aes(y=Index,fill=Age),width=1,position='fill',colour='black')+
  #scale_fill_discrete()+
  scale_fill_grey(start = 1, end = 0.25)+
  xlim(data$minYear-0.5,data$maxYear+0.5)

# survey indices in proportion (e.g. WGIDEEPS) !!! For now this code only works with one survey in proportions
if(PropSurveySwitch==1){ # test that data is provided for the surveys in proportions
  # Age blocks used in the model for data in proportion
  AgeBlocks=as.data.frame(data$AgeBlocks)
  AgeBlocks$AgeBlocks=(AgeBlocks$maxAge-AgeBlocks$minAge+1)
  AgeBlocks$valid=1
  AgeBlocks=rbind(AgeBlocks[1,],AgeBlocks)
  AgeBlocks$AgeBlocks[1]=AgeBlocks$minAge[1]
  AgeBlocks$valid[1]=0
  nageblocks=dim(AgeBlocks)[1]
  SurveyProps=as.data.frame(data$SurveyProps)
  minageblock=min(SurveyProps$AgeBlock)
  maxageblock=max(SurveyProps$AgeBlock)
  SurveyProps$AgeBlock=as.factor(SurveyProps$AgeBlock)
  SurveyProps$Survey=as.factor(SurveyProps$Survey)
  SurveyProps$Survey=mapvalues(SurveyProps$Survey, from = as.character(1:length(surveysProp)), to = surveysProp)
  SurveyProps$AgeBlock=factor(SurveyProps$AgeBlock,levels=levels(SurveyProps$AgeBlock)[nageblocks:1])
  Labels=paste(as.character(AgeBlocks$minAge[maxageblock]),"+",sep='')
  Labels=c(Labels,paste(as.character(AgeBlocks$minAge[(maxageblock-1):minageblock]),as.character(AgeBlocks$maxAge[(maxageblock-1):minageblock]),sep ='-'))
  
  print(ggplot(data=SurveyProps,aes(x=Year))+
    facet_grid(Survey~.,scales='free_y')+
    theme(panel.grid=element_blank())+
    geom_bar(stat='identity',aes(y=IndexProp,fill=AgeBlock),width=1,position='stack',colour='black')+
    scale_fill_grey(start = 1, end = 0.25,labels=Labels)+
    #scale_fill_grey(start = 1, end = 0.25)+
    xlim(data$minYear-0.5,data$maxYear+0.5))
  
  # Age blocks used in the model for data in proportion
  print(ggplot(data=AgeBlocks)+
    geom_bar(stat='identity',aes(x=0,y=AgeBlocks,fill=valid),width=1,position='stack',colour='black',show.legend=FALSE)+
    theme(aspect.ratio=.1)+
    coord_flip())
} 



# Reshaping of Maturity and Weight-at-age data
M=as.vector(t(data$MaturityAtAge))
W=as.vector(t(data$WeightAtAge))
Y=rep(YearSpan,each=(data$maxAge-data$minAge+1))
A=rep(data$minAge:data$maxAge,times=length(YearSpan))
AtAge=data.frame(Year=Y,Age=A,Maturity=M,Weight=W)

# Maturity-at-age
# ggplot(data=AtAge,aes(x=Age,y=Maturity))+
#   facet_wrap(~Year,nrow=4)+
#   geom_line()+
#   geom_abline(slope=0,intercept=0.5,colour='red')

ggplot(data=AtAge,aes(x=Year,y=Maturity))+
  facet_wrap(~Age,nrow=3)+
  geom_line()+
  theme(axis.text.x=element_text(angle=90))

# Weight-at-age
# ggplot(data=AtAge,aes(x=Age,y=Weight))+
#   facet_wrap(~Year,nrow=4)+
#   geom_point()

ggplot(data=AtAge,aes(x=Year,y=Weight))+
  facet_wrap(~Age,nrow=3)+
  geom_line()+
  theme(axis.text.x=element_text(angle=90))

dev.off()
