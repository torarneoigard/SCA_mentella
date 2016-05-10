# SCA_mentella_Do_It_All

#####################
## INITIALISATION  ##
#####################

## LIBRARIES
require(TMB)                                # Load TMB library
require(ggplot2)                            # Load ggplot2 library
require(gtable)                             # Load gtable library
require(reshape)                            # load the reshape library
require(gamlss)                             # load the gamlss library

## CLEAR ALL
cat("\014")                                 # cear console
rm(list=ls())                               # clear workspace
graphics.off()                              # clear graphical windows

## OPTIONS
YearSpan=1992:2015                          # set the range of year over which the model is run
REswitch=1                                  # Switch for running the model with random effects on the recruits (NA1)
                                            # 0=fixed effects, 1=random effects
## Which surveys to include
surveys <- c("Winter","Ecosystem","Russian")
#surveys <- c("Ecosystem","Russian")


## Working Directory
#setwd("~/Documents/Work/Redfish/ICES/AFWG2016/mentella/SCA_mentella_TMB_2016/CurrentVersion") #<- to be adjusted to individual machines

#####################
## DATA PREP       ##
#####################
source('./SCA_mentella_data.R')

#####################
## DATA PLOTS       ##
#####################
source('./SCA_mentella_dataplots.R',print.eval = TRUE)

#####################
## MODEL RUN       ##
#####################
source('./SCA_mentella_model.R')

#####################
## MODEL PLOTS     ##
#####################
source('./SCA_mentella_plots.R',print.eval = TRUE)
