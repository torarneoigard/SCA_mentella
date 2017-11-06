# SCA_mentella_Do_It_All

#####################
## INITIALISATION  ##
#####################

## LIBRARIES
require(TMB)                                # Load TMB library
require(plyr)                               # Load plyr library
require(tidyr)                              # Load tidyr library
require(ggplot2)                            # Load ggplot2 library
require(gtable)                             # Load gtable library
require(reshape)                            # load the reshape library
require(gamlss)                             # load the gamlss library
require(grid)

## CLEAR ALL
cat("\014")                                 # cear console
rm(list=ls())                               # clear workspace
graphics.off()                              # clear graphical windows

YearSpan=1992:2016                          # set the range of year over which the model is run

## OPTIONS
RElogNA=1                                   # Recruits (NA1): 0=fixed effect, 1=random effect
REDemFishMort = 1                           # Demersal fleet fishing mortality: 0=fixed effect, 1=random effect
REPelFishMort = 1                           # Pelagic fleet fishing mortality: 0=fixed effect, 1=random effect
REDemFishSel = 1                            # Demersal fishing selctivity: 0=fixed effect, 1=random effect
PropSurveySwitch = 1                        # Turn on and off surveys with proportions 0 is off and 1 is on


#
# Which surveys to include
surveys <- c("Winter","Ecosystem","Russian")
#surveys <- c("Winter","Ecosystem")
#surveys <- c("Ecosystem","Russian")
#surveys <- c("Ecosystem")

surveysProp <- c("WGIDEEPS")

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

#####################
## PROJECTION TABLE##
#####################
source('./SCA_mentella_projections.R',print.eval = TRUE)
