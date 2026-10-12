#########################################################################
## Malaria vector ODE model                                            ##
## John Marshall (john.marshall@berkeley.edu)                          ##
## 30/May/2015                                                         ##
## Adapted and extended by SK starting Feb 2016                        ##
#########################################################################
## Coding vector model from White et al. (2011) as a system of ODEs    ##
## Adapted by SK on Feb 08, 2016 for extension                         ##
#########################################################################
## Quick-run script with the default parameters (Driver.R is the full   ##
## driver). It draws the same figures as the Shiny app.                 ##
#########################################################################

#------------------------------------------------------------------------
########### LOAD LIBRARIES AND SOURCES ##################################
rm(list=ls())
library(deSolve)
library(ggplot2)
source("ODEModel.R")
source("ODEMosquitoParameters.R")
source("Plots.R")
source("ODEControlMeasuresParameters.R")
source("ODETransmissionParameters.R")
source("ODEInterventions.R")
source("ODEModelOutput.R")
source("FeedingCycle.R")

########### INIT MODEL ##################################################
theta = getTheta()
initState = calculateInitialState(theta)

########### RUN THE IVM ODEs ############################################
IVM_traj = runODE(360,1,initState,theta,"lsoda")

########### PLOT RESULTS ################################################
barChartMosquitoDemographicsShiny(IVM_traj)
plotTrajectoryShiny(IVM_traj)
plotEIRShiny(IVM_traj)
#------------------------------------------------------------------------
