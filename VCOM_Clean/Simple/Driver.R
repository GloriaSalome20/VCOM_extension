###############################################################################
# Driver.R  --  VCOM driver script
#
#   Sources the model files, sets up the scenario (species, run time,
#   intervention coverages and start times), runs the ODE model and draws
#   the standard plots.
#
#   How to use:  1. set the working directory to this folder
#                2. edit section 1 (scenario) and, if needed, the species in section 2;
#                   the climate model is chosen in section 3
#                3. run the script (or source("Driver.R"))
#
#   Files used (all in this folder):
#     ODEModel.R                      ODE model, initial state, solver wrapper
#     ODEMosquitoParameters.R         species / life-cycle parameters and getTheta()
#     ISIMIP3b_climate.R              daily rainfall, temperature and humidity of Tanga from ISIMIP3b climate models
#                                     (reads isimip3b_tanga_climatology.csv, made by buildClimatology(), section 3)
#     ODETemperatureParameters.R      temperature- and humidity-dependent parameters (only if switched on, section 5)
#     ODEControlMeasuresParameters.R  intervention coverage and start times
#     ODETransmissionParameters.R     transmission parameters
#     ODEInterventions.R              impact of each intervention
#     FeedingCycle.R                  feeding-cycle impact of the interventions
#     ODEModelOutput.R                EIR, vectorial capacity, R0, ...
#     Plots.R, multiplot.R            figures
#     Parsers.R                       import of parameter files
#     KL_gamma_temp.R                 seasonal carrying capacity K(t) (only if the seasonal K is switched on, section 4;
#                                     made from the rainfall and temperature of the ISIMIP3b climatology and written
#                                     to carrying_capacity_tanga_gamma_isimip_<model>.csv)
###############################################################################

# 0. SETUP -------------------------------------------------------------------

# Clear all stored objects
rm(list=ls())

library(shiny)
library(deSolve)
library(ggplot2)
#library(gridGraphics)   # Disable this after testing to avoid messing up Shiny
#library(shinyjs)
#library(shinythemes)

source("ODEModel.R")
source("ODEMosquitoParameters.R")
source("ODETemperatureParameters.R")
source("ODEControlMeasuresParameters.R")
source("ODETransmissionParameters.R")
source("ODEInterventions.R")
source("FeedingCycle.R")
source("ODEModelOutput.R")
source("Plots.R")
source("multiplot.R")
source("Parsers.R")


# 1. SCENARIO ----------------------------------------------------------------

##NOTE: Different start time for LAR and BIO will crash at the max time entered - to be fixed

# Simulation length: runs per day up to this end time (days)
INITIAL_MODELRUNTIME_VALUE   = 1765

# Burn-in: years the model is run before day 0, without interventions, so that the scenario starts in the model's own
# seasonal cycle and not at the constant-environment equilibrium, which is far too high on 1 January (see ISIMIP_CHANGES.md).
# The interventions still start at the days set below, counted from day 0. The burn-in years are dropped from the output,
# and model day 0 is still 1 January. 0 = no burn-in (start at the equilibrium). One year is enough: the starting state
# is the same for 1 to 4 years.
BURN_IN_YEARS = 1
BURN_IN_DAYS  = 365*BURN_IN_YEARS

# Mosquito species parameters (the one used is chosen in section 2)
MOSQUITO_PARAMETERS_1 = getAnGambiaeParameters()
MOSQUITO_PARAMETERS_2 = getAnArabiensisParameters()
MOSQUITO_PARAMETERS_3 = getAnFunestusParameters()

# Interventions: enter the coverage and the day the intervention is switched on.
#   To turn it ON  - enter the required time, which must be < the max model run time
#   To turn it OFF - enter INITIAL_MODELRUNTIME_VALUE + 1

# Source reduction
INITIAL_SRE_COVERAGE = 0.0
INITIAL_SRE_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# Larviciding
INITIAL_LAR_COVERAGE = 0.7
INITIAL_LAR_TIME     = 20

# Biological control
INITIAL_BIO_COVERAGE = 1
INITIAL_BIO_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# ATSB
INITIAL_ATSB_COVERAGE = .50
INITIAL_ATSB_TIME     = 20

# Space spraying
INITIAL_SSP_COVERAGE = .00
INITIAL_SSP_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# Odor traps
INITIAL_OBT_COVERAGE = .00
INITIAL_OBT_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# LLINs
INITIAL_ITN_COVERAGE = .50
INITIAL_ITN_TIME     = 20

# IRS
INITIAL_IRS_COVERAGE = 0.8
INITIAL_IRS_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# House modification
INITIAL_HOU_COVERAGE = 0.0
INITIAL_HOU_TIME     = 20

# Personal protection measure
INITIAL_PPM_COVERAGE = 0.90
INITIAL_PPM_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# Cattle - systemic
INITIAL_ECS_COVERAGE = 0.50
INITIAL_ECS_TIME     = 20

# Cattle - topical
INITIAL_ECT_COVERAGE = 0.50
INITIAL_ECT_TIME     = INITIAL_MODELRUNTIME_VALUE + 1

# Resting and ovipositing - ovitraps (assuming same coverage for ATSB and SSP)
INITIAL_OVI_COVERAGE = 0.80
INITIAL_OVI_TIME     = INITIAL_MODELRUNTIME_VALUE + 1


# 2. PARAMETERS (theta) ---------------------------------------------------------

## Intervention parameters
INTERVENTION_PARAMETERS = getInterventionsParameters(
  SREcov=INITIAL_SRE_COVERAGE,   time_SRE_on=INITIAL_SRE_TIME + BURN_IN_DAYS,    # Source reduction
  LARcov=INITIAL_LAR_COVERAGE,   time_LAR_on=INITIAL_LAR_TIME + BURN_IN_DAYS,    # Larviciding
  BIOcov=INITIAL_BIO_COVERAGE,   time_BIO_on=INITIAL_BIO_TIME + BURN_IN_DAYS,    # Biological control
  ATSBcov=INITIAL_ATSB_COVERAGE, time_ATSB_on=INITIAL_ATSB_TIME + BURN_IN_DAYS,  # ATSB
  SSPcov=INITIAL_SSP_COVERAGE,   time_SSP_on=INITIAL_SSP_TIME + BURN_IN_DAYS,    # Space spraying
  OBTcov=INITIAL_OBT_COVERAGE,   time_OBT_on=INITIAL_OBT_TIME + BURN_IN_DAYS,    # Odor traps
  ITNcov=INITIAL_ITN_COVERAGE,   time_ITN_on=INITIAL_ITN_TIME + BURN_IN_DAYS,    # LLINs
  IRScov=INITIAL_IRS_COVERAGE,   time_IRS_on=INITIAL_IRS_TIME + BURN_IN_DAYS,    # IRS
  HOUcov=INITIAL_HOU_COVERAGE,   time_HOU_on=INITIAL_HOU_TIME + BURN_IN_DAYS,    # House modification
  PPMcov=INITIAL_PPM_COVERAGE,   time_PPM_on=INITIAL_PPM_TIME + BURN_IN_DAYS,    # Personal protection measure
  ECScov=INITIAL_ECS_COVERAGE,   time_ECS_on=INITIAL_ECS_TIME + BURN_IN_DAYS,    # Cattle - systemic
  ECTcov=INITIAL_ECT_COVERAGE,   time_ECT_on=INITIAL_ECT_TIME + BURN_IN_DAYS,    # Cattle - topical
  OVIcov=INITIAL_OVI_COVERAGE,   time_OVI_on=INITIAL_OVI_TIME + BURN_IN_DAYS)    # Resting & ovipositing - ovitraps (same for ATSB, SSP)

## Pass in the intervention parameters (with the updated coverage) and the species.
## Species used: MOSQUITO_PARAMETERS_1 = An. gambiae, _2 = An. arabiensis, _3 = An. funestus
theta <<- getTheta(interventionParameters=INTERVENTION_PARAMETERS, speciesSpecificParameters=MOSQUITO_PARAMETERS_3)


# 3. CLIMATE (ISIMIP3b) ---------------------------------------------------------

## Rainfall, temperature and humidity come from the bias-adjusted output of a global climate model (ISIMIP3b, historical
## run) for the Tanga grid cell, instead of the CRU monthly temperature and rain.csv. isimip3b_tanga_climatology.csv is
## made once (the daily values of the cell come from ISIMIP3b_get_tanga.py, the table from buildClimatology() in
## ISIMIP3b_climate.R, see ISIMIP_CHANGES.md): for each day of the year the average over 1995-2014, smoothed over 15 days,
## 365 days, Day 1 = 1 January = model day 0.
ISIMIP_MODEL = "ensemble"    # "ensemble": mean of the five models, or one model:
                             # "GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"
source("ISIMIP3b_climate.R")
CLIMATE <- loadClimate(ISIMIP_MODEL)        # data frame: day, tas (deg C), pr (mm/day), hurs (%)
describeClimate(CLIMATE, ISIMIP_MODEL)


# 4. SEASONAL LARVAL CARRYING CAPACITY -----------------------------------------

## The larval carrying capacity K can follow rainfall (with a lag) and temperature instead of staying constant.
## KL_gamma_temp.R makes K(t) (gamma-weighted rainfall memory x Briere temperature response) from the rainfall and
## temperature of CLIMATE (section 3), and writes it to K_FILE: 365 daily values, where Day 1 is model day 0 and the
## year repeats every 365 days.
SEASONAL_CARRYING_CAPACITY = TRUE    # FALSE: constant K, as in the cleaned model without seasonality
RAIN_PADDING = "circular"            # rain before Day 1 for the rainfall memory: "circular" = the end of the same cycle (the year repeats),
                                     #   "zero" = none, as in the CRU version (see KL_gamma_temp.R, section 4)
K_FILE  = paste0("carrying_capacity_tanga_gamma_isimip_", climateTag(ISIMIP_MODEL), if (RAIN_PADDING == "zero") "_zeropad", ".csv")   # one file per climate model
K_SCALE = "relative"                 # "relative": K(t) = K0 * K_file(t)/max(K_file); K0 is the equilibrium K of the model,
                                     #             so the peak of K(t) equals K0 and the shape comes from the file
                                     # "absolute": K(t) = K_file(t) exactly as written (its maximum is L_eq in KL_gamma_temp.R)

if (SEASONAL_CARRYING_CAPACITY) {
  ## K_FILE is made again when it does not exist or the climate table is newer (delete K_FILE after editing KL_gamma_temp.R)
  if (!file.exists(K_FILE) || file.mtime(K_FILE) < file.mtime(ISIMIP_CLIMATE_FILE)) source("KL_gamma_temp.R")
  carryingCapacity <<- makeSeasonalCarryingCapacity(read.csv(K_FILE)$K, K_SCALE)   # picked up by IVM_ode()
}


# 5. TEMPERATURE- AND HUMIDITY-DEPENDENT PARAMETERS ----------------------------

## Some life-cycle parameters can follow the daily temperature of CLIMATE (ODETemperatureParameters.R). List the groups to
## switch on; character(0) keeps all parameters constant, as in the cleaned model.
##   "adult_mortality"    muV                 adult death rate (and the humidity modifier S(h) if USE_HUMIDITY)
##   "biting_rate"        tau1, tau2, f0      biting rate b(T): the gonotrophic cycle tau1 + tau2 becomes 1/b(T)
##   "eip"                durEV               extrinsic incubation period n(T) = 111/(T - 16)
##   "aquatic_mortality"  muEL, muLL, muPL    egg, larva and pupa mortality
##   "development"        durEL, durLL, durPL egg, larva and pupa development time
##   "oviposition"        beta                eggs per female per day
TEMPERATURE_PARAMETERS = c("adult_mortality", "biting_rate", "eip", "aquatic_mortality", "development", "oviposition")
USE_HUMIDITY  = TRUE     # TRUE: the humidity modifier S(h) uses hurs of CLIMATE (S = 1 above 42 %, so it only acts on days drier than that);
                         # FALSE: humidity not limiting
T_REFERENCE   = NULL     # deg C at which the model starts at its equilibrium; NULL = the temperature on day 0

if (length(TEMPERATURE_PARAMETERS) > 0) {
  temperatureYear <- makeDailyCycle(CLIMATE$tas)          # the same 365 days as K(t): Day 1 = model day 0, repeating
  if (USE_HUMIDITY && all(is.na(CLIMATE$hurs))) message("no humidity (hurs) in the climate table: humidity is not limiting")
  humidityYear    <- if (USE_HUMIDITY && !all(is.na(CLIMATE$hurs))) makeDailyCycle(CLIMATE$hurs) else NULL
  ## The equilibrium (initial state, NV_eq, lambdaV, K0) uses the parameters at the reference temperature ...
  T_start <- if (is.null(T_REFERENCE)) temperatureYear(0) else T_REFERENCE
  theta <<- temperatureDependentTheta(theta, T_start, if (is.null(humidityYear)) NULL else humidityYear(0), TEMPERATURE_PARAMETERS)
  ## ... and during the run IVM_ode() takes the parameters at the temperature of each day
  timeVaryingParameters <<- makeTemperatureDependentParameters(temperatureYear, humidityYear, TEMPERATURE_PARAMETERS)
}


# 6. RUN THE MODEL -------------------------------------------------------------

## Initialize the model
## (<<- is intentional: calculateInitialState() adds lambdaV and NV_eq to the global theta)
initState <<- calculateInitialState(theta)

## Run the model
## To compare scenarios (1 intervention, 2 interventions, ...) call runODE() again with a different theta, with the same
## run time (INITIAL_MODELRUNTIME_VALUE + BURN_IN_DAYS) and the burn-in rows dropped as below
IVM_traj <- runODE(INITIAL_MODELRUNTIME_VALUE + BURN_IN_DAYS,1,initState,theta,"daspk")
IVM_traj <- IVM_traj[IVM_traj$time >= BURN_IN_DAYS, ]         # drop the burn-in years: model day 0 = 1 January
IVM_traj$time <- IVM_traj$time - BURN_IN_DAYS
rownames(IVM_traj) <- NULL
IVM_traj <<- IVM_traj


# 7. OUTPUTS -------------------------------------------------------------------

## Histogram - mosquitoes in each life stage
barChartMosquitoDemographics(IVM_traj)

## ggplot - vectors
plotTrajectory(IVM_traj)

## ggplot - humans
plotTrajectoryHumans(IVM_traj)

## ggplot - EIR only
plotEIR(IVM_traj)

## Optional plots
#plotTrajectoryALL(IVM_traj)         # all mosquito compartments
#plotTrajectoryMalesMosq(IVM_traj)   # male mosquitoes
#plotEIR_VC_R0(IVM_traj)             # compute and plot EIR, VC, R0


## Notes ##########################################################################
# TODO (Sam): check the impact of mosquito-proofed housing
# To debug, put  browser()  where you want the code to stop.
