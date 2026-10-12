########################################################################
# ODEMosquitoParameters.R
# Builds the parameter vector `theta` used by the ODE model from:
#   1. species-specific parameters   (An. gambiae / arabiensis / funestus)
#   2. mosquito life-cycle parameters
#   3. intervention parameters       (see ODEControlMeasuresParameters.R)
#   4. transmission parameters       (see ODETransmissionParameters.R)
########################################################################

######################################################################################
#------------------------------------------------------------------------
############ THETA FACADE ###############################################
#------------------------------------------------------------------------
# Names (and order) of the parameters that make up theta, grouped by source.
LIFECYCLE_PARAMETER_NAMES <- c(
  "beta", "muEL", "muLL", "muPL", "durEL", "durLL", "durPL", "durEV", "gamma", "tau1", "tau2"
)
INTERVENTION_COVERAGE_NAMES <- c(
  "ITNcov", "IRScov", "ECScov", "ECTcov", "HOUcov", "SPRcov", "PPMcov",
  "OBTcov", "SREcov", "BIOcov", "LARcov", "ATSBcov", "SSPcov", "OVIcov"
)
INTERVENTION_TIME_NAMES <- c(
  "time_ITN_on", "time_IRS_on", "time_HOU_on", "time_ECS_on", "time_ECT_on", "time_SPR_on", "time_PPM_on",
  "time_OBT_on", "time_BIO_on", "time_SRE_on", "time_LAR_on", "time_OVI_on", "time_ATSB_on", "time_SSP_on"
)
SPECIES_PARAMETER_NAMES <- c(
  "muV", "Q0", "phiB", "phiI", "rOVI", "sOVI", "rFOG", "sFOG", "rLAR", "sLAR", "rODO", "sODO",
  "rBIO", "sBIO", "rSRE", "sSRE", "rITN", "sITN", "rIRS", "sIRS", "rHOU", "sHOU", "rSPR", "sSPR",
  "rPPM", "sPPM", "rECS", "sECS", "rECT", "sECT", "dIRS", "dHOU", "aOBT", "eSRE",
  "fBIO", "fLAR", "fATSB", "fSSP", "fOVI", "eOVI"
)
TRANSMISSION_PARAMETER_NAMES <- c("f0", "epsilon0", "recRate", "NH_eq", "bV", "bh")

# Picks the named entries of `parameters`, in the given order (errors if a name is missing)
pickParameters = function(parameters, parameterNames){
  vapply(parameterNames, function(name) parameters[[name]], numeric(1))
}

getTheta = function(
  #. getTheta: Facade function to return theta from the selected parameters

  # Move this to driver
  # Species: getAnGambiaeParameters(), getAnArabiensisParameters() or getAnFunestusParameters()
  speciesSpecificParameters=getAnGambiaeParameters(),
  interventionParameters=getInterventionsParameters(),
  additionalTransmissionParameters=getAdditionalTransmissionParameters(),
  mosquitoLifeCycleParameters=getMosquitoLifeCycleParameters()
){
  c(
    pickParameters(mosquitoLifeCycleParameters,      LIFECYCLE_PARAMETER_NAMES),
    pickParameters(interventionParameters,           INTERVENTION_COVERAGE_NAMES),
    pickParameters(interventionParameters,           INTERVENTION_TIME_NAMES),
    pickParameters(speciesSpecificParameters,        SPECIES_PARAMETER_NAMES),
    pickParameters(additionalTransmissionParameters, TRANSMISSION_PARAMETER_NAMES)
  )
}
######################################################################################
#------------------------------------------------------------------------
############ SPECIES SPECIFIC PARAMETERS ################################
#------------------------------------------------------------------------
#* muV:     Adult mosquito daily mortality
#* Q0:      Human blood index
#* phiB:    Proportion of bites on a person while they are in bed
#* phiI:    Proportion of bites on a person while they are indoors
#* rSPR:    Probability of mosquito repeating a feeding attempt due to spatial repellent
#* sSPR:    Probability of mosquito feeding successfully in presence of spatial repellent
#* rPPM:    Probability of mosquito repeating a feeding attempt due to personal protection measures
#* sPPM:    Probability of mosquito feeding successfully in presence of personal protection measures
#* rIRS:    Probability of mosquito repeating a feeding attempt due to IRS
#* sIRS:    Probability of mosquito feeding and surviving in presence of IRS
#* rITN:    Probability of mosquito repeating a feeding attempt due to ITN
#* sITN:    Probability of mosquito feeding and surviving in presence of ITN
#* rECT:    Probability of mosquito repeating a feeding attempt due to insecticide (topical) treated cattle
#* sECT:    Probability of mosquito feeding and surviving in presence of insecticide (topical) treated cattle
#* rECS:    Probability of mosquito repeating a feeding attempt due to insecticide (systemic) treated cattle
#* sECS:    Probability of mosquito feeding and surviving in presence of insecticide (systemic) treated cattle
#* dHOU:    Death due to encountering mosquito proofed housing
#* dIRS:    Death due to encountering IRS treated house
#* aOBT:    Availability of one odor baited trap in relation to one human
#* eSRE:    Effectiveness of the impact of source reduction on K
#* eOVI:    Effectiveness of the ovitrap
#* fATSB    Factor allowing for increased mosquito death due to ATSB
#* fLAR     Factor allowing for increased mosquito death due to larviciding
#* fBIO     Factor allowing for increased mosquito death due to biological control
#* fSSP     Factor allowing for increased mosquito death due to space spraying
#* fOVI     Factor allowing for increased mosquito death due to ovitraps
#------------------------------------------------------------------------
getAnGambiaeParameters = function(){
  #. getAnGambiaeParameters: Returns the list of Anopheles gambiae specific parameters
  c(
    muV = 1/7.6, Q0 = 0.92,   phiB = 0.89, phiI = 0.97,
    rOVI = 0.00, sOVI = 0.00, rFOG = 0.00, sFOG = 0.00,
    rLAR = 0.00, sLAR = 0.00, rBIO = 0.00, sBIO = 0.00,
    rSRE = 0.00, sSRE = 0.00, rIRS = 0.60, sIRS = 0.00,
    rITN = 0.56, sITN = 0.03, rECS = 0.00, sECS = 0.10,
    rHOU = 0.56, sHOU = 0.03, rODO = 0.00, sODO = 0.00,
    rSPR = 0.82, sSPR = 0.12, rECT = 0.30, sECT = 0.50,
    rPPM = 0.69, sPPM = 0.31, dHOU = 0.41, dIRS = 0.40,
    aOBT = 2.00, eSRE = 0.50, fATSB =12.2/7.6, fSSP = 1.50,
    fOVI = 1.50, fLAR = 55.5, fBIO = 38.1, eOVI = 0.7
  )
}
getAnArabiensisParameters = function(){
  #. getAnArabiensisParameters: Returns the list of Anopheles arabiensis specific parameters
  c(
    # Previous values:  Phi = 0.96, phiB = 90 (modified -> get the actual value)
    #                   muV = 1/7.6, Q0 = 0.71, phiB = 0.50, phiI = 0.50,
    muV = 1/7.6, Q0 = 0.71,   phiB = 0.90, phiI = 0.96,
    rOVI = 0.00, sOVI = 0.00, rFOG = 0.00, sFOG = 0.00,
    rLAR = 0.00, sLAR = 0.00, rBIO = 0.00, sBIO = 0.00,
    rSRE = 0.00, sSRE = 0.00, rIRS = 0.60, sIRS = 0.00,
    rITN = 0.48, sITN = 0.39, rECS = 0.90, sECS = 0.0,
    rHOU = 0.48, sHOU = 0.39, rODO = 0.00, sODO = 0.00,
    rSPR = 0.02, sSPR = 0.15, rECT = 0.50, sECT = 0.50,
    rPPM = 0.69, sPPM = 0.31, dHOU = 0.13, dIRS = 0.40,
    aOBT = 2.00, eSRE = 0.50, fATSB =12.2/7.6, fSSP = 1.50,
    fOVI = 1.50, fLAR = 55.5, fBIO = 38.1, eOVI = 0.7
  )
}
getAnFunestusParameters = function(){
  #. getAnFunestusParameters: Returns the list of Anopheles funestus specific parameters
  c(
    muV = 1/8.9, Q0 = 0.94,   phiB = 0.90, phiI = 0.98,
    rOVI = 0.00, sOVI = 0.00, rFOG = 0.00, sFOG = 0.00,
    rLAR = 0.00, sLAR = 0.00, rBIO = 0.00, sBIO = 0.00,
    rSRE = 0.00, sSRE = 0.00, rIRS = 0.63, sIRS = 0.00,
    rITN = 0.56, sITN = 0.03, rECS = 0.60, sECS = 0.50,
    rHOU = 0.56, sHOU = 0.03, rODO = 0.00, sODO = 0.00,
    rSPR = 0.82, sSPR = 0.12, rECT = 0.60, sECT = 0.50,
    rPPM = 0.69, sPPM = 0.31, dHOU = 0.41, dIRS = 0.37,
    aOBT = 2.00, eSRE = 0.5,  fATSB = 12.2/8.9, fSSP = 1.50,
    fOVI = 1.50, fLAR = 55.5, fBIO = 38.1, eOVI = 0.7
  )
}
######################################################################################
#------------------------------------------------------------------------
############ MOSQUITO LIFE CYCLE PARAMETERS #############################
#------------------------------------------------------------------------
#* beta:    Number of eggs laid per day per female mosquito
#* muEL:    Early larval instar daily mortality
#* muLL:    Late larval instar daily mortality
#* muPL:    Pupal daily mortality
#* durEL:   Duration of early instar stage
#* durLL:   Duration of late instar stage
#* durPL:   Duration of pupal stage
#* durEV:   Duration of latent period in mosquito (days)
#* gamma:   Effect of density-dependence on late instars relative to early instars
#* tau1:    Time spent foraging for a blood meal at 0% ITN coverage
#* tau2:    Time spent resting and ovipositing by a mosquito
#------------------------------------------------------------------------
getMosquitoLifeCycleParameters = function(
  beta=21.19,muEL=0.034,muLL=0.035,muPL=0.25,
  durEL=6.64,durLL=6.64,durPL=0.64,durEV=10,gamma=13.25,
  tau1=0.68,tau2=2.32
){
  #. getMosquitoLifeCycleParameters: Returns the list of mosquito life cycle parameters
  c(
    beta=beta,
    muEL=muEL,muLL=muLL,muPL=muPL,
    durEL=durEL,durLL=durLL,durPL=durPL,durEV=durEV,
    gamma=gamma,tau1=tau1,tau2=tau2
  )
}
######################################################################################
