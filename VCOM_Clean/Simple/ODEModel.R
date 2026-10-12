########################################################################
# ODEModel.R
# The main ODE model of the mosquito population (with a human SIR component).
#
#   runODE()                  wrapper that integrates the model with deSolve::dede()
#   calculateInitialState()   equilibrium initial conditions
#   IVM_ode()                 model definition (right-hand side of the ODE system)
########################################################################

# Position of SV (susceptible vectors) in the state vector; used to read the lagged state.
# State order: EL, EL_LAR, EL_BIO, EL_LAR_BIO, LL, LL_LAR, LL_BIO, LL_LAR_BIO, PL, SV, EV, IV, SH, IH, Mm
# If you add or reorder states, make sure this index still points to SV.
SV_INDEX <- 10

# Value used instead of an exact zero for the larval coverage shares (avoids dividing by zero)
ZERO_COVERAGE_GUARD <- 2*exp(-16)

# Larval carrying capacity over time ---------------------------------------------------------------
# IVM_ode() computes the baseline (equilibrium) carrying capacity K0 and then asks
#   carryingCapacity(time, K0)
# for the carrying capacity K at that time. The default below keeps K constant (K = K0), which is
# what the model did before and what main.R and the Shiny app still use.
# Driver.R can replace it with the seasonal carrying capacity made by KL_gamma_temp.R
# (rainfall memory x temperature), using makeSeasonalCarryingCapacity().
if (!exists("carryingCapacity")) carryingCapacity <- function(time, K0) K0

makeDailyCycle = function(daily_values){
  #. makeDailyCycle: builds a function of the model time from a table of daily values
  #.   (the first value is model day 0, linear between days, repeating every length(daily_values) days)
  n        <- length(daily_values)
  interp   <- approxfun(x = 0:n, y = c(daily_values, daily_values[1]))   # day n joins back to day 0
  function(time) interp(time %% n)
}

makeSeasonalCarryingCapacity = function(K_daily, scale = c("relative", "absolute")){
  #. makeSeasonalCarryingCapacity: builds carryingCapacity(time, K0) from a table of daily values K_daily
  #.   (the first value is model day 0, and the table repeats every length(K_daily) days).
  #.   scale = "relative": K(time) = K0 * K_daily(time)/max(K_daily), i.e. the peak of K(t) is the model's K0
  #.   scale = "absolute": K(time) = K_daily(time), the values exactly as in the table
  scale   <- match.arg(scale)
  K_cycle <- makeDailyCycle(K_daily)
  if (scale == "relative") function(time, K0) K0*K_cycle(time)/max(K_daily)
  else                     function(time, K0) K_cycle(time)
}

# Mosquito life-cycle parameters over time -----------------------------------------------------------
# IVM_ode() asks timeVaryingParameters(time, theta) for the parameters beta, muEL, muLL, muPL, muV,
# durEL, durLL, durPL, durEV, tau1 and tau2 at the current time. The default returns theta unchanged,
# i.e. constant parameters, which is what the model did before and what main.R and the Shiny app still use.
# Driver.R can replace it with the temperature- and humidity-dependent parameters of ODETemperatureParameters.R.
# The equilibrium (calculateInitialState(), NV_eq, lambdaV) and the baseline carrying capacity K0 always use
# the theta that is passed to the solver, so they do not change with time.
if (!exists("timeVaryingParameters")) timeVaryingParameters <- function(time, theta) theta

######################################################################################
runODE = function(tMax, tIncrement, initState, theta, method){
  #. runODE: Main ODE wrapper for simulating the mosquito population (runs from day 0 to tMax)
  times <- seq(0, tMax, by = tIncrement)
  IVM_traj <- data.frame(dede(y = initState, times = times, parms = theta, func = IVM_ode, method = method))
  return(IVM_traj)
}
######################################################################################
computeOmega = function(beta, muEL, muLL, muPL, muV, durEL, durLL, durPL, gamma){
  #. computeOmega: ratio used for the equilibrium larval densities and the carrying capacity K
  b_omega <- gamma*muLL/muEL - durEL/durLL + (gamma-1)*muLL*durEL
  omega <- -0.5*b_omega + sqrt(0.25*b_omega^2 + gamma*beta*muLL*durEL/(2*muEL*muV*durLL*(1+durPL*muPL)))
  return(omega)
}
######################################################################################
calculateInitialState = function(theta){
  #. calculateInitialState: Calculates the initial conditions of a system given the theta parameters
  #  SIDE EFFECT: adds lambdaV and NV_eq to the GLOBAL theta (via <<-); IVM_ode() reads them from there.
  beta <- theta[["beta"]]
  muEL <- theta[["muEL"]]
  muLL <- theta[["muLL"]]
  muPL <- theta[["muPL"]]
  muV <- theta[["muV"]]
  durEL <- theta[["durEL"]]
  durLL <- theta[["durLL"]]
  durPL <- theta[["durPL"]]
  durEV <- theta[["durEV"]]
  gamma <- theta[["gamma"]]
  Q0 <- theta[["Q0"]]
  f0 <- theta[["f0"]]
  epsilon0 <- theta[["epsilon0"]]
  recRate <- theta[["recRate"]]
  NH_eq <- theta[["NH_eq"]]
  bV <- theta[["bV"]]
  bh <- theta[["bh"]]

  omega <- computeOmega(beta, muEL, muLL, muPL, muV, durEL, durLL, durPL, gamma)
  a0 <- Q0*f0 # Human biting rate at equilibrium

  # Humans at equilibrium
  iH_eq <- (epsilon0 *bh)/(recRate+(epsilon0*bh))
  IH_eq <- iH_eq*NH_eq
  SH_eq = NH_eq* (1-iH_eq)

  lambdaV <- a0*iH_eq*bV # Force of infection in mosquitoes at equilibrium
  theta["lambdaV"] <<- lambdaV # Include vector force of infection in vector of parameters (theta)

  # Mosquitoes at equilibrium
  iV_eq <- lambdaV*exp(-muV*durEV)/(lambdaV + muV)
  sV_eq <- iV_eq*muV/(lambdaV*exp(-muV*durEV))
  eV_eq <- 1 - sV_eq - iV_eq

  NV_eq <- epsilon0*NH_eq/(iV_eq*a0)
  theta["NV_eq"] <<- NV_eq # Include equilibrium vector population size in vector of parameters (theta)

  EL_eq <- 2*omega*muV*durLL*(1 + muPL*durPL)*NV_eq
  LL_eq <- 2*muV*durLL*(1 + muPL*durPL)*NV_eq
  PL_eq <- 2*muV*durPL*NV_eq
  SV_eq <- sV_eq*NV_eq
  EV_eq <- eV_eq*NV_eq
  IV_eq <- iV_eq*NV_eq

  # The compartments with larvicide / biological control (EL_LAR, EL_BIO, ...) start at 0.
  # Males (Mm) start at SV_eq.
  initState <- c(EL = EL_eq,EL_LAR = 0,EL_BIO = 0,EL_LAR_BIO = 0,
                 LL = LL_eq,LL_LAR = 0,LL_BIO = 0,LL_LAR_BIO = 0,
                 PL = PL_eq,SV = SV_eq,EV = EV_eq,IV = IV_eq,
                 SH = SH_eq, IH = IH_eq, Mm = SV_eq)

  return(initState)
}
######################################################################################
IVM_ode <- function(time, state, theta){
  #. IVM_ode: ODE Model definition

  ## ---- Parameters (mosquito life cycle) ------------------------------------------
  beta <- theta[["beta"]]   # Eggs laid per day by female mosquito
  muEL <- theta[["muEL"]]   # Early instar stage daily mortality
  muLL <- theta[["muLL"]]   # Late instar stage daily mortality
  muPL <- theta[["muPL"]]   # Pupal stage daily mortality
  muV <- theta[["muV"]]     # Adult mosquito daily mortality
  Q0 <- theta[["Q0"]]       # Human blood index
  phiB <- theta[["phiB"]]   # Proportion of bites on a person while they are in bed
  phiI <- theta[["phiI"]]   # Proportion of bites on a person while they are indoors
  durEL <- theta[["durEL"]] # Duration of early instar stage (days)
  durLL <- theta[["durLL"]] # Duration of late instar stage (days)
  durPL <- theta[["durPL"]] # Duration of pupal stage (days)
  durEV <- theta[["durEV"]] # Duration of latent period in mosquito (days)
  gamma <- theta[["gamma"]] # Effect of density-dependence on late instars relative to early instars
  tau1 <- theta[["tau1"]]   # Time spent foraging for a blood meal (no ITNs) (days)
  tau2 <- theta[["tau2"]]   # Time spent resting and ovipositing (days)
  NV_eq <- theta[["NV_eq"]] # Number of female mosquitoes at equilibrium
  recRate <- theta[["recRate"]]
  lambdaV <- theta[["lambdaV"]] # Force of infection in vectors at equilibrium

  ## ---- Parameters (interventions) -----------------------------------------------------------
  # LLINs and IRS
  ITNcov <- theta[["ITNcov"]]               # ITN coverage
  IRScov <- theta[["IRScov"]]               # IRS coverage
  time_ITN_on <- theta[["time_ITN_on"]]     # When ITNs are applied (days)
  time_IRS_on <- theta[["time_IRS_on"]]     # When IRS is applied (days)
  rITN <- theta[["rITN"]]                   # Probability of mosquito repeating a feeding attempt due to ITN
  sITN <- theta[["sITN"]]                   # Probability of mosquito feeding and surviving in presence of ITNs
  rIRS <- theta[["rIRS"]]                   # Probability of mosquito repeating a feeding attempt due to IRS
  sIRS <- theta[["sIRS"]]                   # Probability of mosquito feeding and surviving in presence of IRS
  dIRS <- theta[["dIRS"]]                   # Death rate after encountering IRS
  dHOU <- theta[["dHOU"]]                   # Death rate after encountering mosquito proofed housing

  # Mosquito proofed housing
  HOUcov <- theta[["HOUcov"]]               # Proportion of houses that are mosquito proofed
  time_HOU_on <- theta[["time_HOU_on"]]     # When house modification is applied (days)
  rHOU <- theta[["rHOU"]]                   # Probability of mosquito repeating a feeding attempt due to mosquito proofed housing
  sHOU <- theta[["sHOU"]]                   # Probability of mosquito feeding and surviving in presence of mosquito proofed housing

  # Endectocide-treated cattle (either topical or systemic)
  ECScov <- theta[["ECScov"]]               # Endectocide (e.g., ivermectin) coverage - systemic
  time_ECS_on <- theta[["time_ECS_on"]]     # When cattle are treated (systemic) (days)
  sECS <- theta[["sECS"]]                   # Probability of mosquito feeding and surviving in presence of systemic-treated cattle
  ECTcov <- theta[["ECTcov"]]               # Endectocide (e.g., ivermectin) coverage - topical
  time_ECT_on <- theta[["time_ECT_on"]]     # When cattle are treated (topical) (days)
  sECT <- theta[["sECT"]]                   # Probability of mosquito feeding and surviving in presence of topical-treated cattle
  rECT <- theta[["rECT"]]                   # Probability of mosquito repeating a feeding attempt due to topical-treated cattle

  # Protecting humans outdoors - spatial repellents (SPR) and personal protection measures (PPM)
  SPRcov <- theta[["SPRcov"]]               # Proportion of a defined space with spatial repellents
  time_SPR_on <- theta[["time_SPR_on"]]     # When SPR is on
  rSPR <- theta[["rSPR"]]                   # Probability of mosquito repeating a feeding attempt due to SPR
  sSPR <- theta[["sSPR"]]                   # Probability of mosquito feeding in presence of SPR
  PPMcov <- theta[["PPMcov"]]               # Proportion of people outside with PPM
  time_PPM_on <- theta[["time_PPM_on"]]     # When PPM is on
  rPPM <- theta[["rPPM"]]                   # Probability of mosquito repeating a feeding attempt due to PPM
  sPPM <- theta[["sPPM"]]                   # Probability of mosquito feeding in presence of PPM

  # Odor baited traps
  time_OBT_on <- theta[["time_OBT_on"]]     # When OBT is on
  aOBT <- theta[["aOBT"]]                   # Availability of one trap in relation to one human
  OBTcov <- theta[["OBTcov"]]               # Ratio of traps to humans (i.e., coverage)

  # Source reduction
  time_SRE_on <- theta[["time_SRE_on"]]     # When source reduction is on
  eSRE <- theta[["eSRE"]]                   # Effectiveness of source reduction
  SREcov <- theta[["SREcov"]]               # Proportion of aquatic habitats covered by source reduction

  # Host seeking - ATSB and space spraying
  time_ATSB_on <- theta[["time_ATSB_on"]]   # When ATSB is on
  fATSB <- theta[["fATSB"]]                 # Factor allowing for increased death rate due to ATSB
  ATSBcov <- theta[["ATSBcov"]]             # Proportion of vegetation sprayed with ATSB
  time_SSP_on <- theta[["time_SSP_on"]]     # When space spraying is on
  fSSP <- theta[["fSSP"]]                   # Factor allowing for increased death rate due to space spraying
  SSPcov <- theta[["SSPcov"]]               # Proportion of mosquito habitats in a defined area sprayed with insecticide

  # Resting and ovipositing - ATSB, space spraying, ovitraps (ATSB and SSP are defined above)
  time_OVI_on <- theta[["time_OVI_on"]]     # When ovitraps are on
  fOVI <- theta[["fOVI"]]                   # Factor allowing for increased death rate due to ovitraps
  eOVI <- theta[["eOVI"]]                   # Effectiveness of the ovitrap
  OVIcov <- theta[["OVIcov"]]               # Proportion of ovitraps

  # Larviciding and biological control
  # NOTE (unchanged from the original): fBIO and fLAR are both read from theta[["fOVI"]],
  # not from theta[["fBIO"]] / theta[["fLAR"]].
  time_BIO_on <- theta[["time_BIO_on"]]     # When biological control is on
  fBIO <- theta[["fOVI"]]                   # Factor allowing for increased death rate due to biological control
  BIOcov <- theta[["BIOcov"]]               # Proportion of aquatic habitats covered by biological control
  time_LAR_on <- theta[["time_LAR_on"]]     # When larvicide is on
  fLAR <- theta[["fOVI"]]                   # Factor allowing for increased death rate due to larvicide
  LARcov <- theta[["LARcov"]]               # Proportion of aquatic habitats covered by larvicide

  ## ---- States (definitions added by SK) ---------------------------------------------------
  EL         <- state[["EL"]]          # Early instar stage
  EL_LAR     <- state[["EL_LAR"]]      # Early instar stage in presence of larviciding
  EL_BIO     <- state[["EL_BIO"]]      # Early instar stage in presence of biological control
  EL_LAR_BIO <- state[["EL_LAR_BIO"]]  # Early instar stage in presence of both
  LL         <- state[["LL"]]          # Late instar stage
  LL_LAR     <- state[["LL_LAR"]]      # Late instar stage in presence of larviciding
  LL_BIO     <- state[["LL_BIO"]]      # Late instar stage in presence of biological control
  LL_LAR_BIO <- state[["LL_LAR_BIO"]]  # Late instar stage in presence of both
  PL         <- state[["PL"]]          # Pupal stage
  SV         <- state[["SV"]]          # Susceptible vectors
  EV         <- state[["EV"]]          # Latent period
  IV         <- state[["IV"]]          # Infectious vectors
  Mm         <- state[["Mm"]]          # Adult male mosquitoes
  SH         <- state[["SH"]]          # Susceptible humans (SIR)
  IH         <- state[["IH"]]          # Infectious humans (SIR)

  ## ---- Derived parameters (SK - IMPORTANT: add to the write up) ------------------------------
  NV <- SV + EV + IV # Total mosquito population size
  omega <- computeOmega(beta, muEL, muLL, muPL, muV, durEL, durLL, durPL, gamma)
  K0 <- 2*NV_eq*muV*durLL*(1 + muPL*durPL)*gamma*(omega+1)/(omega/(muLL*durEL) - 1/(muLL*durLL) - 1) # Baseline (equilibrium) larval carrying capacity
  K  <- carryingCapacity(time, K0)  # Larval carrying capacity at this time (= K0 unless Driver.R switches on the seasonal K)

  ## ---- Life-cycle parameters at this time (constants unless Driver.R switches on the temperature-dependent ones) --
  # K0 above is computed from the baseline theta on purpose: the equilibrium does not move with the temperature.
  theta_t <- timeVaryingParameters(time, theta)
  beta  <- theta_t[["beta"]]
  muEL  <- theta_t[["muEL"]]
  muLL  <- theta_t[["muLL"]]
  muPL  <- theta_t[["muPL"]]
  muV   <- theta_t[["muV"]]
  durEL <- theta_t[["durEL"]]
  durLL <- theta_t[["durLL"]]
  durPL <- theta_t[["durPL"]]
  durEV <- theta_t[["durEV"]]
  tau1  <- theta_t[["tau1"]]
  tau2  <- theta_t[["tau2"]]

  ## ---- Aquatic habitats - impact of source management on K -----------------------------------
  K_sr    <<- impactSourceReduction(time,eSRE,SREcov,time_SRE_on,K)
  K = K_sr   # Update K

  ## ---- Feeding cycle in presence of interventions ----------------------------------------------
  # (the equations of the gonotrophic cycle live in FeedingCycle.R)
  feedingCycleImpact = impactFeedingCycleParameters(
    time = time, beta = beta, tau1 = tau1, tau2 = tau2,
    time_ATSB_on = time_ATSB_on, ATSBcov = ATSBcov, time_SSP_on = time_SSP_on, SSPcov = SSPcov,
    fSSP = fSSP, fATSB = fATSB, muV = muV,
    Q0 = Q0, aOBT = aOBT, OBTcov = OBTcov, time_OBT_on = time_OBT_on,
    time_ITN_on = time_ITN_on, ITNcov = ITNcov, time_IRS_on = time_IRS_on, IRScov = IRScov,
    HOUcov = HOUcov, time_HOU_on = time_HOU_on,
    rITN = rITN, sITN = sITN, rIRS = rIRS, rHOU = rHOU, sIRS = sIRS, sHOU = sHOU,
    phiB = phiB, phiI = phiI, dHOU = dHOU, dIRS = dIRS,
    time_SPR_on = time_SPR_on, SPRcov = SPRcov, time_PPM_on = time_PPM_on, PPMcov = PPMcov,
    rSPR = rSPR, rPPM = rPPM, sSPR = sSPR, sPPM = sPPM,
    time_ECS_on = time_ECS_on, ECScov = ECScov, time_ECT_on = time_ECT_on, ECTcov = ECTcov,
    rECT = rECT, sECS = sECS, sECT = sECT,
    time_OVI_on = time_OVI_on, OVIcov = OVIcov, fOVI = fOVI, eOVI = eOVI)

  muVCom  = feedingCycleImpact[1]
  betaCom = feedingCycleImpact[2]
  a_theta = feedingCycleImpact[3]

  # Impact of ovitraps on egg laying
  OVIcov_t <- activeCoverage(time, OVIcov, time_OVI_on)
  betaCom = betaCom *(1-eOVI*OVIcov_t)

  ## ---- Other outputs: EIR, R0 and vectorial capacity ---------------------------------------------
  # NOTE (unchanged): the transmission values come from getAdditionalTransmissionParameters()
  # with its default arguments, not from theta.
  transmissionValues = getAdditionalTransmissionParameters()
  bv = transmissionValues[["bV"]]
  NH = transmissionValues[["NH_eq"]]
  bh = transmissionValues[["bh"]]

  EIR <- computeEIR(a_theta, IV, NH)              # Entomological inoculation rate
  R0 <- computeRO(a_theta,muVCom, NV,bv,bh,NH, durEV)   # Basic reproductive number
  VC <- computeVC(a_theta, NV,NH,muVCom, durEV)          # Vectorial capacity

  ## ---- Larviciding and biological control: share of habitats under each treatment ----------------
  larBio <- impactLarvacidingBiologicalControl(time,time_LAR_on, LARcov,time_BIO_on,BIOcov,fLAR,fBIO)
  cLAR       <- larBio$cLAR         # Coverage with larviciding only
  cBIO       <- larBio$cBIO         # Coverage with biological control only
  c_LAR_BIO  <- larBio$c_LAR_BIO    # Both applied
  c0_LAR_BIO <- larBio$c0_LAR_BIO   # Neither applied
  f_LAR_BIO  <- larBio$f_LAR_BIO

  # Avoid a share being exactly 0
  if (cLAR == 0 )       cLAR       = ZERO_COVERAGE_GUARD
  if (cBIO == 0 )       cBIO       = ZERO_COVERAGE_GUARD
  if (c_LAR_BIO == 0 )  c_LAR_BIO  = ZERO_COVERAGE_GUARD
  if (c0_LAR_BIO == 0 ) c0_LAR_BIO = ZERO_COVERAGE_GUARD

  ## ---- Lagged susceptible vectors (latent period) -------------------------------------------------
  if(time < durEV){
    SVLag <- SV
  }else{
    lagStates <- lagvalue(time-durEV)
    SVLag <- lagStates[SV_INDEX]
  }

  ## ---- ODEs: aquatic stages -----------------------------------------------------------------------
  # Habitats with no larvicide / biological control
  dEL <- c0_LAR_BIO*betaCom*NV - muEL*(1 + ( (EL+LL)/(K*c0_LAR_BIO) ))*EL - EL/durEL
  # Removed c0_LAR_BIO per suggestions
  dLL <- EL/durEL - muLL*(1 + gamma*((EL+LL)/(K*c0_LAR_BIO)))*LL - LL/durLL
  dPL <- (LL+LL_LAR+LL_BIO+LL_LAR_BIO)/durLL - muPL*PL - PL/durPL

  # Habitats with larviciding and/or biological control
  larviciding_on <- (time > time_LAR_on) && (cLAR > 0)
  biocontrol_on  <- (time > time_BIO_on) && (cBIO > 0)
  both_on        <- larviciding_on && biocontrol_on

  if (larviciding_on) {
    dEL_LAR <- cLAR*betaCom*NV - fLAR*muEL*(1 + ((EL_LAR+LL_LAR)/(cLAR*K)))*EL_LAR - EL_LAR/durEL
  } else {
    dEL_LAR <- 0
  }
  if (biocontrol_on) {
    dEL_BIO <- cBIO*betaCom*NV - fBIO*muEL*(1 + ((EL_BIO+LL_BIO)/(cBIO*K)))*EL_BIO - EL_BIO/durEL
  } else {
    dEL_BIO <- 0
  }
  if (both_on) {
    dEL_LAR_BIO <- c_LAR_BIO*betaCom*NV - f_LAR_BIO*muEL*(1 + ((EL_LAR_BIO+LL_LAR_BIO)/(c_LAR_BIO*K)))*EL_LAR_BIO - EL_LAR_BIO/durEL
  } else {
    dEL_LAR_BIO <- 0
  }

  # NOTE (unchanged): the inflow to every treated late-instar compartment is EL/durEL
  # (from the untreated early instars), as in the original code.
  if (larviciding_on) {
    dLL_LAR <- EL/durEL - fLAR*muLL*(1 + gamma*((EL_LAR+LL_LAR)/(K*cLAR)))*LL_LAR - LL_LAR/durLL
  } else {
    dLL_LAR <- 0
  }
  if (biocontrol_on) {
    dLL_BIO <- EL/durEL - fBIO*muLL*(1 + gamma*((EL_BIO+LL_BIO)/(K*cBIO)))*LL_BIO - LL_BIO/durLL
  } else {
    dLL_BIO <- 0
  }
  if (both_on) {
    dLL_LAR_BIO <- EL/durEL - f_LAR_BIO*muLL*(1 + gamma*((EL_LAR_BIO+LL_LAR_BIO)/(K*c_LAR_BIO)))*LL_LAR_BIO - LL_LAR_BIO/durLL
  } else {
    dLL_LAR_BIO <- 0
  }

  ## ---- ODEs: adult females (0.5 -> only female mosquitoes) ---------------------------------------
  dSV <- 0.5*PL/durPL - lambdaV*SV - muVCom*SV
  dEV <- lambdaV*SV - lambdaV*SVLag*exp(-muVCom*durEV) - muVCom*EV
  dIV <- lambdaV*SVLag*exp(-muVCom*durEV) - muVCom*IV

  ## ---- ODEs: adult males ------------------------------------------------------------------------
  muMm = muV
  f_SwS = 11
  muMmCom = muMm*f_SwS
  dMm <- 0.5*PL/durPL - muMmCom*Mm

  ## ---- ODEs: SIR model for humans ---------------------------------------------------------------
  dSH <- -(NV/NH)*a_theta*(IV/NV)*bh*SH +  recRate*IH
  dIH <-  (NV/NH)*a_theta*(IV/NV)*bh*SH -  recRate*IH

  # If you add more states, make sure SV_INDEX (top of this file) still points to SV
  return(list(c(dEL,dEL_LAR,dEL_BIO,dEL_LAR_BIO, dLL,dLL_LAR,dLL_BIO,dLL_LAR_BIO, dPL, dSV, dEV, dIV,dSH,dIH,dMm), EIR=EIR, VC = VC, R0 = R0))
}
######################################################################################
