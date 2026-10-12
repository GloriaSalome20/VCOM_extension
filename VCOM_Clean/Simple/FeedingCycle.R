########################################################################
# FeedingCycle.R
# Computes, in the presence of interventions:
#        1) Probability of surviving the feeding period
#        2) Probability of surviving resting and oviposition
#        3) Probability of surviving one day
#        4) Mortality rate of female mosquitoes (muVCom)
#        5) Eggs laid per day per female mosquito (betaCom)
#        6) Human biting rate (a_theta)
# (This single file replaces the identical FeedingCycle.R / ODEFeedingCycle.R pair.)
########################################################################

# Could just pass theta instead of the long list of arguments.
# e_ov and c0 are re-computed inside the function; they are only kept in the argument list
# so that existing calls keep working.
impactFeedingCycleParameters = function(time,beta,tau1,tau2,e_ov=NULL,time_ATSB_on,ATSBcov,time_SSP_on,SSPcov,fSSP,fATSB,muV,
                                        Q0,aOBT,OBTcov,time_OBT_on,time_ITN_on,ITNcov,time_IRS_on,IRScov,HOUcov,
                                        time_HOU_on,rITN,sITN,rIRS,rHOU,sIRS,sHOU, phiB, phiI,dHOU,dIRS,
                                        time_SPR_on,SPRcov,time_PPM_on,PPMcov,rSPR,rPPM,sSPR,sPPM,
                                        c0=NULL,time_ECS_on,ECScov,time_ECT_on,ECTcov,rECT,sECS,sECT,
                                        time_OVI_on,OVIcov,fOVI,eOVI){
  #. impactFeedingCycleParameters: compute the impact of interventions on the feeding cycle parameters
  #  Returns c(muVCom, betaCom, a_theta)
  #  The intermediate results are also stored in the global environment (<<-) for inspection/debugging.

  delta <- 1/(tau1+tau2)                  # Inverse of gonotrophic cycle without interventions
  e_ov <- beta*(exp(muV/delta)-1)/muV     # Number of eggs per oviposition per mosquito

  ##**********Host seeking - ATSB and space spraying ****##############################
  muV_1_Com <<- impactATSBSpaceSpraying(time,time_ATSB_on,ATSBcov,time_SSP_on,SSPcov,fSSP,fATSB,muV)

  ##************************Host seeking - Odor baited traps********************
  impactOdorT <<- impactOdorBaitedTraps(time,Q0,aOBT,OBTcov,time_OBT_on)

  Q0_t_h = impactOdorT[1]
  Q0_t_c = impactOdorT[2]   # Proportion to cattle affected by traps

  ##*************************Human - Indoor protection *********************************##
  impactIndoor <<- impactIndoorProtection(time,time_ITN_on,ITNcov,time_IRS_on,IRScov,HOUcov,time_HOU_on,
                                          rITN,sITN,rIRS,rHOU,sIRS,sHOU, Q0_t_h,Q0, phiB, phiI,dHOU,dIRS,time_OBT_on)

  zCom_Human_Indoor = impactIndoor[1]
  wCom_Human_Indoor = impactIndoor[2]
  c0                = impactIndoor[3]  # Share of the population without indoor interventions (see eqn)

  ##*************************Human - Outdoor protection *********************************##
  impactOutdoor <<- impactOutdoorProtection(time,time_SPR_on,SPRcov,time_PPM_on,PPMcov,rSPR,rPPM,sSPR,sPPM,Q0_t_h,Q0,phiI,c0,time_OBT_on)

  zCom_Human_Outdoor = impactOutdoor[1]
  wCom_Human_Outdoor = impactOutdoor[2]

  ##*************************** Cattle *********************************************##
  impactCattle = impactInsecticideTreatedCattle(time,time_ECS_on,ECScov,time_ECT_on,ECTcov,rECT,sECS,sECT,Q0_t_c,Q0,time_OBT_on)
  zCom_Cattle = impactCattle[1]
  wCom_Cattle = impactCattle[2]

  ##*************************** Resting & Ovipositing *********************************************##
  muV_2_Com = impactRestingOvipositing(time,time_OVI_on,OVIcov,time_ATSB_on,ATSBcov,SSPcov,time_SSP_on,
                                       fOVI,fATSB,fSSP,muV,eOVI)

  ######**************************Computing overall impact***********************************#####
  ## zCom: Probability of a mosquito being repelled : SAM CHECK THIS
  zCom <- zCom_Cattle + zCom_Human_Outdoor + zCom_Human_Indoor
  # (alternative without outdoor protection: zCom <- zCom_Cattle + zCom_Human_Indoor)

  # deltaCom: Inverse of gonotrophic cycle length with ITNs & IRS
  deltaCom <- 1/(tau1/(1-zCom) + tau2)

  # Probability that a surviving mosquito succeeds in feeding during a single attempt - only on humans
  wCom_human = wCom_Human_Outdoor + wCom_Human_Indoor

  ### wCom: Probability that a surviving mosquito succeeds in feeding during a single attempt
  wCom <- wCom_Cattle + wCom_Human_Outdoor + wCom_Human_Indoor
  # (alternative without outdoor protection: wCom <- wCom_Cattle + wCom_Human_Indoor)

  # Probability of surviving the feeding period in the presence of an intervention
  p10 <- exp(-muV_1_Com*tau1)

  # Probability of surviving the first attempt, second, and so on - eqn 3 (SK)
  p1Com <- p10*wCom/(1 - zCom*p10)

  # Probability of resting (SK)
  p2 <- exp(-muV_2_Com*tau2)

  # Probability of surviving one day
  pCom <- (p1Com*p2)^deltaCom
  # Female mosquito mortality rate due to ITN and IRS (SK)
  muVCom <- -log(pCom)
  # betaCom: Eggs laid per day by female mosquitoes in presence of ITNs & IRS
  betaCom <- e_ov*muVCom/(exp(muVCom/deltaCom) - 1)

  # Length of the gonotrophic cycle
  f_theta = computeLengthGonotrophicycle(deltaCom)

  # Human blood index
  HBI_com = computeHBI(wCom_human,wCom)

  # Human biting rate
  a_theta =  computeHumanBitingRate(f_theta,HBI_com)

  impactFeedingCycle <- c(muVCom,betaCom,a_theta)

  return(impactFeedingCycle)
}
