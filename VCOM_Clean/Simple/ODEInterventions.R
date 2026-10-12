########################################################################
# ODEInterventions.R
# Impact of every intervention on the mosquito life cycle.
#
#   Aquatic habitats   : source reduction (SRE), larviciding (LAR), biological control (BIO)
#   Host seeking       : ATSB, space spraying (SSP), odor-baited traps (OBT)
#   Humans, indoors    : LLINs (ITN), IRS, mosquito-proofed housing (HOU)
#   Humans, outdoors   : spatial repellents (SPR), personal protection measures (PPM)
#   Cattle             : insecticide-treated cattle (ECS systemic, ECT topical)
#   Resting/ovipositing: ovitraps (OVI), ATSB, space spraying
#
# Notation used below
#   XXXcov_t  coverage of intervention XXX at the current time (0 before time_XXX_on)
#   cXXX      share of the population covered by XXX only
#   cA_B      share covered by both A and B
#   c0        share covered by none of them
#   zCom      probability that a mosquito is repelled and repeats its feeding attempt
#   wCom      probability that a surviving mosquito succeeds in feeding in one attempt
########################################################################


##*********************** Shared helpers ****************************##

# Coverage that is in force at `time`: zero until the intervention is switched on
activeCoverage = function(time, coverage, time_on){
  if (time > time_on) coverage else 0
}

# Splits the population between two overlapping interventions A and B
splitTwoCoverages = function(a, b){
  list(
    onlyA   = a - a*b,           # covered by A only
    onlyB   = b - a*b,           # covered by B only
    both    = a*b,               # covered by both
    neither = 1 - a - b + a*b    # covered by neither
  )
}


##**********Aquatic habitats - impact of source management ****##############################
impactSourceReduction = function(time,eSRE,SREcov,time_SRE_on,K){
  #. impactSourceReduction: compute updated K due to source management
  SREcov_t <- activeCoverage(time, SREcov, time_SRE_on)

  K_sr = (1-eSRE*SREcov_t)*K
  return(K_sr)
}


##**********Aquatic habitats - impact of larviciding and biological control ****##############################
impactLarvacidingBiologicalControl = function(time,time_LAR_on, LARcov,time_BIO_on,BIOcov,fLAR,fBIO){
  #. impactLarvacidingBiologicalControl: compute the share of aquatic habitats under larviciding only (cLAR),
  #  biological control only (cBIO), both (c_LAR_BIO) or neither (c0_LAR_BIO), and the combined factor f_LAR_BIO
  LARcov_t <- activeCoverage(time, LARcov, time_LAR_on)
  BIOcov_t <- activeCoverage(time, BIOcov, time_BIO_on)

  cov <- splitTwoCoverages(LARcov_t, BIOcov_t)

  list(
    cLAR       = cov$onlyA,
    cBIO       = cov$onlyB,
    c_LAR_BIO  = cov$both,
    c0_LAR_BIO = cov$neither,
    f_LAR_BIO  = fLAR*fBIO
  )
}


##**********Host Seeking - ATSB and Space Spraying ****##############################
impactATSBSpaceSpraying = function(time,time_ATSB_on,ATSBcov,time_SSP_on,SSPcov,fSSP,fATSB,muV){
  #. impactATSBSpaceSpraying: compute the impact of ATSB and SS on P1 (host-seeking mortality)

  # NOTE (unchanged from the original): both branches of the original `if (time > time_ATSB_on)`
  # returned 0, so ATSB has no effect on host-seeking mortality here (it acts in
  # impactRestingOvipositing() instead). To switch it on use:
  #   ATSBcov_t <- activeCoverage(time, ATSBcov, time_ATSB_on)
  ATSBcov_t <- 0
  SSPcov_t  <- activeCoverage(time, SSPcov, time_SSP_on)

  cov <- splitTwoCoverages(ATSBcov_t, SSPcov_t)   # onlyA = ATSB, onlyB = SSP

  # Impact on mosquito mortality
  muV_1 = muV*(cov$onlyA*fATSB + cov$onlyB*fSSP + cov$both*fATSB*fSSP + cov$neither)

  return(muV_1)
}


##**********The impact of Odor baited traps ****##############################
impactOdorBaitedTraps = function(time,Q0,aOBT,OBTcov,time_OBT_on){
  #. impactOdorBaitedTraps: compute the impact of odor baited traps in reducing Q0
  #  aOBT   = availability of one trap in relation to one human
  #  OBTcov = coverage, i.e. ratio of traps to humans
  #  (the switch-on time is applied by the callers, which use these values only after time_OBT_on)

  # Old formula:
  #   Q0_t_h = Q0*(1/(1+aOBT*OBTcov_t))        # impact on Q0 for humans
  #   Q0_t_c = (1-Q0)*(1/(1+aOBT*OBTcov_t))    # impact on 1-Q0, the proportion going to cattle

  # Formula edited with NC comments:
  Nh = theta[["NH_eq"]]   # Number of humans (read from the global theta)
  # Change these later to make them easier to change
  a_c = 1                 # Availability of cattle in relation to humans
  r_c = 0.5               # Ratio of cattle to humans
  Q0_t_h = Nh/(Nh+(a_c*r_c*Nh)+(aOBT*OBTcov*Nh))              # Impact on Q0 for humans
  Q0_t_c = (a_c*r_c*Nh)/(Nh+(a_c*r_c*Nh)+(aOBT*OBTcov*Nh))    # Impact on proportion going to cattle
  impactOdor <- c(Q0_t_h,Q0_t_c)

  return(impactOdor)
}


##***********************************Protecting Humans Indoor**********************##
impactIndoorProtection = function(time,time_ITN_on,ITNcov,time_IRS_on,IRScov,HOUcov,time_HOU_on,
                                  rITN,sITN,rIRS,rHOU,sIRS,sHOU, Q0_t_h,Q0, phiB, phiI,dHOU,dIRS,time_OBT_on){
  #. impactIndoorProtection: compute both zcom and wcom for indoor interventions
  # Incorporate odor baited traps
  if (time > time_OBT_on) Q0 <- Q0_t_h

  ITNcov_t <- activeCoverage(time, ITNcov, time_ITN_on)
  IRScov_t <- activeCoverage(time, IRScov, time_IRS_on)
  HOUcov_t <- activeCoverage(time, HOUcov, time_HOU_on)

  ## ---- Share of the population under each combination of ITN / IRS / HOU
  # Coverage with LLINs only
  cITN <- ITNcov_t - ITNcov_t*IRScov_t - ITNcov_t*HOUcov_t + ITNcov_t*IRScov_t*HOUcov_t
  # Coverage with IRS only
  cIRS <- IRScov_t - ITNcov_t*IRScov_t - IRScov_t*HOUcov_t + ITNcov_t*IRScov_t*HOUcov_t
  # Coverage with house modification only
  cHOU <- HOUcov_t - ITNcov_t*HOUcov_t - IRScov_t*HOUcov_t + ITNcov_t*IRScov_t*HOUcov_t
  # LLINs and IRS
  cITN_IRS <- ITNcov_t*IRScov_t - ITNcov_t*IRScov_t*HOUcov_t
  # LLINs and HOU
  cITN_HOU <- ITNcov_t*HOUcov_t - ITNcov_t*IRScov_t*HOUcov_t
  # IRS and HOU
  cIRS_HOU <- IRScov_t*HOUcov_t - ITNcov_t*IRScov_t*HOUcov_t
  # LLIN and IRS and HOU
  cCom <- ITNcov_t*IRScov_t*HOUcov_t
  # None of them
  c0 <- 1 - ITNcov_t - IRScov_t - HOUcov_t + ITNcov_t*IRScov_t + ITNcov_t*HOUcov_t + IRScov_t*HOUcov_t- ITNcov_t*IRScov_t*HOUcov_t

  ## ---- Probability of repeating a feeding attempt when encountering several interventions
  # LLINs and IRS
  r_LLIN_IRS <- rIRS + (1-rIRS)*rITN
  # LLINs and HOU
  r_LLIN_HOU <- rHOU + (1-rHOU)*rITN
  # IRS and HOU
  r_IRS_HOU <- rHOU + (1-rHOU)*rIRS
  # LLINs, IRS and HOU
  rCom <- rHOU + (1-rHOU)*rIRS + (1-rHOU)*(1-rIRS)*rITN

  # Human - new search probability after a mosquito is repelled (SK) by indoor interventions
  zCom_Human <- Q0*(cHOU*phiI*rHOU + cITN*phiB*rITN + cIRS*phiI*rIRS + cITN_HOU*((phiI-phiB)*rHOU+phiB*r_LLIN_HOU)
                    + cIRS_HOU*phiI*r_IRS_HOU + cITN_IRS*((phiI-phiB)*rIRS+phiB*r_LLIN_IRS)
                    + cCom*((phiI-phiB)*r_IRS_HOU + phiB*rCom))

  ## ---- Probability that a surviving mosquito succeeds in feeding during a single attempt
  # Successfully feeding in presence of LLINs and HM
  s_ITN_HOU <- (1-rHOU)*sITN*sHOU
  # Successfully feeding in presence of LLINs and IRS
  s_ITN_IRS <- (1-rIRS)*sITN*sIRS
  # Successfully feeding in presence of IRS and HM
  s_IRS_HOU <- (1-rHOU)*(1-rIRS)*sIRS*sHOU
  # dIRS = death due to IRS and dHOU = death due to HM: required to normalise sIRS and sHOU
  # since the killing effect occurs later.
  # Feeding in presence of LLINs and IRS and HM
  sCom  <- (1-rHOU)*(1-rIRS)*sITN*(sIRS/(sIRS+dIRS))*(sHOU/(sHOU+dHOU))

  # Probability of successfully feeding in presence of one or all of these interventions
  wCom_Human <- Q0*(cHOU*((1-phiI)+phiI*sHOU)+cITN*((1-phiB)+phiB*sITN)+cIRS*((1-phiI)+phiI*sIRS)
                    + cITN_HOU*((1-phiI)+(phiI-phiB)*sHOU+phiB*s_ITN_HOU)
                    + cITN_IRS*((1-phiI)+(phiI-phiB)*sIRS+phiB*s_ITN_IRS)
                    + cIRS_HOU*((1-phiI)+(phiI )*s_IRS_HOU)
                    + cCom*((1-phiI)+(phiI-phiB)*s_IRS_HOU+phiB*sCom)   )

  # Extracting zcom and wcom
  impactIndoor <- c(zCom_Human,wCom_Human,c0)

  return(impactIndoor)
}


##*********Protecting humans outdoor **********************************************##
impactOutdoorProtection = function(time,time_SPR_on,SPRcov,time_PPM_on,PPMcov,rSPR,rPPM,sSPR,sPPM,Q0_t_h,Q0,phiI,c0,time_OBT_on){
  #. impactOutdoorProtection: compute both zcom and wcom for outdoor protection
  #  (spatial repellents SPR and personal protection measures PPM)
  #  c0 = share of the population with no indoor intervention (from impactIndoorProtection)

  # Incorporate odor baited traps
  if (time > time_OBT_on) Q0 <- Q0_t_h

  SPRcov_t <- activeCoverage(time, SPRcov, time_SPR_on)
  PPMcov_t <- activeCoverage(time, PPMcov, time_PPM_on)

  cov <- splitTwoCoverages(SPRcov_t, PPMcov_t)
  cSPR         <- cov$onlyA     # spatial repellent only
  cPPM         <- cov$onlyB     # personal protection measure only
  cCom_Outdoor <- cov$both      # both applied
  c0_Outdoor   <- cov$neither   # neither applied

  # Repeating due to encountering any
  r_SPR_PPM <- rSPR + (1-rSPR)*rPPM

  # Human - new search probability after a mosquito is repelled by SPR and/or PPM
  zCom_Outdoor <- Q0*c0*(cSPR*rSPR+cPPM*rPPM+cCom_Outdoor*r_SPR_PPM)

  # Successfully feeding in presence of both interventions
  s_SPR_PPM <- (1-rSPR)*sSPR*sPPM

  # (older version: wCom_Outdoor <- Q0*(1-phiI)*(c0_Outdoor+cSPR*sSPR+cPPM*sPPM+cCom_Outdoor*s_SPR_PPM))
  wCom_Outdoor <- Q0*c0*(c0_Outdoor+cSPR*sSPR+cPPM*sPPM+cCom_Outdoor*s_SPR_PPM)

  impactOutdoor <- c(zCom_Outdoor,wCom_Outdoor)

  return(impactOutdoor)
}


##**********************Protection by insecticide treated cattle **********************##
impactInsecticideTreatedCattle = function(time,time_ECS_on,ECScov,time_ECT_on,ECTcov,rECT,sECS,sECT,Q0_t_c,Q0,time_OBT_on){
  #. impactInsecticideTreatedCattle: compute both zcom and wcom for insecticide treated cattle (systemic and topical)

  # Incorporate odor baited traps by changing the proportion of bites on cattle when OBT is on
  if (time > time_OBT_on) { Q0_t_C <- Q0_t_c } else { Q0_t_C <- 1 - Q0 }

  ECScov_t <- activeCoverage(time, ECScov, time_ECS_on)
  ECTcov_t <- activeCoverage(time, ECTcov, time_ECT_on)

  cov <- splitTwoCoverages(ECScov_t, ECTcov_t)
  cECS        <- cov$onlyA     # systemic insecticide only
  cECT        <- cov$onlyB     # topical insecticide only
  cCom_Cattle <- cov$both      # both applied
  c0_Cattle   <- cov$neither   # neither applied

  # Repeating due to encountering insecticide (topical applied (only)) treated cattle.
  # Fix this with JM - need to consider repellency from cattle:
  #   zCom_Cattle <- (1 - Q0)*(cECT*rECT+cCom_Cattle*rECT)
  #   zCom_Cattle <- Q0_t_C*(cECT*rECT+cCom_Cattle*rECT)
  zCom_Cattle = 0

  # Cattle treated with endectocide (systemic and/or topical applied)
  wCom_Cattle <- Q0_t_C*(c0_Cattle+cECS*sECS+cECT*sECT+cCom_Cattle*sECS*sECT)

  impactCattle <- c(zCom_Cattle,wCom_Cattle)

  return(impactCattle)
}


##**************************Targeting Resting, ovipositing mosquitoes**********************##
impactRestingOvipositing = function(time,time_OVI_on,OVIcov,time_ATSB_on,ATSBcov,SSPcov,time_SSP_on,
                                    fOVI,fATSB,fSSP,muV,eOVI){
  #. impactRestingOvipositing: compute the impact of ATSB, Space Spraying, and ovitraps on muV_2

  ##*************************Resting and ovipositing *********************************##
  OVIcov_t  <- activeCoverage(time, OVIcov*eOVI, time_OVI_on)   # ovitrap coverage x trap effectiveness
  ATSBcov_t <- activeCoverage(time, ATSBcov, time_ATSB_on)
  SSPcov_t  <- activeCoverage(time, SSPcov, time_SSP_on)

  # Coverage with ovitraps only
  cOVI <- OVIcov_t - OVIcov_t*ATSBcov_t - OVIcov_t*SSPcov_t + OVIcov_t*ATSBcov_t*SSPcov_t
  # Coverage with ATSB only
  cATSB <- ATSBcov_t - ATSBcov_t*OVIcov_t - ATSBcov_t*SSPcov_t + OVIcov_t*ATSBcov_t*SSPcov_t
  # Coverage with space spraying only
  cSSP <- SSPcov_t - OVIcov_t*SSPcov_t - ATSBcov_t*SSPcov_t + OVIcov_t*ATSBcov_t*SSPcov_t
  # OVI and ATSB
  cOVI_ATSB <- OVIcov_t*ATSBcov_t - OVIcov_t*ATSBcov_t*SSPcov_t
  # OVI and SSP
  cOVI_SSP <- OVIcov_t*SSPcov_t - OVIcov_t*ATSBcov_t*SSPcov_t
  # ATSB and SSP
  cATSB_SSP <- ATSBcov_t*SSPcov_t - OVIcov_t*ATSBcov_t*SSPcov_t
  # OVI and ATSB and SSP
  cCom_RestingOvipositing <- OVIcov_t*ATSBcov_t*SSPcov_t
  # None of them
  # NOTE (unchanged): the third pairwise term below is ATSBcov_t*OVIcov_t, i.e. it repeats the
  # OVI x ATSB term; an exact inclusion-exclusion would use ATSBcov_t*SSPcov_t here.
  c0_RestingOvipositing <- 1 - OVIcov_t - ATSBcov_t - SSPcov_t + OVIcov_t*ATSBcov_t + OVIcov_t*SSPcov_t + ATSBcov_t*OVIcov_t- OVIcov_t*ATSBcov_t*SSPcov_t

  # Change this to accommodate the effectiveness of the trap??
  muV_2_Com = muV*((cOVI*fOVI) + cATSB*fATSB + cSSP*fSSP + cOVI_ATSB*fATSB*(cOVI*fOVI) + cOVI_SSP *(cOVI*fOVI)*fSSP + cATSB_SSP*fATSB*fSSP
                   + cCom_RestingOvipositing *(cOVI*fOVI)*fATSB*fSSP + c0_RestingOvipositing)

  return(muV_2_Com)
}


### Notes for Sam ####
# 1. Check intersection of the three interventions - ITNcov_t*IRScov_t*HOUcov_t
