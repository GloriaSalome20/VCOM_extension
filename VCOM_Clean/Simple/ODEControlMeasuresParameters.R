########################################################################
# ODEControlMeasuresParameters.R
# Intervention parameters: coverage of each control measure and the day
# (time_XXX_on) on which it is switched on.
########################################################################

#------------------------------------------------------------------------
# INTERVENTION PARAMETERS
#------------------------------------------------------------------------
# For every intervention XXX:  XXXcov = coverage (0-1),  time_XXX_on = start day.
#
#* OVI   Ovitraps
#* FOG   Fogging
#* LAR   Larvicide
#* BIO   Biological control
#* SRE   Source reduction
#* IRS   Indoor residual spraying
#* ITN   Insecticide-treated nets (LLINs)
#* ECS   Endectocide in cattle, systemic
#* ECT   Endectocide in cattle, topical
#* HOU   House modification (mosquito-proofed housing)
#* OBT   Odor-baited traps (coverage = ratio of traps to humans)
#* SPR   Spatial repellents
#* PPM   Personal protection measures
#* ATSB  Attractive toxic sugar baits
#* SSP   Space spraying
#------------------------------------------------------------------------
getInterventionsParameters = function(
  OVIcov=0,   time_OVI_on=0,
  FOGcov=0,   time_FOG_on=0,
  LARcov=0,   time_LAR_on=0,
  BIOcov=0,   time_BIO_on=0,
  SREcov=0,   time_SRE_on=0,
  IRScov=0,   time_IRS_on=0,
  ITNcov=0,   time_ITN_on=0,
  ECScov=0,   time_ECS_on=0,
  ECTcov=0,   time_ECT_on=0,
  HOUcov=0,   time_HOU_on=0,
  OBTcov=0,   time_OBT_on=0,
  SPRcov=0,   time_SPR_on=0,
  PPMcov=0,   time_PPM_on=0,
  ATSBcov=0,  time_ATSB_on=0,
  SSPcov=0,   time_SSP_on=0
){
  #. getInterventionsParameters: Returns the intervention coverages and the
  #  time (day) at which each intervention is switched on (SK)
  c(
    OVIcov=OVIcov,   time_OVI_on=time_OVI_on,
    FOGcov=FOGcov,   time_FOG_on=time_FOG_on,
    LARcov=LARcov,   time_LAR_on=time_LAR_on,
    BIOcov=BIOcov,   time_BIO_on=time_BIO_on,
    SREcov=SREcov,   time_SRE_on=time_SRE_on,
    IRScov=IRScov,   time_IRS_on=time_IRS_on,
    ITNcov=ITNcov,   time_ITN_on=time_ITN_on,
    ECScov=ECScov,   time_ECS_on=time_ECS_on,
    ECTcov=ECTcov,   time_ECT_on=time_ECT_on,
    HOUcov=HOUcov,   time_HOU_on=time_HOU_on,
    OBTcov=OBTcov,   time_OBT_on=time_OBT_on,
    SPRcov=SPRcov,   time_SPR_on=time_SPR_on,
    PPMcov=PPMcov,   time_PPM_on=time_PPM_on,
    ATSBcov=ATSBcov, time_ATSB_on=time_ATSB_on,
    SSPcov=SSPcov,   time_SSP_on=time_SSP_on
  )
}
