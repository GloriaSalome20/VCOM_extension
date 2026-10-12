########################################################################
# ODETransmissionParameters.R
# Disease-transmission parameters (human side of the model).
########################################################################

#------------------------------------------------------------------------
# ADDITIONAL TRANSMISSION PARAMETERS
#------------------------------------------------------------------------
#* f0:        Daily biting rate by mosquitoes on animals and humans
#* epsilon0:  Daily entomological inoculation rate
#* recRate:   Human recovery rate
#* NH_eq:     Equilibrium human population size
#* bV:        Probability of transmission from human to vector per infectious bite
#* bh:        Transmission efficiency from an infectious mosquito to an
#*            uninfected, susceptible human
#------------------------------------------------------------------------
getAdditionalTransmissionParameters = function(
  # Low / medium / high transmission is set through epsilon0
  # (e.g. Dickson's change: epsilon0 = 100/365)
  f0 = 1/3, epsilon0 = 50/365, recRate = 1/50,
  NH_eq = 2000, bV = 0.05, bh = 0.5   # bV = originally 0.05
){
  #. getAdditionalTransmissionParameters: Returns the transmission parameters
  c(
    f0 = f0, epsilon0 = epsilon0,
    recRate = recRate, NH_eq = NH_eq,
    bV = bV, bh = bh
  )
}
