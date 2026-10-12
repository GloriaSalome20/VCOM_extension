########################################################################
# ODEModelOutput.R
# Functions that compute model outputs (gonotrophic cycle, biting rate,
# vectorial capacity, EIR, R0, ...).
# Need to test them - then delete this function (note from the original)
########################################################################

######################################################################################
computeLengthGonotrophicycle = function(deltaCom){
  #. computeLengthGonotrophicycle: Length of the gonotrophic cycle in presence of interventions
  f_theta <- 1/deltaCom
  return(f_theta)
}

computeHumanBitingRate = function(f_theta,HBI_com){
  #. computeHumanBitingRate: Human biting rate per mosquito in presence of interventions
  a_theta <- HBI_com/f_theta
  return(a_theta)
}
######################################################################################
computeMosqDensity = function(NV,NH){
  #. computeMosqDensity: Mosquito density in presence of interventions
  Mdensity <- NV/NH
  return(Mdensity)
}
######################################################################################
computeHBI = function(wCom_human,wCom){
  #. computeHBI: Human blood index
  #  HBI_com <- (Pfi+pfo)/pf
  HBI_com <- (wCom_human)/wCom
  return(HBI_com)
}

######################################################################################
# NOTE: computeVC() and computeRO() read recRate from the global `theta`, and durEV too unless the
# incubation period `tilda` is given (IVM_ode() gives it, because it can change with the temperature).
computeVC = function(a_theta, NV,NH,muVCom, tilda = theta[["durEV"]]){
  #. computeVC: Vectorial capacity (tilda = incubation period)
  VC <- ((NV/NH) * a_theta^2 *exp(-muVCom*tilda))/muVCom
  return(VC)
}
######################################################################################
computeEIR = function(a_theta, IV, NH){
  #. computeEIR: Entomological inoculation rate
  EIR <- 365*(IV*a_theta)/NH    # 365 per year
  return(EIR)
}

######################################################################################
computeRO = function(a_theta,muVCom, NV,bv,bh,NH, tilda = theta[["durEV"]]){
  #. computeRO: Basic reproduction number (tilda = incubation period)
  R_o <- (NV/NH)*a_theta^2*bv*(bh*exp(-muVCom*tilda)/(theta[["recRate"]]*muVCom))
  return(R_o)
}
