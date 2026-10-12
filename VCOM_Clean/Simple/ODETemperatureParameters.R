########################################################################
# ODETemperatureParameters.R
# Temperature- and humidity-dependent mosquito parameters
# (thesis Sec. "Temperature- and Humidity-Dependent Mosquito Parameters").
#
# Each function takes the temperature T (deg C) and returns one of the quantities of the
# material. temperatureDependentTheta() puts them into the parameter vector theta, group by
# group, and makeTemperatureDependentParameters() turns them into the function that
# IVM_ode() calls at every time step (timeVaryingParameters(), see ODEModel.R).
#
#   group               what it sets in theta            function of the material
#   adult_mortality     muV                              Martens (1997) g(T), p = exp(-1/g(T)); humidity S(h)
#   biting_rate         tau1, tau2 (and f0)              b(T), quadratic (Kalula et al. 2023)
#   eip                 durEV                            n(T) = k/(T - T0)
#   aquatic_mortality   muEL, muLL, muPL                 c0 + c1 (T - 25)^2   (Bayoh 2003)
#   development         durEL, durLL, durPL              1/[q T (T - T0) sqrt(Tm - T)]  (Briere)
#   oviposition         beta                             q T (T - T0) sqrt(Tm - T) x eggs per batch
#
# Everything else in theta (gamma, Q0, the intervention parameters, bV, bh, ...) stays constant.
########################################################################

# 1. PARAMETERS OF THE MATERIAL ---------------------------------------------------------------

# Adult mortality: g(T) = a0 + a1 T + a2 T^2 (Martens 1997)
ADULT_MORTALITY_QUADRATIC <- c(a0 = -4.4, a1 = 1.31, a2 = -0.03)
# Biting rate (per day): b(T) = a0 + a1 T + a2 T^2 (Kalula et al. 2023)
BITING_RATE_QUADRATIC     <- c(a0 = -0.322, a1 = 0.027, a2 = -0.00014)

# Humidity: S(h) is 0 up to HUMIDITY_MIN (%), 1 from HUMIDITY_MAX (%), linear in between (Yamana & Eltahir 2013)
HUMIDITY_MIN <- 5
HUMIDITY_MAX <- 42

# Extrinsic incubation period: n(T) = k/(T - T0) days (Macdonald 1957, Detinova 1962, Gething 2011)
EIP_K  <- 111
EIP_T0 <- 16

# Aquatic-stage mortality: c0 + c1 (T - Tref)^2 (Bayoh 2003)
AQUATIC_REFERENCE_TEMPERATURE  <- 25
AQUATIC_MORTALITY_COEFFICIENTS <- list(
  egg   = c(c0 = 0.010, c1 = 0.0020),
  larva = c(c0 = 0.020, c1 = 0.0030),
  pupa  = c(c0 = 0.015, c1 = 0.0025)
)

# Briere function q T (T - T0) sqrt(Tm - T) between T0 and Tm, 0 outside (Mordecai 2013, Beck-Johnson 2013)
BRIERE_COEFFICIENTS <- list(
  oviposition = c(q = 0.000200, T0 = 10, Tm = 40),   # oviposition rate
  egg         = c(q = 0.000230, T0 = 10, Tm = 40),   # egg development rate, 1/d_E
  larva       = c(q = 0.000103, T0 = 12, Tm = 36),   # larva development rate, 1/d_L
  pupa        = c(q = 0.000464, T0 = 12, Tm = 36)    # pupa development rate, 1/d_P
)

# NOT in the material. The oviposition Briere function is a rate per day (about 0.3), while the model's
# beta is the number of EGGS per female per day (21.19 in the October model). The multiplier 63 is the one
# used for beta in your earlier script (ODEModel_season_rounds_choice_draft2_noseason.R); at 27.3 deg C it
# gives back the model's beta = 21.19.
EGGS_PER_BATCH <- 63

# NOT in the material. Keeps the model finite where a function is zero or undefined (outside the thermal
# limits): the largest adult mortality (per day) and the longest duration (days) that are used.
MAX_ADULT_MORTALITY <- 10
MAX_DURATION        <- 100

# The groups above, in the order of the table at the top of this file
TEMPERATURE_GROUPS <- c("adult_mortality", "biting_rate", "eip", "aquatic_mortality", "development", "oviposition")


# 2. THE FUNCTIONS OF THE MATERIAL ------------------------------------------------------------

quadratic = function(T, coefficients){
  #. quadratic: a0 + a1 T + a2 T^2
  coefficients[[1]] + coefficients[[2]]*T + coefficients[[3]]*T^2
}

briereRate = function(T, q, T0, Tm){
  #. briereRate: q T (T - T0) sqrt(Tm - T) for T0 < T < Tm, 0 outside
  ifelse(T > T0 & T < Tm, q*T*(T - T0)*sqrt(pmax(Tm - T, 0)), 0)
}

humiditySurvival = function(h){
  #. humiditySurvival: S(h), the humidity modifier of the daily survival (h = relative humidity in %)
  pmin(pmax((h - HUMIDITY_MIN)/(HUMIDITY_MAX - HUMIDITY_MIN), 0), 1)
}

adultMortality = function(T, h = NULL){
  #. adultMortality: daily adult death rate muV(T, h) = -log p(T, h), with p(T, h) = exp(-1/g(T)) S(h).
  #.   g(T) is the quadratic of the material (about the expected adult lifespan in days). The material writes
  #.   the daily survival as exp(-g(T)/12) S(h), which does not follow Martens (1997); exp(-1/g(T)) is used here.
  #.   Without humidity (h = NULL) S = 1. Where g(T) <= 0 (outside about 3.7 to 40 deg C) or S = 0 the
  #.   mortality is MAX_ADULT_MORTALITY.
  g  <- quadratic(T, ADULT_MORTALITY_QUADRATIC)
  mu <- ifelse(g > 1/MAX_ADULT_MORTALITY, 1/g, MAX_ADULT_MORTALITY)
  if (!is.null(h)) mu <- mu - log(humiditySurvival(h))
  pmin(mu, MAX_ADULT_MORTALITY)
}

bitingRate = function(T){
  #. bitingRate: b(T), bites per mosquito per day (the inverse of the gonotrophic cycle), not below 0
  pmax(quadratic(T, BITING_RATE_QUADRATIC), 0)
}

gonotrophicCycle = function(T){
  #. gonotrophicCycle: 1/b(T), days from one blood meal to the next (tau1 + tau2)
  pmin(1/bitingRate(T), MAX_DURATION)
}

extrinsicIncubationPeriod = function(T){
  #. extrinsicIncubationPeriod: n(T) = k/(T - T0) days, at most MAX_DURATION (n is undefined for T <= T0)
  pmin(EIP_K/pmax(T - EIP_T0, 1e-9), MAX_DURATION)
}

aquaticMortality = function(T, stage = c("egg", "larva", "pupa")){
  #. aquaticMortality: daily mortality c0 + c1 (T - 25)^2 of the egg, larval or pupal stage
  coefficients <- AQUATIC_MORTALITY_COEFFICIENTS[[match.arg(stage)]]
  coefficients[["c0"]] + coefficients[["c1"]]*(T - AQUATIC_REFERENCE_TEMPERATURE)^2
}

developmentDuration = function(T, stage = c("egg", "larva", "pupa")){
  #. developmentDuration: d(T) = 1/[Briere rate] days of the egg, larva or pupa stage, at most MAX_DURATION
  b    <- BRIERE_COEFFICIENTS[[match.arg(stage)]]
  rate <- briereRate(T, b[["q"]], b[["T0"]], b[["Tm"]])
  pmin(1/pmax(rate, 1/MAX_DURATION), MAX_DURATION)
}

eggsPerFemalePerDay = function(T){
  #. eggsPerFemalePerDay: beta(T) = EGGS_PER_BATCH x Briere oviposition rate (see EGGS_PER_BATCH)
  b <- BRIERE_COEFFICIENTS[["oviposition"]]
  EGGS_PER_BATCH*briereRate(T, b[["q"]], b[["T0"]], b[["Tm"]])
}


# 3. FROM TEMPERATURE TO THE PARAMETER VECTOR theta -----------------------------------------------

temperatureDependentTheta = function(theta, T, h = NULL, use = TEMPERATURE_GROUPS){
  #. temperatureDependentTheta: returns theta with the parameters of the groups in `use` set to their
  #.   values at temperature T (and relative humidity h, if given). The values are absolute: they replace
  #.   the constants, so applying the function again to its own result changes nothing.
  stopifnot(all(use %in% TEMPERATURE_GROUPS))
  if ("adult_mortality" %in% use) theta[["muV"]] <- adultMortality(T, h)
  if ("biting_rate" %in% use) {
    # The gonotrophic cycle tau1 + tau2 becomes 1/b(T); the share spent foraging (tau1) is kept
    cycle_ratio  <- gonotrophicCycle(T)/(theta[["tau1"]] + theta[["tau2"]])
    theta[["tau1"]] <- theta[["tau1"]]*cycle_ratio
    theta[["tau2"]] <- theta[["tau2"]]*cycle_ratio
    theta[["f0"]]   <- bitingRate(T)   # daily biting rate, used for the equilibrium (calculateInitialState)
  }
  if ("eip" %in% use) theta[["durEV"]] <- extrinsicIncubationPeriod(T)
  if ("aquatic_mortality" %in% use) {
    theta[["muEL"]] <- aquaticMortality(T, "egg")
    theta[["muLL"]] <- aquaticMortality(T, "larva")
    theta[["muPL"]] <- aquaticMortality(T, "pupa")
  }
  if ("development" %in% use) {
    theta[["durEL"]] <- developmentDuration(T, "egg")
    theta[["durLL"]] <- developmentDuration(T, "larva")
    theta[["durPL"]] <- developmentDuration(T, "pupa")
  }
  if ("oviposition" %in% use) theta[["beta"]] <- eggsPerFemalePerDay(T)
  theta
}

makeTemperatureDependentParameters = function(temperature, humidity = NULL, use = TEMPERATURE_GROUPS){
  #. makeTemperatureDependentParameters: builds timeVaryingParameters(time, theta) for IVM_ode().
  #.   temperature and humidity are functions of the model time (see makeDailyCycle() in ODEModel.R);
  #.   humidity = NULL means that humidity is not limiting (S = 1).
  force(temperature); force(humidity); force(use)
  function(time, theta){
    temperatureDependentTheta(theta, temperature(time), if (is.null(humidity)) NULL else humidity(time), use)
  }
}


# 4. CHECK PLOT --------------------------------------------------------------------------------

plotTemperatureFunctions = function(temperatures = seq(5, 42, by = 0.25)){
  #. plotTemperatureFunctions: the functions of the material against temperature (a check, not used by the model)
  op <- par(mfrow = c(3, 3), mar = c(4, 4.5, 2, 1)); on.exit(par(op))
  panel = function(y, main, ylab, ylim = NULL) plot(temperatures, y, type = "l", main = main, xlab = "Temperature (deg C)", ylab = ylab, ylim = ylim)
  panel(adultMortality(temperatures),                    "Adult mortality",             "muV (per day)",  c(0, 1))
  panel(bitingRate(temperatures),                        "Biting rate",                 "b(T) (per day)")
  panel(extrinsicIncubationPeriod(temperatures),         "Extrinsic incubation period", "n(T) (days)",    c(0, 40))
  panel(aquaticMortality(temperatures, "egg"),           "Egg mortality",               "per day",        c(0, 0.3))
  panel(aquaticMortality(temperatures, "larva"),         "Larva mortality",             "per day",        c(0, 0.3))
  panel(aquaticMortality(temperatures, "pupa"),          "Pupa mortality",              "per day",        c(0, 0.3))
  panel(eggsPerFemalePerDay(temperatures),               "Eggs per female per day",     "beta(T)")
  panel(developmentDuration(temperatures, "larva"),      "Larva development",           "d_L(T) (days)",  c(0, 40))
  panel(developmentDuration(temperatures, "pupa"),       "Pupa development",            "d_P(T) (days)",  c(0, 40))
}
