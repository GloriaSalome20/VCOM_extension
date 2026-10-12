########################################################################
# SeasonalSetup.R
# The seasonal parts of the model as functions, for the Shiny app (ui.r and server.r).
#
# Driver.R sets up the climate, the seasonal carrying capacity, the temperature-dependent
# parameters and the burn-in in its sections 3 to 6, with the choices written as settings at
# the top. The app needs the same steps but with the choices coming from the sidebar, so they
# are repeated here as functions. They do the same thing as Driver.R (the app and the Driver
# give the same trajectory for the same choices, see ISIMIP_CHANGES.md); if you change one of
# the steps in Driver.R, change it here too.
#
#   setSeasonality()   climate model -> K(t) and the temperature-dependent parameters
#                      (sets carryingCapacity() and timeVaryingParameters(), the two hooks of ODEModel.R)
#   runVCOM()          one complete run: species, EIR, interventions, seasonality, burn-in
#
# Files used: isimip3b_tanga_climatology.csv (the climate), carrying_capacity_tanga_gamma_isimip_<model>.csv
# (K(t) of each climate model, made by KL_gamma_temp.R; it is made again here if one is missing),
# ISIMIP3b_climate.R, ODETemperatureParameters.R.
########################################################################

source("ISIMIP3b_climate.R")
source("ODETemperatureParameters.R")


# 1. CHOICES OF THE SIDEBAR ------------------------------------------------------------------

SEASONAL_CLIMATE_CHOICES <- c("Mean of the five models (ensemble)" = "ensemble",
                              setNames(ISIMIP_MODEL_NAMES, ISIMIP_MODEL_NAMES))

SEASONAL_TEMPERATURE_CHOICES <- c(
  "Adult mortality (and humidity)"       = "adult_mortality",
  "Biting rate"                          = "biting_rate",
  "Extrinsic incubation period"          = "eip",
  "Egg, larva and pupa mortality"        = "aquatic_mortality",
  "Egg, larva and pupa development time" = "development",
  "Eggs laid per female"                 = "oviposition")
stopifnot(setequal(SEASONAL_TEMPERATURE_CHOICES, TEMPERATURE_GROUPS))

# The interventions of the app (FOG is switched off in the app; ATSB and SSP are not in the app)
APP_INTERVENTIONS <- c("OVI", "LAR", "BIO", "SRE", "ITN", "IRS", "ECS", "ECT", "HOU", "OBT", "SPR", "PPM")


# 2. CLIMATE, K(t) AND TEMPERATURE-DEPENDENT PARAMETERS ----------------------------------------------

seasonalKFile <- function(model, padding = "circular"){
  #. seasonalKFile: name of the K(t) file of one climate model (the same name as in Driver.R)
  paste0("carrying_capacity_tanga_gamma_isimip_", climateTag(model), if (padding == "zero") "_zeropad", ".csv")
}

loadSeasonalK <- function(model, climate, padding = "circular"){
  #. loadSeasonalK: the 365 daily values of K(t) of one climate model. The file is made by KL_gamma_temp.R when it does not
  #.   exist or when the climate table is more than a minute newer (as in Driver.R; the minute allows for files copied together).
  k_file <- seasonalKFile(model, padding)
  stale  <- !file.exists(k_file) ||
            as.numeric(difftime(file.mtime(ISIMIP_CLIMATE_FILE), file.mtime(k_file), units = "secs")) > 60
  if (stale) {
    message("making ", k_file, " with KL_gamma_temp.R")
    env <- new.env(parent = globalenv())          # the settings KL_gamma_temp.R reads (set in Driver.R in the Driver)
    env$ISIMIP_MODEL <- model; env$CLIMATE <- climate; env$RAIN_PADDING <- padding; env$K_FILE <- k_file
    made <- tryCatch({ capture.output(suppressMessages(source("KL_gamma_temp.R", local = env))); TRUE },
                     error = function(e) { message("KL_gamma_temp.R failed: ", conditionMessage(e)); FALSE })
    if (!made && !file.exists(k_file)) stop(k_file, " does not exist and KL_gamma_temp.R could not make it (see the message above)")
  }
  read.csv(k_file)$K
}

setSeasonality <- function(model = "ensemble", seasonalK = TRUE, temperatureGroups = TEMPERATURE_GROUPS,
                           useHumidity = TRUE, kScale = "relative", padding = "circular"){
  #. setSeasonality: sets the two hooks of ODEModel.R for the chosen climate (the same as sections 3 to 5 of Driver.R) and
  #.   returns what runVCOM() needs to start the equilibrium at the temperature of day 0:
  #.   list(climate, T_start, h_start, temperatureGroups). The hooks are global, as IVM_ode() reads them from the global environment.
  stopifnot(model %in% SEASONAL_CLIMATE_CHOICES, all(temperatureGroups %in% TEMPERATURE_GROUPS), kScale %in% c("relative", "absolute"))
  G <- globalenv()
  climate <- loadClimate(model)                    # data frame: day, tas (deg C), pr (mm/day), hurs (%)

  ## K(t): the seasonal curve, or the constant K0 of the model before
  if (seasonalK) assign("carryingCapacity", makeSeasonalCarryingCapacity(loadSeasonalK(model, climate, padding), kScale), envir = G)
  else           assign("carryingCapacity", function(time, K0) K0, envir = G)

  ## Temperature- and humidity-dependent parameters: the same 365 days as K(t), Day 1 = model day 0, repeating
  T_start <- h_start <- NULL
  if (length(temperatureGroups) > 0) {
    temperatureYear <- makeDailyCycle(climate$tas)
    humidityYear    <- if (useHumidity && !all(is.na(climate$hurs))) makeDailyCycle(climate$hurs) else NULL
    T_start <- temperatureYear(0)                  # the equilibrium (initial state, NV_eq, lambdaV, K0) is at the temperature of day 0
    h_start <- if (is.null(humidityYear)) NULL else humidityYear(0)
    assign("timeVaryingParameters", makeTemperatureDependentParameters(temperatureYear, humidityYear, temperatureGroups), envir = G)
  } else {
    assign("timeVaryingParameters", function(time, theta) theta, envir = G)     # constant parameters
  }
  list(climate = climate, T_start = T_start, h_start = h_start, temperatureGroups = temperatureGroups)
}


# 3. ONE COMPLETE RUN --------------------------------------------------------------------------------

runVCOM <- function(species = "GAM", EIR = 100, days = 365,
                    coverage = setNames(rep(0, length(APP_INTERVENTIONS)), APP_INTERVENTIONS),
                    start    = setNames(rep(0, length(APP_INTERVENTIONS)), APP_INTERVENTIONS),
                    climateModel = "ensemble", seasonalK = TRUE, temperatureGroups = TEMPERATURE_GROUPS,
                    burnInYears = 1){
  #. runVCOM: runs the model as server.r always did, plus the seasonal choices, and leaves the objects the app needs in the
  #.   global environment (MOSQUITO_PARAMETERS, TRANSMISSION_PARAMETERS, INTERVENTION_PARAMETERS, theta, initState, IVM_traj,
  #.   RUN_SUMMARY). `coverage` and `start` are named by intervention code (APP_INTERVENTIONS); the start days are counted from
  #.   day 0 = 1 January, after the burn-in. An intervention with coverage 0 is left off (see below). Returns the trajectory (invisibly).
  stopifnot(species %in% c("GAM", "ARA", "FUN"), burnInYears >= 0, burnInYears == round(burnInYears))
  G <- globalenv()
  burnInDays <- 365*burnInYears
  ## interventions that are not given have coverage 0 and start day 0
  fillInterventions <- function(x){
    stopifnot(all(names(x) %in% APP_INTERVENTIONS))
    full <- setNames(rep(0, length(APP_INTERVENTIONS)), APP_INTERVENTIONS); full[names(x)] <- x; full
  }
  coverage <- fillInterventions(coverage)
  start    <- fillInterventions(start)
  ## An intervention with coverage 0 is not switched on at all (its start day is moved after the end of the run, as Driver.R does
  ## for the interventions that are off). It should change nothing, but the odor-baited traps formula is not neutral at coverage 0:
  ## from its start day it replaces the human blood index Q0 of the species by 1/(1 + 0.5) (cattle only), which lowers the EIR,
  ## the vectorial capacity and R0 by about 28 % for An. gambiae in one step (see ISIMIP_CHANGES.md).
  start[coverage <= 0] <- days + 1

  season <- setSeasonality(climateModel, seasonalK, temperatureGroups)

  ## Parameters (theta): species, transmission (EIR), interventions (start days moved by the burn-in)
  assign("MOSQUITO_PARAMETERS", switch(species,
    "GAM" = getAnGambiaeParameters(),
    "ARA" = getAnArabiensisParameters(),
    "FUN" = getAnFunestusParameters()), envir = G)
  assign("TRANSMISSION_PARAMETERS", getAdditionalTransmissionParameters(epsilon0 = as.numeric(EIR)/365), envir = G)
  coverage_args <- setNames(as.list(as.numeric(coverage[APP_INTERVENTIONS])),           paste0(APP_INTERVENTIONS, "cov"))
  time_args     <- setNames(as.list(as.numeric(start[APP_INTERVENTIONS]) + burnInDays), paste0("time_", APP_INTERVENTIONS, "_on"))
  assign("INTERVENTION_PARAMETERS", do.call(getInterventionsParameters, c(coverage_args, time_args)), envir = G)
  assign("theta", getTheta(speciesSpecificParameters = get("MOSQUITO_PARAMETERS", envir = G),
                           interventionParameters = get("INTERVENTION_PARAMETERS", envir = G),
                           additionalTransmissionParameters = get("TRANSMISSION_PARAMETERS", envir = G)), envir = G)
  if (length(season$temperatureGroups) > 0)
    assign("theta", temperatureDependentTheta(get("theta", envir = G), season$T_start, season$h_start, season$temperatureGroups), envir = G)

  ## Initial state (calculateInitialState() adds lambdaV and NV_eq to the GLOBAL theta, which is what the solver then uses)
  assign("initState", calculateInitialState(get("theta", envir = G)), envir = G)

  ## Run; the burn-in years are run without interventions and dropped, so day 0 is 1 January in the model's own seasonal state
  traj <- runODE(days + burnInDays, 1, get("initState", envir = G), get("theta", envir = G), "daspk")
  traj <- traj[traj$time >= burnInDays, ]
  traj$time <- traj$time - burnInDays
  rownames(traj) <- NULL
  assign("IVM_traj", traj, envir = G)

  assign("RUN_SUMMARY", describeRun(species, EIR, days, season, seasonalK, burnInYears, climateModel), envir = G)
  invisible(traj)
}

describeRun <- function(species, EIR, days, season, seasonalK, burnInYears, climateModel){
  #. describeRun: one short text for the page, saying what was run
  clim <- season$climate
  groups <- season$temperatureGroups
  paste0(
    "Shown: ", switch(species, "GAM" = "An. gambiae", "ARA" = "An. arabiensis", "FUN" = "An. funestus"),
    ", EIR level ", EIR, ", ", days, " days from 1 January",
    if (burnInYears > 0) paste0(", after ", burnInYears, if (burnInYears == 1) " year" else " years", " of burn-in without interventions") else ", no burn-in (starts at the equilibrium of a constant environment)",
    ". Climate: ISIMIP3b ", if (climateModel == "ensemble") "mean of five models" else climateModel,
    sprintf(" for Tanga (mean %.1f C, %.0f mm rain a year). ", mean(clim$tas), sum(clim$pr)),
    "Larval carrying capacity: ", if (seasonalK) "seasonal K(t) from rainfall and temperature. " else "constant. ",
    "Temperature-dependent parameters: ",
    if (length(groups) == 0) "none (constant)." else paste0(length(groups), " of ", length(TEMPERATURE_GROUPS), " groups on (",
      paste(names(SEASONAL_TEMPERATURE_CHOICES)[match(groups, SEASONAL_TEMPERATURE_CHOICES)], collapse = ", "), ")."))
}
