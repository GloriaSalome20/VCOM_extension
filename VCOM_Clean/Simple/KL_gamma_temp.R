################################################################################
##  KL_gamma_temp.R
##
##  Recomputes the seasonal larval carrying capacity K_L(t) for Tanga,
##  Tanzania using:
##
##    1. Gamma-distributed rainfall memory kernel (Beck-Johnson et al. 2013)
##         R_tilde(t) = sum_{tau=0}^{tau_max}  w_tau * R(t - tau)
##         w_tau  ~  Gamma(a, b)   (discrete, normalised to sum = 1)
##       This produces a LAG -- weight peaks ~(a-1)*b days before today,
##       reflecting the delay between rainfall and available breeding sites.
##
##    2. Briere-type temperature suitability (Mordecai et al. 2013)
##         g(T) = c_T * T * (T - T_min) * sqrt(T_max - T)   for T_min < T < T_max
##         g(T) = 0                                           otherwise
##       Parameterised for Anopheles gambiae larval development rate.
##
##    K_L(t) = K0 * R_tilde(t) * g_norm(T(t))
##
##  Replaces the White et al. (2011) exponential kernel previously used in
##  carrying_capacity_tanga_rain_temp.csv.
##
##  Inputs (must exist in the working directory = this folder):
##    isimip3b_tanga_climatology.csv -- 365 days of rainfall (mm/day) and temperature (deg C) of
##                                      one ISIMIP3b climate model, or of the mean of the models
##                                      (made by buildClimatology() in ISIMIP3b_climate.R)
##    The model is chosen with ISIMIP_MODEL (set it before sourcing this file; default "ensemble").
##    This replaces rain.csv and monthly_tmp_2019_2020.rds / Temperature_daily_converted.R (CRU).
##
##  Output:
##    K_FILE (default carrying_capacity_tanga_gamma_isimip.csv) -- 365-row K_L with columns:
##                                           Day, K, rain_memory, temp_suitability
##
##  Run from this folder as working directory.
##
##  ISIMIP3b version: only sections 1 and 2 (the inputs) and the output name in section 8 differ
##  from the CRU version; the kernel, the Briere function and the scaling are unchanged.
##
##  Author: Gloria Salome Shirima  (IHI / NM-AIST)
##  Reference: Beck-Johnson et al. (2013) PLoS ONE; Mordecai et al. (2013) Ecol Lett
################################################################################

## ---- 0. Packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## ---- 1. Daily rainfall from ISIMIP3b ---------------------------------------
## Mean annual cycle of daily rainfall (mm/day) of the chosen climate model for the Tanga cell
source("ISIMIP3b_climate.R")
if (!exists("ISIMIP_MODEL")) ISIMIP_MODEL <- "ensemble"
if (!exists("CLIMATE"))      CLIMATE <- loadClimate(ISIMIP_MODEL)   # day, tas, pr, hurs (set in Driver.R)
R_vec   <- CLIMATE$pr[1:365]
R_vec   <- pmax(R_vec, 0)
cat(sprintf("Rainfall loaded (ISIMIP3b, %s): %d days, mean = %.2f mm/day, max = %.2f mm/day\n",
            ISIMIP_MODEL, length(R_vec), mean(R_vec), max(R_vec)))

## ---- 2. Daily temperature from ISIMIP3b -------------------------------------
## Same climate model, same 365 days: Day 1 = 1 January
T_vec <- CLIMATE$tas[1:365]
cat(sprintf("Temperature loaded (ISIMIP3b, %s): %d days, mean = %.1f C, range = [%.1f, %.1f] C\n",
            ISIMIP_MODEL, length(T_vec), mean(T_vec), min(T_vec), max(T_vec)))

## ---- 3. Gamma kernel parameters (Beck-Johnson et al. 2013) -----------------
## Shape a and scale b control where the peak weight falls:
##   peak lag = (a-1)*b  days before today
##   a=3, b=4 => peak at 8 days ago  (moderate lag)
##   a=4, b=4 => peak at 12 days ago (longer lag)
##   a=3, b=3 => peak at 6 days ago  (shorter lag)
## tau_max: maximum days of memory included (set >= 3*(a-1)*b)

a_gamma  <- 3      # gamma shape
b_gamma  <- 4      # gamma scale (days)
tau_max  <- 21     # memory window (days)

## Discrete gamma weights (tau = 0, 1, ..., tau_max)
tau_grid     <- 0:tau_max
gamma_weights <- (tau_grid / b_gamma)^(a_gamma - 1) * exp(-tau_grid / b_gamma)
## Avoid NaN at tau=0 when a<1 (not an issue here since a=3)
gamma_weights[is.nan(gamma_weights)] <- 0
gamma_weights <- gamma_weights / sum(gamma_weights)  # normalise to sum = 1

cat(sprintf("\nGamma kernel: a=%g, b=%g, tau_max=%d days\n",
            a_gamma, b_gamma, tau_max))
cat(sprintf("  Peak weight at lag: %d days\n", tau_grid[which.max(gamma_weights)]))
cat(sprintf("  Weights sum to: %.6f\n", sum(gamma_weights)))

## ---- 4. Compute R_tilde(t) -- gamma-weighted rainfall memory ---------------
## The memory of Day 1 needs the rain of the tau_max days before Day 1. RAIN_PADDING (set in Driver.R) chooses them:
##   "circular": the last tau_max days of the cycle, because the climatology repeats every year and Day 1 follows Day 365
##               (the model repeats the year, so K(365) and K(1) then join smoothly)
##   "zero":     no rain before Day 1, as in the CRU version, where one real year of rain was used
if (!exists("RAIN_PADDING")) RAIN_PADDING <- "circular"
stopifnot(RAIN_PADDING %in% c("circular", "zero"))
R_padded <- c(if (RAIN_PADDING == "zero") rep(0, tau_max) else tail(R_vec, tau_max), R_vec)   # length = 365 + tau_max
cat(sprintf("Rain before Day 1: %s\n", RAIN_PADDING))

R_tilde <- numeric(365)
for (t in 1:365) {
  ## tau=0 => today (index t + tau_max), tau=tau_max => most lagged
  R_window <- R_padded[(t + tau_max):(t + tau_max - tau_max)]
  ## reverse: R_window[1] = R(t), R_window[tau_max+1] = R(t - tau_max)
  R_tilde[t] <- sum(gamma_weights * R_window)
}

cat(sprintf("\nR_tilde range: [%.4f, %.4f] mm/day (weighted)\n",
            min(R_tilde), max(R_tilde)))
cat(sprintf("R_tilde peak day: %d  (rainfall peak day: %d)  => lag: %d days\n",
            which.max(R_tilde), which.max(R_vec),
            which.max(R_tilde) - which.max(R_vec)))

## ---- 5. Briere temperature suitability (Mordecai et al. 2013) --------------
## For An. gambiae larval development rate (MDR):
##   T_min = 14.58 C,  T_max = 34.61 C,  c_T = 0.000111 (normalised out below)
## We use the shape only; c_T cancels after normalisation.

T_min_briereE <- 14.58    # lower thermal limit (°C)
T_max_briereE <- 34.61    # upper thermal limit (°C)

briere <- function(T, T_min, T_max) {
  g <- ifelse(T > T_min & T < T_max,
              T * (T - T_min) * sqrt(T_max - T),
              0)
  pmax(g, 0)
}

g_T_raw <- briere(T_vec, T_min_briereE, T_max_briereE)

## Normalise to [0, 1] so K0 has a simple biological meaning
g_T_norm <- g_T_raw / max(g_T_raw, na.rm = TRUE)
g_T_norm[is.nan(g_T_norm)] <- 0

cat(sprintf("\nBriere g(T): peak at T=%.1f C (analytical optimum)\n",
            T_vec[which.max(g_T_norm)]))
cat(sprintf("  g(T) days > 0.9 * max: %d days\n",
            sum(g_T_norm > 0.9, na.rm = TRUE)))

## ---- 6. Normalise R_tilde to [0, 1] ----------------------------------------
R_tilde_norm <- R_tilde / max(R_tilde, na.rm = TRUE)

## ---- 7. Combined K_L(t) = K0 * R_tilde(t) * g(T(t)) ----------------------
##  K0 calibrated to match L_eq (equilibrium larval carrying capacity)
##  These values should be tuned to your calibrated mosquito population.
L_eq          <- 120000     # target maximum K_L (larvae/ha) -- match your calibration
constant_frac <- 0.05       # minimum K_L as fraction of L_eq (dry season floor)

K_raw <- R_tilde_norm * g_T_norm
## Rescale: K_L ranges from constant_frac*L_eq (dry season) to L_eq (peak)
K_L <- K_raw / max(K_raw, na.rm = TRUE) * L_eq * (1 - constant_frac) +
       L_eq * constant_frac

cat(sprintf("\nK_L range: [%.0f, %.0f]\n", min(K_L), max(K_L)))
cat(sprintf("K_L peak day  : %d\n", which.max(K_L)))
cat(sprintf("R_tilde peak  : %d\n", which.max(R_tilde)))
cat(sprintf("Rainfall peak : %d\n", which.max(R_vec)))
cat(sprintf("g(T) peak day : %d\n", which.max(g_T_norm)))

## ---- 8. Save output ---------------------------------------------------------
out_df <- data.frame(
  Day              = 1:365,
  K                = K_L,
  rain_memory      = R_tilde,
  temp_suitability = g_T_norm
)

K_OUT <- if (exists("K_FILE")) K_FILE else "carrying_capacity_tanga_gamma_isimip.csv"   # K_FILE is set in Driver.R
write.csv(out_df, K_OUT, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", K_OUT))