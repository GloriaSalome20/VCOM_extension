########################################################################
# ISIMIP3b_climate.R
# Daily climate of Tanga from ISIMIP3b bias-adjusted climate-model output.
#
#   isimip3b_tanga_climatology.csv      model, day, tas (deg C), pr (mm/day), hurs (%)
#
# Day 1 is 1 January and is model day 0 (as K(t) and the temperature already assumed); Feb 29 is left out.
#
# What it does, in two steps
#   1. buildClimatology()  reads the daily files of the historical ISIMIP3b run in the folder isimip3b_tanga/
#                          (one value per day at the grid cell of Tanga), averages every day of the year over the
#                          years CLIMATOLOGY_YEARS, smooths the cycle, and writes the table above for each climate
#                          model and for the mean of the models ("ensemble"). Run once, after the data are downloaded.
#   2. loadClimate()       reads the table back for one model. Driver.R and KL_gamma_temp.R use this.
#
# Two kinds of input file are understood (both from the ISIMIP repository):
#   *.csv  the "select point" output of the repository: lines "YYYY-MM-DD,value", no header. Small (about 80 KB per
#          decade) and needs no extra package. This is what ISIMIP3b_get_tanga.py downloads.
#   *.nc   NetCDF files (global, or cut out around Tanga); needs the R package ncdf4. The grid cell nearest to
#          TARGET_LON/TARGET_LAT that has data is used.
#
# Why a climatology and not a particular year: a climate model has its own weather. Its 2003 is not the
# real 2003, only the long-term seasonal cycle (the average over many years) is meant to resemble reality.
########################################################################

# 1. SETTINGS (a value already set before this file is sourced is kept) ------------------------

# The five climate models of ISIMIP3b (Lange 2021)
ISIMIP_MODEL_NAMES <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")

if (!exists("ISIMIP_DATA_DIR"))      ISIMIP_DATA_DIR      <- "isimip3b_tanga"                   # folder with the downloaded files
if (!exists("ISIMIP_CLIMATE_FILE"))  ISIMIP_CLIMATE_FILE  <- "isimip3b_tanga_climatology.csv"   # table made by buildClimatology()

# Target point. Tanga city is at 5.07 S, 39.10 E.  For the CSV files the repository has already picked the 0.5 degree cell of the point;
# for NetCDF files the cell containing the point is used, or the nearest cell with data if that one is sea.
if (!exists("TARGET_LON")) TARGET_LON <- 39.10
if (!exists("TARGET_LAT")) TARGET_LAT <- -5.07

# Years averaged into the climatology (the historical run ends in 2014). 1995-2014 = the last 20 years.
if (!exists("CLIMATOLOGY_YEARS")) CLIMATOLOGY_YEARS <- c(1995, 2014)

# Width (days, odd) of the circular moving average applied to the 365-day mean cycle. The average over
# 20 years is still rough for daily rainfall; 1 = no smoothing.
if (!exists("SMOOTH_DAYS")) SMOOTH_DAYS <- 15

ISIMIP_VARIABLES   <- c("tas", "pr", "hurs")     # mean air temperature, precipitation, relative humidity
DAYS_BEFORE_MONTH  <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)   # 365-day year


# 2. READING ONE FILE ---------------------------------------------------------------------------------

isimipCalendar <- function(time, units, calendar = "standard"){
  #. isimipCalendar: year and day of the 365-day year (1 January = 1, Feb 29 = NA) of the time axis of a NetCDF file.
  #.   units are "<days|hours|seconds> since YYYY-MM-DD ..." (ISIMIP: "days since 1601-1-1 00:00:00").
  m <- regmatches(units, regexec("^\\s*([A-Za-z]+)\\s+since\\s+(-?[0-9]+)-([0-9]+)-([0-9]+)", units))[[1]]
  if (length(m) != 5) stop("cannot read the time units: ", units)
  per_day <- c(days = 1, day = 1, hours = 24, hour = 24, minutes = 1440, minute = 1440, seconds = 86400, second = 86400)[tolower(m[2])]
  if (is.na(per_day)) stop("unknown time unit: ", m[2])
  days   <- floor(time/per_day + 1e-6)
  origin <- as.integer(m[3:5])
  calendar <- tolower(calendar)
  if (calendar %in% c("noleap", "365_day", "no_leap")) {
    # every year has 365 days: count the days from the origin directly
    n0   <- origin[1]*365 + DAYS_BEFORE_MONTH[origin[2]] + origin[3] - 1
    n    <- n0 + days
    year <- n %/% 365
    doy  <- n %% 365 + 1
  } else if (calendar %in% c("standard", "gregorian", "proleptic_gregorian", "")) {
    date  <- as.Date(sprintf("%04d-%02d-%02d", origin[1], origin[2], origin[3])) + days
    year  <- as.integer(format(date, "%Y"))
    month <- as.integer(format(date, "%m"))
    day   <- as.integer(format(date, "%d"))
    doy   <- DAYS_BEFORE_MONTH[month] + day
    doy[month == 2 & day == 29] <- NA
  } else stop("calendar not supported: ", calendar)
  data.frame(year = year, doy = doy)
}

convertIsimipUnits <- function(x, variable, units){
  #. convertIsimipUnits: tas to deg C, pr to mm/day, hurs to %, from the units of the file
  u <- tolower(gsub("\\s+", "", units))
  if (variable == "tas") {
    if (u %in% c("k", "kelvin") || (u == "" && median(x, na.rm = TRUE) > 100)) x <- x - 273.15
  } else if (variable == "pr") {
    if (grepl("kg", u) || (u == "" && median(x, na.rm = TRUE) < 0.01)) x <- x*86400    # kg m-2 s-1 = mm/s
  } else if (variable == "hurs") {
    if (u %in% c("1", "fraction") || (u == "" && max(x, na.rm = TRUE) <= 1.0001)) x <- x*100
  }
  x
}

isimipVariable <- function(file){
  #. isimipVariable: tas, pr or hurs, from the file name (".._historical_tas_global_daily_2001_2010.nc")
  b <- basename(file)
  v <- regmatches(b, regexpr(paste0("_(", paste(ISIMIP_VARIABLES, collapse = "|"), ")_"), b))
  if (length(v) == 0) NA_character_ else gsub("_", "", v)
}

cellFromReadme <- function(file, variable){
  #. cellFromReadme: longitude and latitude of the grid cell that the repository used for a "select point" CSV file.
  #.   The README.txt next to the file lists the command (cdo ... -selindexbox,i,i,j,j ...); the ISIMIP grid has
  #.   0.5 degree cells, i counted from 179.75 W and j from 89.75 N. NA if the README is not there.
  cell <- c(lon = NA_real_, lat = NA_real_)
  readme <- file.path(dirname(file), "README.txt")
  yrs <- regmatches(basename(file), regexec("_([0-9]{4})_([0-9]{4})\\.csv$", basename(file)))[[1]]
  if (!file.exists(readme) || length(yrs) != 3) return(cell)
  lines <- grep(paste0("_", variable, "_global_daily_", yrs[2], "_", yrs[3], "\\.nc"), readLines(readme, warn = FALSE), value = TRUE)
  m <- if (length(lines)) regmatches(lines[1], regexec("selindexbox,([0-9]+),([0-9]+),([0-9]+),([0-9]+)", lines[1]))[[1]] else character(0)
  if (length(m) == 5 && m[2] == m[3] && m[4] == m[5])
    cell <- c(lon = -179.75 + (as.numeric(m[2]) - 1)*0.5, lat = 89.75 - (as.numeric(m[4]) - 1)*0.5)
  cell
}

readIsimipCsv <- function(file){
  #. readIsimipCsv: a "select point" CSV of the ISIMIP repository: lines "YYYY-MM-DD,value", no header.
  #.   Units are those of the ISIMIP3b files (tas in K, pr in kg m-2 s-1, hurs in %).
  variable <- isimipVariable(file)
  if (is.na(variable)) stop(basename(file), ": not a tas, pr or hurs file")
  d <- read.csv(file, header = FALSE, col.names = c("date", "value"), stringsAsFactors = FALSE)
  date <- as.Date(d$date)
  if (anyNA(date)) stop(basename(file), ": cannot read the dates")
  month <- as.integer(format(date, "%m")); day <- as.integer(format(date, "%d"))
  doy <- DAYS_BEFORE_MONTH[month] + day; doy[month == 2 & day == 29] <- NA
  value <- suppressWarnings(as.numeric(d$value)); value[abs(value) > 1e19] <- NA
  units <- c(tas = "K", pr = "kg m-2 s-1", hurs = "%")[[variable]]
  list(variable = variable, cell = cellFromReadme(file, variable),
       data = data.frame(year = as.integer(format(date, "%Y")), doy = doy, value = convertIsimipUnits(value, variable, units)))
}

readIsimipCell <- function(file, lon = TARGET_LON, lat = TARGET_LAT){
  #. readIsimipCell: daily values of the grid cell nearest to (lon, lat) that has data, from one ISIMIP NetCDF file.
  #.   Returns list(variable, cell = c(lon, lat), data = data.frame(year, doy, value)).
  if (!requireNamespace("ncdf4", quietly = TRUE)) stop("the package ncdf4 is needed to read NetCDF files")
  nc <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(nc))
  variable <- intersect(names(nc$var), ISIMIP_VARIABLES)
  if (length(variable) != 1) stop(basename(file), ": expected one of ", paste(ISIMIP_VARIABLES, collapse = ", "),
                                  " but found ", paste(names(nc$var), collapse = ", "))
  v     <- nc$var[[variable]]
  dnames <- vapply(v$dim, function(d) d$name, "")          # dimensions in the order of the array returned by ncvar_get
  x     <- ncdf4::ncvar_get(nc, variable, collapse_degen = FALSE)    # keeps dimensions of length one
  x[abs(x) > 1e19] <- NA                                    # ISIMIP fill value 1e20
  units <- ncdf4::ncatt_get(nc, variable, "units")$value
  if (!is.character(units)) units <- ""

  i_time <- grep("^time", dnames); i_lon <- grep("^lon", dnames); i_lat <- grep("^lat", dnames)
  if (length(i_time) != 1) stop(basename(file), ": no time dimension")
  tvals  <- ncdf4::ncvar_get(nc, dnames[i_time])
  tunits <- ncdf4::ncatt_get(nc, dnames[i_time], "units")$value
  tcal   <- ncdf4::ncatt_get(nc, dnames[i_time], "calendar")$value
  cal    <- isimipCalendar(tvals, tunits, if (is.character(tcal)) tcal else "standard")

  if (length(i_lon) == 1 && length(i_lat) == 1) {
    lons <- as.numeric(ncdf4::ncvar_get(nc, dnames[i_lon])); lats <- as.numeric(ncdf4::ncvar_get(nc, dnames[i_lat]))
    x <- aperm(x, c(i_lon, i_lat, i_time))                  # lon, lat, time
    has_data <- apply(x, c(1, 2), function(s) mean(is.finite(s)) > 0.99)
    if (!any(has_data)) stop(basename(file), ": no grid cell with data")
    dist <- outer(((lons - lon)*cos(lat*pi/180))^2, (lats - lat)^2, "+")   # squared distance in degrees
    dist[!has_data] <- Inf
    k <- arrayInd(which.min(dist), dim(dist))
    series <- x[k[1], k[2], ]
    cell   <- c(lon = lons[k[1]], lat = lats[k[2]])
  } else {                                                   # a file that holds a single point
    series <- as.numeric(x)
    cell   <- c(lon = NA_real_, lat = NA_real_)
  }
  list(variable = variable, cell = cell,
       data = data.frame(cal, value = convertIsimipUnits(series, variable, units)))
}

readIsimipSeries <- function(file, lon = TARGET_LON, lat = TARGET_LAT){
  #. readIsimipSeries: one file (.csv or .nc) -> list(variable, cell, data = data.frame(year, doy, value))
  if (grepl("\\.csv$", file, ignore.case = TRUE)) readIsimipCsv(file) else readIsimipCell(file, lon, lat)
}


# 3. FROM THE FILES TO THE 365-DAY CLIMATOLOGY ------------------------------------------------------

isimipFileTable <- function(data_dir = ISIMIP_DATA_DIR){
  #. isimipFileTable: the ISIMIP3b historical files (.csv or .nc) in a folder and its subfolders: model, variable, years.
  #.   Names look like  gfdl-esm4_r1i1p1f1_w5e5_historical_tas_global_daily_2001_2010.nc  (NetCDF) or
  #.   gfdl-esm4_r1i1p1f1_w5e5_historical_tas_lon39.1lat-5.07_daily_2001_2010.csv        (select point)
  files <- list.files(data_dir, pattern = "\\.(nc|csv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("no .csv or .nc files in '", data_dir, "'. Get them with ISIMIP3b_get_tanga.py (see ISIMIP_CHANGES.md).")
  rows <- lapply(files, function(f){
    b <- basename(f)
    model <- ISIMIP_MODEL_NAMES[match(sub("^([A-Za-z0-9-]+)_r[0-9]+i[0-9]+p[0-9]+f[0-9]+_.*$", "\\1", tolower(b)), tolower(ISIMIP_MODEL_NAMES))]
    variable <- isimipVariable(f)
    yrs <- regmatches(b, regexec("_([0-9]{4})_([0-9]{4})\\.(nc|csv)$", b, ignore.case = TRUE))[[1]]
    if (is.na(model) || is.na(variable) || length(yrs) != 4 || !grepl("_historical_", b)) return(NULL)
    data.frame(file = f, model = model, variable = variable, from = as.integer(yrs[2]), to = as.integer(yrs[3]))
  })
  tab <- do.call(rbind, rows)
  if (is.null(tab)) stop("none of the files in '", data_dir, "' is a historical ISIMIP3b file of tas, pr or hurs")
  tab[order(tab$model, tab$variable, tab$from), ]
}

circularSmooth <- function(x, window){
  #. circularSmooth: moving average of odd width `window` over a cycle (day 365 is followed by day 1)
  if (window <= 1) return(x)
  if (window %% 2 == 0) stop("SMOOTH_DAYS must be odd")
  n <- length(x); shift <- -(window %/% 2):(window %/% 2)
  rowMeans(sapply(shift, function(k) x[((seq_len(n) - 1 + k) %% n) + 1]))
}

writeClimatology <- function(cycles, info, years = CLIMATOLOGY_YEARS, lon = TARGET_LON, lat = TARGET_LAT,
                             smooth_days = SMOOTH_DAYS, out_file = ISIMIP_CLIMATE_FILE){
  #. writeClimatology: smooths the 365-day mean cycles, adds the mean of the models ("ensemble") and writes
  #.   out_file (model, day, tas, pr, hurs) and a *_info.csv (cell, years, smoothing, number of files).
  #.   cycles[[model]] = list(tas =, pr =, hurs =), each 365 values (a variable may be missing);
  #.   info[[model]]   = list(cell = c(lon, lat), files = n)
  out <- lapply(names(cycles), function(model){
    cy <- cycles[[model]]
    get <- function(v) if (is.null(cy[[v]])) rep(NA_real_, 365) else circularSmooth(as.numeric(cy[[v]]), smooth_days)
    data.frame(model = model, day = 1:365, tas = get("tas"), pr = get("pr"), hurs = get("hurs"))
  })
  names(out) <- names(cycles)
  ens <- data.frame(model = "ensemble", day = 1:365,
                    sapply(c("tas", "pr", "hurs"), function(v){
                      m <- sapply(out, function(o) o[[v]]); if (all(is.na(m))) rep(NA_real_, 365) else rowMeans(m, na.rm = TRUE) }))
  clim <- rbind(do.call(rbind, out), ens); rownames(clim) <- NULL
  write.csv(clim, out_file, row.names = FALSE)
  info_df <- do.call(rbind, lapply(names(info), function(model)
    data.frame(model = model, cell_lon = info[[model]]$cell[["lon"]], cell_lat = info[[model]]$cell[["lat"]], target_lon = lon, target_lat = lat,
               years_from = years[1], years_to = years[2], smooth_days = smooth_days, files = info[[model]]$files)))
  write.csv(info_df, sub("\\.csv$", "_info.csv", out_file), row.names = FALSE)
  message("Saved ", out_file, " (", length(unique(clim$model)), " series: ", paste(unique(clim$model), collapse = ", "), ")")
  invisible(clim)
}

buildClimatology <- function(data_dir = ISIMIP_DATA_DIR, years = CLIMATOLOGY_YEARS, lon = TARGET_LON, lat = TARGET_LAT,
                             smooth_days = SMOOTH_DAYS, out_file = ISIMIP_CLIMATE_FILE){
  #. buildClimatology: reads the files, makes the 365-day mean cycle of tas, pr and hurs for each climate model
  #.   and for their mean ("ensemble"), and writes out_file (and a *_info.csv with the cell and years used).
  tab  <- isimipFileTable(data_dir)
  need <- seq(years[1], years[2])
  cycles <- list(); info <- list()
  for (model in unique(tab$model)) {
    cy <- list(); cells <- list(); nfiles <- 0
    for (variable in ISIMIP_VARIABLES) {
      sel <- tab[tab$model == model & tab$variable == variable & tab$to >= years[1] & tab$from <= years[2], ]
      if (nrow(sel) == 0) { message(model, ": no ", variable, " files for ", years[1], "-", years[2], " (left empty)"); next }
      parts <- lapply(sel$file, readIsimipSeries, lon = lon, lat = lat)
      d <- do.call(rbind, lapply(parts, function(p) p$data))
      d <- d[!is.na(d$doy) & d$year >= years[1] & d$year <= years[2], ]
      counts <- table(factor(d$year, levels = need))
      if (any(counts != 365)) stop(model, " ", variable, ": the files do not hold every day of ", years[1], "-", years[2],
                                   " (years with missing or extra days: ", paste(names(counts)[counts != 365], collapse = ", "), ")")
      if (any(!is.finite(d$value))) stop(model, " ", variable, ": missing values in the cell used")
      cy[[variable]]    <- as.numeric(tapply(d$value, factor(d$doy, levels = 1:365), mean))    # mean over the years, per day
      cells[[variable]] <- parts[[1]]$cell
      nfiles <- nfiles + nrow(sel)
    }
    cycles[[model]] <- cy
    info[[model]]   <- list(cell = cells[[1]], files = nfiles)
  }
  writeClimatology(cycles, info, years, lon, lat, smooth_days, out_file)
}


# 4. USING THE CLIMATOLOGY --------------------------------------------------------------------------

climateTag <- function(model) gsub("[^a-z0-9]+", "-", tolower(model))   # for file names, e.g. "mpi-esm1-2-hr"

loadClimate <- function(model = "ensemble", file = ISIMIP_CLIMATE_FILE){
  #. loadClimate: the 365-day climate of one model (or "ensemble") as a data frame: day, tas (deg C), pr (mm/day), hurs (%)
  if (!file.exists(file)) stop("'", file, "' does not exist. Download the ISIMIP3b files (ISIMIP3b_get_tanga.py), then run ",
                               "source('ISIMIP3b_climate.R'); buildClimatology()   (see ISIMIP_CHANGES.md)")
  all <- read.csv(file, stringsAsFactors = FALSE)
  if (!model %in% all$model) stop("model '", model, "' is not in ", file, ". Available: ", paste(unique(all$model), collapse = ", "))
  clim <- all[all$model == model, c("day", "tas", "pr", "hurs")]
  clim <- clim[order(clim$day), ]; rownames(clim) <- NULL
  if (nrow(clim) != 365 || any(is.na(clim$tas)) || any(is.na(clim$pr))) stop(model, ": needs 365 days of tas and pr")
  if (min(clim$tas) < -10 || max(clim$tas) > 50) warning(model, ": temperatures outside -10 to 50 deg C, check the units")
  clim
}

describeClimate <- function(clim, label = ""){
  #. describeClimate: a few numbers for the log
  cat(sprintf("Climate %s: T mean %.1f C (range %.1f to %.1f), rainfall %.0f mm/year (peak day %d), humidity %s\n",
              label, mean(clim$tas), min(clim$tas), max(clim$tas), sum(clim$pr), which.max(clim$pr),
              if (all(is.na(clim$hurs))) "not available" else sprintf("%.0f%% (min %.0f%%)", mean(clim$hurs), min(clim$hurs))))
}


# 5. CHECK PLOT: THE MODELS AGAINST CRU ------------------------------------------------------------

plotClimateModels <- function(file = ISIMIP_CLIMATE_FILE, pdf_file = "ISIMIP3b_vs_CRU.pdf", cru_dir = "cru_reference"){
  #. plotClimateModels: annual cycle of temperature, rainfall and humidity of each model and of the ensemble,
  #.   against the CRU monthly values and rain.csv of the earlier version (if cru_dir still holds them).
  all <- read.csv(file, stringsAsFactors = FALSE)
  models <- setdiff(unique(all$model), "ensemble"); cols <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#AA3377")[seq_along(models)]
  mid <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349); ndays <- diff(c(DAYS_BEFORE_MONTH, 365))
  f_t <- file.path(cru_dir, "monthly_tmp_2019_2020.rds"); f_p <- file.path(cru_dir, "monthly_pre_2019_2020.rds"); f_r <- file.path(cru_dir, "rain.csv")
  cru_t <- if (file.exists(f_t)) rowMeans(matrix(readRDS(f_t), nrow = 12)) else NULL                 # deg C, mean of 2019 and 2020
  cru_p <- if (file.exists(f_p)) rowMeans(matrix(readRDS(f_p), nrow = 12))/ndays else NULL           # mm/day
  cru_r <- if (file.exists(f_r)) circularSmooth(pmax(read.csv(f_r)[[1]][1:365], 0), 15) else NULL     # mm/day, 15-day mean
  panel <- function(variable, ylab, main, extra = NULL){
    plot(NA, xlim = c(1, 365), ylim = range(c(all[[variable]], extra), na.rm = TRUE), xlab = "Day of year", ylab = ylab, main = main)
    for (i in seq_along(models)) lines(1:365, all[[variable]][all$model == models[i]], col = cols[i], lwd = 1)
    lines(1:365, all[[variable]][all$model == "ensemble"], col = "black", lwd = 2.5)
  }
  pdf(pdf_file, width = 8, height = 9); on.exit(dev.off())
  par(mfrow = c(3, 1), mar = c(4, 4.5, 2, 1))
  panel("tas", "Temperature (deg C)", "Temperature: ISIMIP3b climatology", c(cru_t, max(all$tas, na.rm = TRUE) + 0.9))   # headroom for the legend
  if (!is.null(cru_t)) points(mid, cru_t, pch = 19, col = "grey40")
  legend("top", ncol = 3, bty = "n", cex = 0.7, legend = c(models, "ensemble", if (!is.null(cru_t)) "CRU 2019-2020 (monthly)"),
         col = c(cols, "black", if (!is.null(cru_t)) "grey40"), lwd = c(rep(1, length(models)), 2.5, if (!is.null(cru_t)) NA),
         pch = c(rep(NA, length(models) + 1), if (!is.null(cru_t)) 19))
  panel("pr", "Rainfall (mm/day)", "Rainfall: ISIMIP3b climatology", c(cru_p, cru_r))
  if (!is.null(cru_p)) points(mid, cru_p, pch = 19, col = "grey40")
  if (!is.null(cru_r)) lines(1:365, cru_r, col = "grey60", lty = 2)
  legend("topright", bty = "n", cex = 0.7, legend = c("CRU monthly (2019-2020 mean)", "rain.csv (15-day mean)"), col = c("grey40", "grey60"), pch = c(19, NA), lty = c(NA, 2))
  if (!all(is.na(all$hurs))) panel("hurs", "Relative humidity (%)", "Humidity: ISIMIP3b climatology")
  invisible(NULL)
}
