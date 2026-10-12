########################################################################
# Plots.R
# All plotting functions of the model.
#   1. Standard ggplot figures          (used by Driver.R and by the Shiny downloads)
#   2. Extra / diagnostic ggplot figures
#   3. Manuscript figures
#   4. Interactive plotly figures        (used by the Shiny app)
# (Replaces the plotting half of ODEAuxiliaryFunctions.R and the old Plots.R,
#  which defined several of these functions twice.)
#
# Note: in every ggplot below, `color = State` in the main aes() only supplies the
# legend title "State"; the lines get their data and colour from the geom_line() layers.
########################################################################

######################################################################################
# Helpers shared by the bar charts of the mosquito life stages
######################################################################################
DEMOGRAPHICS_STAGE_CODES <- c("EL","LL","PL","SV","EV","IV")

# Total number of mosquitoes in each life stage over the whole simulation
stageTotals = function(IVM_traj){
  unname(vapply(DEMOGRAPHICS_STAGE_CODES, function(stage) sum(IVM_traj[stage]), numeric(1)))
}

######################################################################################
# 1. STANDARD GGPLOT FIGURES
######################################################################################
plotTrajectory = function(IVM_traj, log_y = FALSE){
  #. plotTrajectory: Plots the evolution of the dynamics of the system (susceptible, exposed, infected mosquitoes)
  #  log_y = TRUE gives the log10 y-axis used by the Shiny app (see plotTrajectoryShiny)
  ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    theme_grey(base_size = 18, base_family = "") +
    geom_line(aes(y = IV, col = "Infected"), size = 1.75) +
    geom_line(aes(y = SV, col = "Susceptible"), size = 1.75) +
    geom_line(aes(y = EV, col = "Exposed"), size = 1.75) +
    (if (log_y) scale_y_log10() else ylim(0, NA)) +
    labs(x = "Time (days)", y = "Number of mosquitoes")
}

plotTrajectoryHumans = function(IVM_traj){
  #. plotTrajectoryHumans: Plots the susceptible and infected humans
  ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    theme_grey(base_size = 18, base_family = "") +
    geom_line(aes(y = IH, col = "Infected"), size = 1.75) +
    geom_line(aes(y = SH, col = "Susceptible"), size = 1.75) +
    ylim(0, NA) +
    labs(x = "Time (days)", y = "Number of humans")
}

# One output variable of the trajectory as a single coloured line
plotSingleSeries = function(IVM_traj, column, line_colour, y_label, title = NULL, base_size = 18){
  ggplot(IVM_traj, aes(x = time, y = EIR, color = State)) +
    theme_grey(base_size = base_size, base_family = "") +
    geom_line(aes(y = .data[[column]], col = column), size = 1.75, colour = line_colour) +
    ylim(0, NA) +
    labs(title = title, x = "Time (days)", y = y_label)
}

plotEIR = function(IVM_traj, line_colour = 'blue', y_label = "EIR", title = "An. Gambiae Baseline EIR = 100", base_size = 18){
  #. plotEIR: Plots the entomological inoculation rate (EIR)
  #  The title is fixed text: change it (title = ...) to match the species / EIR you simulate.
  plotSingleSeries(IVM_traj, "EIR", line_colour, y_label, title, base_size)
}

plotVC = function(IVM_traj){
  #. plotVC: Plots the vectorial capacity
  plotSingleSeries(IVM_traj, "VC", 'purple', "Vectorial Capacity (VC)")
}

plotR0 = function(IVM_traj){
  #. plotR0: Plots the basic reproduction rate
  plotSingleSeries(IVM_traj, "R0", "blue", "Basic Reproduction Rate (R0)")
}

barChartMosquitoDemographics = function(IVM_traj){
  #. barChartMosquitoDemographics: Bar chart (base graphics) with the amount of mosquitoes in each life stage
  barplot(stageTotals(IVM_traj), main='', xlab='Time (days)',ylab='Number of Mosquitoes',names.arg=DEMOGRAPHICS_STAGE_CODES)
}

######################################################################################
# Versions of the standard figures used by the Shiny app (downloads / PDF report)
######################################################################################
plotTrajectoryShiny = function(IVM_traj){
  #. plotTrajectoryShiny: plotTrajectory() with a log10 y-axis
  plotTrajectory(IVM_traj, log_y = TRUE)
}

plotEIRShiny = function(IVM_traj){
  #. plotEIRShiny: plotEIR() in magenta, with the full y-axis label, no title and the default (smaller, size 11) text
  plotEIR(IVM_traj, line_colour = 'magenta', y_label = "Entomological Inoculation Rate (EIR)", title = NULL, base_size = 11)
}

barChartMosquitoDemographicsShiny = function(IVM_traj, stage_labels = DEMOGRAPHICS_STAGE_CODES){
  #. barChartMosquitoDemographicsShiny: ggplot bar chart (log10 y-axis) with the amount of mosquitoes in each life stage
  traj_dat <- data.frame(stage = stage_labels, sum = stageTotals(IVM_traj))
  traj_dat$stage <- factor(traj_dat$stage,levels=traj_dat$stage)
  ggplot(data=traj_dat,aes(stage,sum,fill=stage)) +
    geom_bar(stat="identity",colour="black") +
    theme_grey(base_size = 18, base_family = "") +
    scale_y_log10() +
    guides(fill=FALSE) +
    labs(x="Mosquito Stage",y="Number of Mosquitoes")
}

######################################################################################
# 2. EXTRA / DIAGNOSTIC GGPLOT FIGURES
######################################################################################
plotTrajectoryMalesMosq = function(IVM_traj){
  #. plotTrajectoryMalesMosq: Plots the male mosquitoes
  ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    geom_line(aes(y = Mm, col = "Mm"), size = 1.2) +
    labs(x = "Time (days)", y = "Number of Male Mosquitoes")
}

plotTrajectoryALL = function(IVM_traj){
  #. plotTrajectoryALL: Plots all mosquito compartments (log10 y-axis)
  ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    geom_line(aes(y = EL, col = "Early Instar"), size = 1.2) +
    geom_line(aes(y = LL, col = "Late Instar"), size = 1.2) +
    geom_line(aes(y = PL, col = "Pupae"), size = 1.2) +
    geom_line(aes(y = IV, col = "Infected"), size = 1.2) +
    geom_line(aes(y = SV, col = "Susceptible"), size = 1.2) +
    geom_line(aes(y = EV, col = "Exposed"), size = 1.2) +
    geom_line(aes(y = EL_LAR, col = "EL_LAR"), size = 1.2) +
    geom_line(aes(y = EL_BIO, col = "EL_BIO"), size = 1.2) +
    geom_line(aes(y = EL_LAR_BIO, col = "EL_LAR_BIO"), size = 1.2) +
    geom_line(aes(y = LL_LAR, col = "LL_LAR"), size = 1.2) +
    geom_line(aes(y = LL_BIO, col = "LL_BIO"), size = 1.2) +
    geom_line(aes(y = LL_LAR_BIO, col = "LL_LAR_BIO"), size = 1.2) +
    labs(x = "Time (days)", y = "Number of mosquitoes") +
    scale_y_log10() +
    theme_bw()
}

plotTrajectoryDEBUG <- function(traj,log=FALSE){
  #. plotTrajectoryDEBUG: Used only to debug (every column of the trajectory, optionally on a log axis)
  require(reshape2)
  IVM_melt <- melt(traj,id.vars="time")
  p <- ggplot(data=IVM_melt,aes(x=time,y=value,group=variable,colour=variable,linetype=variable)) +
    geom_line(size=1.15) +
    theme_bw()
  if(log){
    p <- p + scale_y_log10()
  }
  p
}

plotEIR_1 = function(IVM_traj,EIR_1){
  #. plotEIR_1: Plots the EIR vector `EIR_1` against the time of IVM_traj
  ggplot(IVM_traj, aes(x = time, y = EIR, color = State)) +
    geom_line(aes(y = EIR_1, col = "EIR"), size = 1.75) +
    labs(x = "Time (days)", y = " EIR")
}

plotEIR_VC_R0 = function(IVM_traj){
  #. plotEIR_VC_R0: Plots mosquito numbers, EIR, VC and R0 dynamics of the system in one panel

  # Mosquito density
  p1 <- ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    geom_line(aes(y = SV+EV+IV, col = "NV"), size = 1.2) +
    labs(x = "Time (days)", y = "Number of total mosquitoes")

  # Infectious vectors
  p2 <- ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    geom_line(aes(y = IV, col = "IV"), size = 1.2) +
    labs(x = "Time (days)", y = "Number of infectious mosquitoes")

  # Total and infectious vectors
  p3 <- ggplot(IVM_traj, aes(x = time, y = IVM_traj, color = State)) +
    geom_line(aes(y = SV+EV+IV, col = "NV"), size = 1.2) +
    geom_line(aes(y = IV, col = "IV"), size = 1.2) +
    labs(x = "Time (days)", y = "Number of infections, mosquitoes")

  p4 <- ggplot(IVM_traj, aes(x = time, y = EIR, color = State)) +
    geom_line(aes(y = IVM_traj[["EIR"]], col = "EIR"), size = 1.2) +
    labs(x = "Time (days)", y = " EIR")

  p5 <- ggplot(IVM_traj, aes(x = time, y = EIR, color = State)) +
    geom_line(aes(y = IVM_traj[["R0"]], col = "R0"), size = 1.2) +
    labs(x = "Time (days)", y = " R0")

  p6 <- ggplot(IVM_traj, aes(x = time, y = EIR, color = State)) +
    geom_line(aes(y = IVM_traj[["VC"]], col = "VC"), size = 1.2) +
    labs(x = "Time (days)", y = " VC")

  multiplot(p1, p2, p3, p4,p5,p6, cols=2)
}

barChartMosquitoDemographics_slwu <- function(IVM_traj){
  #. barChartMosquitoDemographics_slwu: Demographics bar chart with the full life-stage names
  barChartMosquitoDemographicsShiny(IVM_traj,
    stage_labels = c("Early Instar","Late Instar","Pupae","Susceptible","Exposed","Infected"))
}

######################################################################################
# 3. MANUSCRIPT FIGURES
#    (use multiplot() from multiplot.R and stringr::str_wrap)
######################################################################################
plot_MS_Figures = function(IVM_traj_Control,IVM_traj_LLIN_50,IVM_traj_LLIN_80,IVM_traj_LLIN_90,theta,INITIAL_MODELRUNTIME_VALUE){
  #. plot_MS_Figures - produce some figures for the manuscript

  # Mosquito population
  NV_Control = IVM_traj_Control[["SV"]]+IVM_traj_Control[["EV"]]+IVM_traj_Control[["IV"]]
  NV_LLIN_80 = IVM_traj_LLIN_80[["SV"]]+IVM_traj_LLIN_80[["EV"]]+IVM_traj_LLIN_80[["IV"]]
  NV_LLIN_50 = IVM_traj_LLIN_50[["SV"]]+IVM_traj_LLIN_50[["EV"]]+IVM_traj_LLIN_50[["IV"]]
  NV_LLIN_90 = IVM_traj_LLIN_90[["SV"]]+IVM_traj_LLIN_90[["EV"]]+IVM_traj_LLIN_90[["IV"]]

  # Just mosquito density
  p <- ggplot(IVM_traj_Control, aes(x = time, y = IVM_traj_Control, color = State)) +
    geom_line(aes(y = NV_Control, col = "Mosquito-Control"), size = 1.2) +
    geom_line(aes(y = NV_LLIN_50, col = "Mosquito-LLIN-50"), size = 1.2) +
    geom_line(aes(y = NV_LLIN_80, col = "Mosquito-LLIN-80"), size = 1.2) +
    geom_line(aes(y = NV_LLIN_90, col = "Mosquito-LLIN-80"), size = 1.2) +   # NOTE: label says -80 as in the original
    labs(title = "An. Arabiensis", x = "Time (days)", y = "Number of mosquitoes")
  p <- p + expand_limits(x = 0, y = 0)
  multiplot(p, cols=1)

  # EIR
  p2 <- ggplot(IVM_traj_Control, aes(x = 1:length(IVM_traj_Control[,1]), y = EIR, color = key)) +
    geom_line(aes(y = IVM_traj_Control[["EIR"]] , col = "EIR - Control "), size = 1.2) +
    geom_line(aes(y = IVM_traj_LLIN_50[["EIR"]], col = "EIR - LLIN - 50 "), size = 1.2) +
    geom_line(aes(y = IVM_traj_LLIN_80[["EIR"]], col = "EIR - LLIN - 80"), size = 1.2) +
    geom_line(aes(y = IVM_traj_LLIN_90[["EIR"]], col = "EIR - LLIN - 90"), size = 1.2) +
    labs(title = "An. Arabiensis", x = "Time (days)", y = "EIR")

  multiplot(p2, cols=1)
}

# Bar chart of the final EIR for a set of intervention combinations
barChartFinalEIR = function(trajectories, stage_labels, title){
  traj_dat <- data.frame(stage = stage_labels,
                         value = vapply(trajectories, function(traj) tail(traj[["EIR"]],1), numeric(1), USE.NAMES = FALSE))
  traj_dat$stage <- factor(traj_dat$stage,levels=traj_dat$stage)

  traj_dat$newx = stringr::str_wrap(traj_dat$stage, width = 10)

  ggplot(data=traj_dat,aes(newx,value,fill=stage)) +
    geom_bar(stat="identity",colour="black") +
    theme_grey(base_size = 18, base_family = "") +
    guides(fill=FALSE) +
    labs(title = title, x="Tools",y="EIR")
}

barChart_MS_Other <- function (IVM_traj_Control,IVM_traj_ITN_LAR_5070,IVM_traj_ITN_LAR_ATSB_508050,IVM_traj_ITN_LAR_ATSB_ECS_50805050 )
{
  barChartFinalEIR(
    trajectories = list(IVM_traj_ITN_LAR_5070,IVM_traj_ITN_LAR_ATSB_508050,IVM_traj_ITN_LAR_ATSB_ECS_50805050),
    stage_labels = c("ITN 50%, LAR 80%","ITN 50%, LAR 80%, ATSB 50%","ITN 50%, LAR 80, ATSB 50%, ECT 50%"),
    title = "An. Funestus  : Baseline EIR = 100")
}

barChart_MS <- function(IVM_traj_Control,IVM_traj_ITN_HOU_5050,IVM_traj_ITN_PPM_5050,IVM_traj_ITN_ECS_5050,
                        IVM_traj_ITN_LAR_ECS_505050,IVM_traj_ITN_LAR_5050,IVM_traj_ITN_PPM_LAR_505050
                        ,IVM_traj_LLIN_ECT_PPM_505050,IVM_traj_ITN_ATSB_5050,IVM_traj_ITN_ATSB_LAR_505050 ){
  barChartFinalEIR(
    trajectories = list(IVM_traj_ITN_HOU_5050,IVM_traj_ITN_PPM_5050,IVM_traj_ITN_ECS_5050,
                        IVM_traj_ITN_LAR_ECS_505050,IVM_traj_ITN_LAR_5050,
                        IVM_traj_ITN_PPM_LAR_505050, IVM_traj_LLIN_ECT_PPM_505050,
                        IVM_traj_ITN_ATSB_5050,IVM_traj_ITN_ATSB_LAR_505050),
    stage_labels = c("ITN 80%, HOU 80%","ITN 80%, PPM 80%","ITN 80%, ECT 80%", "ITN 80%, LAR 80%, ECS 50%","ITN 80%, LAR 80%",
                     "ITN 80%, PPM 80% LAR 80%","ITN 80%, ECT 80%, PPM 80%","ITN 80%, ATSB 80%","ITN 80%, ATSB 80%, LAR 80%"),
    title = "An.Arabiensis  : Baseline EIR = 10")
}

######################################################################################
# 4. INTERACTIVE (PLOTLY) FIGURES FOR THE SHINY APP
######################################################################################
LINE_WIDTH = 4
LINE_STYLE = list(width=LINE_WIDTH)
FONT_LABEL <- list(
  size = 15,
  color = "grey"
)
FONT_TICKS <- list(
  size = 10,
  color = "grey"
)

# Axis styles for plotly (linear, bar chart, logarithmic)
AxisStyle = function(label){
  list(
    title = label,
    tickwidth = 2,
    titlefont = FONT_LABEL,
    tickfont = FONT_TICKS,
    tickcolor = toRGB("grey"),
    rangemode = "tozero",
    linewidth = 4,
    linecolor = "grey",
    mirror = "ticks"
  )
}
AxisStyleBar = function(label){
  # AxisStyle with smaller black tick labels and without rangemode
  style <- AxisStyle(label)
  style$tickfont  <- list(size = 8, color = "black")
  style$rangemode <- NULL
  style
}
AxisStyleLog = function(label){
  # AxisStyle with a logarithmic axis (type inserted right after rangemode)
  append(AxisStyle(label), list(type = "log"), after = 6)
}
AxisStyleLogBar = AxisStyleLog   # identical to AxisStyleLog

plotTrajectoryPlotLy = function(IVM_traj){
  #. plotTrajectoryPlotLy: Plots the evolution of the dynamics of the system
  # NOTE (unchanged): SV is plotted twice, as "Susceptible" and as "Sus".
  plot_ly(x = IVM_traj$time, y = IVM_traj$SV, line = LINE_STYLE, name = "Susceptible") %>%
    add_trace(y = IVM_traj$IV, line = LINE_STYLE, name = "Infected") %>%
    add_trace(y = IVM_traj$SV, line = LINE_STYLE, name = "Sus") %>%
    add_trace(y = IVM_traj$EV, line = LINE_STYLE, name = "Exposed") %>%
    layout(xaxis = AxisStyle("Time"), yaxis = AxisStyleLog("Number of Mosquitoes"))
}
plotTrajectoryHumansPlotLy = function(IVM_traj){
  #. plotTrajectoryHumansPlotLy: Plots the susceptible and infected humans
  plot_ly(x = IVM_traj$time, y = IVM_traj$SH, line = LINE_STYLE, name = "Susceptible") %>%
    add_trace(y = IVM_traj$IH, line = LINE_STYLE, name = "Infected") %>%
    layout(xaxis = AxisStyle("Time"), yaxis = AxisStyle("Number of Humans"))
}
plotEIRPlotLy = function(IVM_traj){
  #. plotEIRPlotLy: Plots the EIR
  plot_ly(x = IVM_traj$time, y = IVM_traj$EIR, line = list(width=LINE_WIDTH,color="magenta"), name = "EIR") %>%
    layout(xaxis = AxisStyle("Time"), yaxis = AxisStyle("Entomological Inoculation Rate (EIR)"))
}
plotVCPlotLy = function(IVM_traj){
  #. plotVCPlotLy: Plots the vectorial capacity
  plot_ly(x = IVM_traj$time, y = IVM_traj$VC, line = list(width=LINE_WIDTH,color="red"), name = "VC", color="cyan") %>%
    layout(xaxis = AxisStyle("Time"), yaxis = AxisStyle("Vectorial Capacity (VC)"))
}
plotR0PlotLy = function(IVM_traj){
  #. plotR0PlotLy: Plots the basic reproduction rate
  plot_ly(x = IVM_traj$time, y = IVM_traj$R0, line = list(width=LINE_WIDTH,color="violet"), name = "R0") %>%
    layout(xaxis = AxisStyle("Time"), yaxis = AxisStyle("Basic Reproduction Rate (R0)"))
}
barChartMosquitoDemographicsPlotLy = function(IVM_traj){
  #. barChartMosquitoDemographicsPlotLy: Bar chart with the amount of mosquitoes in each life stage
  stage_labels <- c("Early Instar","Late Instar","Pupae","Susceptible Adult","Exposed Adult","Infected Adult")
  traj_dat <- data.frame(stage = stage_labels, sum = stageTotals(IVM_traj))
  traj_dat$stage <- factor(traj_dat$stage,levels=traj_dat$stage)
  plot_ly(x = stage_labels, y = traj_dat$sum, type = "bar",
          name = as.list(DEMOGRAPHICS_STAGE_CODES), color = as.list(DEMOGRAPHICS_STAGE_CODES)) %>%
    layout(xaxis = AxisStyleBar(""), yaxis = AxisStyleLogBar("Number of Mosquitoes"), showlegend = FALSE)
}
