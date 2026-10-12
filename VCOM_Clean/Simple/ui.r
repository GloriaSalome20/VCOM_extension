#------------------------------------------------------------------------------#
################################################################################
## Malaria vector ODE model GUI                                               ##
## Hector M. Sanchez C. (sanchez.hmsc@itesm.mx)                               ##
## 02/May/2016                                                                ##
################################################################################
#------------------------------------------------------------------------------#
# NOTE: the <<- assignments below are intentional. Shiny evaluates ui.r in its own
# environment, so the objects that server.r needs (theta, initState, IVM_traj, ...)
# have to be written to the global environment.

################################################################################
# LOAD LIBRARIES AND FILES
################################################################################
library(shiny)
library(deSolve)
library(ggplot2)
library(shinyjs)
#library(shinythemes)
library(plotly)
library(BH)
source("ODEModel.R")
source("ODEMosquitoParameters.R")
source("Plots.R")
source("ODEControlMeasuresParameters.R")
source("ODETransmissionParameters.R")
source("ODEInterventions.R")
source("ODEModelOutput.R")
source("FeedingCycle.R")
source("SeasonalSetup.R")      # climate, seasonal K(t), temperature-dependent parameters, burn-in (runVCOM())

################################################################################
# GLOBAL GUI PARAMETERS ########################################################
importedFile = NULL
BOXES_WIDTH <<- "75px"
COVERAGE_STEP_SIZE = .025
COVERAGE_LABELS_SIZE = 3
COVERAGE_BAR_SIZE = 5
COVERAGE_INIT_SIZE = 4
INITIAL_TIME_VALUE = 365
MAX_SIMULATION_DAYS = 1825     # 5 years (the climate repeats every year)
INITIAL_ITN_COVERAGE = .0
INITIAL_EIR = 100
INITIAL_ITN_TIME = 50
# Seasonality (see SeasonalSetup.R and ISIMIP_CHANGES.md): the defaults are those of Driver.R
INITIAL_CLIMATE_MODEL = "ensemble"                   # ISIMIP3b climate model of Tanga: "ensemble" or one of the five models
INITIAL_SEASONAL_K = TRUE                            # seasonal larval carrying capacity K(t) (FALSE: constant K)
INITIAL_TEMPERATURE_GROUPS = TEMPERATURE_GROUPS      # temperature-dependent parameters that are switched on (character(0): none)
INITIAL_BURN_IN_YEARS = 1                            # years run without interventions before day 0
MAX_BURN_IN_YEARS = 2
ENTER_DOWN_RUN <<- '
  $(document).ready(function(){
    $("body").keydown(function(e){
      if(e.which === 13){$("#buttonRun").click();}
    });
  });'
WORKING_MESSAGE_STYLE <<- "#loadmessage {
    position: fixed;top: 50%;left: 0px;line-height: 80px;height: 100px;
    width: 100%;padding: 5px 0px 5px 0px;text-align: center;font-weight: bold;
    font-size: 200%;color: #FFFFFF;background-color: rgba(100, 200, 255, .75);z-index: 105;
  }"

###################################################################################
# PRIME SYSTEM  ###################################################################
# The model is run once with the defaults of the sidebar, so that the plots are there when the page opens.
# runVCOM() (SeasonalSetup.R) does what server.r does when "Run" is clicked and leaves MOSQUITO_PARAMETERS,
# TRANSMISSION_PARAMETERS, INTERVENTION_PARAMETERS, theta, initState, IVM_traj and RUN_SUMMARY in the global environment.
runVCOM(species = "GAM", EIR = INITIAL_EIR, days = INITIAL_TIME_VALUE,
        coverage = c(ITN = INITIAL_ITN_COVERAGE), start = c(ITN = INITIAL_ITN_TIME),
        climateModel = INITIAL_CLIMATE_MODEL, seasonalK = INITIAL_SEASONAL_K,
        temperatureGroups = INITIAL_TEMPERATURE_GROUPS, burnInYears = INITIAL_BURN_IN_YEARS)
#initialParametersValues <<- importCSVParametersFromDirectory("SetupTemplates/SETUP_MosquitoLifeCycleParameters.csv")
#theta <<- parseImportedCSVParameters(initialParametersValues)
timing <<- seq(0,INITIAL_TIME_VALUE,by=1)

###################################################################################
# One row of the interventions table: label | coverage slider | start-day box
interventionRow = function(label, coverageId, timeId, coverage = 0, time = 0){
  fluidRow(
    column(COVERAGE_LABELS_SIZE,h5(label,align="left")),
    column(COVERAGE_BAR_SIZE,sliderInput(coverageId,NULL,min=0,max=1,value=coverage,step=COVERAGE_STEP_SIZE)),
    column(COVERAGE_INIT_SIZE,numericInput(timeId,NULL,value=time,min=0,max=MAX_SIMULATION_DAYS))
  )
}

###################################################################################
# USER INTERFACE ##################################################################
shinyUI(
  fluidPage(theme = "cerulean.css",
            titlePanel(h1("VCOM: Simple",align="center"),windowTitle="VCOM: Simple"),
            titlePanel(h4("Vector Control Optimization Model",align="center")),
            titlePanel(h4(tags$a(class="btn btn-default", href="http://chipdelmal.github.io/VCOM/", "Home", style="display: block; width: 100%;"),align="center")),
            ###############################################################################
            useShinyjs(),
            tags$script(ENTER_DOWN_RUN),
            tags$head(tags$style(type="text/css",WORKING_MESSAGE_STYLE)),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",tags$div("Working...",id="loadmessage")),
            ###############################################################################
            sidebarLayout(
              sidebarPanel(
                #####################################################################
                titlePanel(h1("Instructions",align="left")),
                helpText("(1) Select the mosquito species."),
                helpText("(2) Select the EIR (entomological inoculation rate) level."),
                helpText("(3) Select the climate model and which parts of the model follow the season."),
                helpText("(4) Select the number of days to simulate and the burn-in."),
                helpText("(5) Run the model! (click the button or hit 'ENTER')"),
                helpText("(6 optional) Download a PDF report of the plots."),
                helpText("(7) Setup the desired interventions and repeat step 5 as required."),
                helpText("(8) Eliminate malaria!"),
                fluidRow(h3("1. Mosquito Selection")),
                radioButtons("radioSpecies",label=NULL,
                             choices=list(
                               "An. gambiae"="GAM",
                               "An. arabiensis"="ARA",
                               "An. funestus"="FUN"
                             ),selected="GAM"),
                #####################################################################
                fluidRow(h3("2. EIR Level Selection")),
                sliderInput("radioEIR",label=NULL,min=0.1,max=1000,value=100,step=10),
                #####################################################################
                fluidRow(h3("3. Climate and Seasonality")),
                helpText("Rainfall, temperature and humidity of Tanga from ISIMIP3b climate models (average year 1995-2014).
                          Day 0 is 1 January and the year repeats."),
                selectInput("climateModel","Climate model:",choices=SEASONAL_CLIMATE_CHOICES,selected=INITIAL_CLIMATE_MODEL),
                checkboxInput("seasonalK","Seasonal larval carrying capacity K(t) (rainfall and temperature)",value=INITIAL_SEASONAL_K),
                checkboxGroupInput("temperatureGroups","Temperature-dependent parameters:",
                                   choices=SEASONAL_TEMPERATURE_CHOICES,selected=INITIAL_TEMPERATURE_GROUPS),
                #####################################################################
                fluidRow(h3("4. Simulation Time")),
                sliderInput("sliderTime","Days to Simulate:",min=1,max=MAX_SIMULATION_DAYS,value=INITIAL_TIME_VALUE),
                sliderInput("burnInYears","Burn-in (years before day 0, without interventions):",
                            min=0,max=MAX_BURN_IN_YEARS,value=INITIAL_BURN_IN_YEARS,step=1),
                helpText("The burn-in runs the model for these years first and drops them, so that day 0 starts in the model's own
                          seasonal state and not at the equilibrium of a constant environment (which is far too high on 1 January).
                          Interventions start on the days set in step 7, counted from day 0. 0 = no burn-in."),
                #####################################################################
                fluidRow(h3("5. Run Model")),
                actionButton("buttonRun","Run",width="100%"),
                #####################################################################
                fluidRow(h3("6. Download PDF and CSV Trace")),
                downloadButton("downloadPlots", 'Download PDF Plots Report'),
                downloadButton("downloadCSVTrace", 'Download CSV Trace'),

                #####################################################################
                fluidRow(h3("7. Interventions")),
                fluidRow(helpText(
                  "Key: Ovitraps (OVI), Larviciding (LAR), Biological Control (BIO), Source Reduction (SRE)
                Insecticide-Treated Nets (ITN), Indoor Residual Spraying (IRS), Endocticide (ECS), Endocticide (ECT), House Modifications (HOU),
                Odor-Baited Traps (OBT), Spatial-Spraying (SPR), Personal Protection Measure (PPM)"
                )),
                fluidRow(
                  column(COVERAGE_LABELS_SIZE,h4("",align="center")),
                  column(COVERAGE_BAR_SIZE,h4("Coverage",align="center")),
                  column(COVERAGE_INIT_SIZE,h4("Time",align="center"))
                ),
                interventionRow("OVI", "OVIcov", "time_OVI_on"),
                # (FOG - fogging - is disabled in the GUI)
                interventionRow("LAR", "LARcov", "time_LAR_on"),
                interventionRow("BIO", "BIOcov", "time_BIO_on"),
                interventionRow("SRE", "SREcov", "time_SRE_on"),
                interventionRow("ITN", "ITNcov", "time_ITN_on", coverage=INITIAL_ITN_COVERAGE, time=INITIAL_ITN_TIME),
                interventionRow("IRS", "IRScov", "time_IRS_on"),
                interventionRow("ECS", "ECScov", "time_ECS_on"),
                interventionRow("ECT", "ECTcov", "time_ECT_on"),
                interventionRow("HOU", "HOUcov", "time_HOU_on"),
                interventionRow("OBT", "OBTcov", "time_OBT_on"),
                interventionRow("SPR", "SPRcov", "time_SPR_on"),
                interventionRow("PPM", "PPMcov", "time_PPM_on")
              ),
              # Plots (right-hand side): static ggplot figures, one under the other,
              # in the same order as the deployed app (https://skiware.shinyapps.io/simple/).
              # Each id is filled in server.r by renderModelPlots(); runSummary says what was run.
              mainPanel(
                tags$p(textOutput("runSummary"),style="font-size:90%;color:#555555;"),
                plotOutput("plotEIR"),
                plotOutput("plotTrajectory"),
                plotOutput("plotDemographics"),
                plotOutput("plotVC"),
                plotOutput("plotR0"),
                plotOutput("plotHuman")
              )
            ),
            tags$div(
              HTML('
              <br>
              <div>
                <center>
                  <img src="berkeley.jpeg" style=height:75px;">
                  <img src="ihi.jpeg" style=height:75px;">
                  <img src="ucsf.png" style=height:75px;">
                  <img src="itesm.png" style=height:75px;">
                </center>
              </div><br><br>
            '),
              helpText("Cite as: Cool paper here"),
              helpText("Contacts: <Model: Samson.Kiware@ucsf.edu, John.Marshall@berkeley.edu> <GUI: sanchez.hmsc@itesm.mx, tarrasch@berkeley.edu>"),
              helpText("CSS theme downloaded from: http://bootswatch.com (MIT licence)")
            )
  )
)
