#------------------------------------------------------------------------------##
#################################################################################
## Hector M. Sanchez C. (sanchez.hmsc@itesm.mx)                               ###
## Malaria vector ODE model GUI                                               ###
## 02/May/2016                                                                ###
#################################################################################
#------------------------------------------------------------------------------##

#################################################################################
#---# MAIN SHINY SERVER APPLICATION #############################################
shinyServer(
  function(input,output,session){

    #############################################################################
    # PLOTS #####################################################################
    # The six plots on the right of the page (they read the global IVM_traj).
    # Static ggplot figures, as in the deployed app (https://skiware.shinyapps.io/simple/);
    # the matching slots (plotOutput) are in ui.r.
    renderModelPlots <- function(){
      output$plotEIR = renderPlot({plotEIR(IVM_traj)})
      output$plotTrajectory = renderPlot({plotTrajectory(IVM_traj)})
      output$plotDemographics = renderPlot({barChartMosquitoDemographics_slwu(IVM_traj)})
      output$plotVC = renderPlot({plotVC(IVM_traj)})
      output$plotR0 = renderPlot({plotR0(IVM_traj)})
      output$plotHuman = renderPlot({plotTrajectoryHumans(IVM_traj)})
      # The interactive plotly versions (plotTrajectoryPlotLy, plotVCPlotLy, plotR0PlotLy, ...) are
      # still in Plots.R; to use one, swap plotOutput -> plotlyOutput in ui.r and renderPlot -> renderPlotly here.
    }

    #############################################################################
    # PRIMING GUI ###############################################################
    output$modelStatus <- renderText({"Waiting"})
    output$runSummary <- renderText({RUN_SUMMARY})      # what the first run (ui.r) used
    # (the download buttons can be disabled at start with shinyjs::disable("downloadCSVTrace"), ...)
    renderModelPlots()

    #############################################################################
    # CLICK EVENTS (console logging) ############################################
    observeEvent(input$buttonTest,{cat("Button event!\n")})
    observeEvent(input$radioSpecies,{cat("Radio event!\n")})
    observeEvent(input$sliderTime,{cat("Slider event!\n")})

    #############################################################################
    # RUN MODEL #################################################################
    # The value of a numeric box, 0 if it was left empty
    valueOrZero <- function(x){ if (is.null(x) || is.na(x)) 0 else x }

    observeEvent(input$buttonRun,{
      output$modelStatus <- renderText({"Calculating. Please wait."})
      cat("Run button event!\n")
      # COVERAGE AND START DAY OF EACH INTERVENTION ------------------------------
      coverage <- setNames(vapply(APP_INTERVENTIONS, function(code) valueOrZero(input[[paste0(code,"cov")]]), numeric(1)), APP_INTERVENTIONS)
      start    <- setNames(vapply(APP_INTERVENTIONS, function(code) valueOrZero(input[[paste0("time_",code,"_on")]]), numeric(1)), APP_INTERVENTIONS)
      # RUN: theta, initial state and the ODE, as in Driver.R (SeasonalSetup.R) -------
      #   climate model, seasonal K(t), temperature-dependent parameters and burn-in come from the sidebar
      result <- tryCatch({
        runVCOM(
          species=input$radioSpecies, EIR=input$radioEIR, days=input$sliderTime,
          coverage=coverage, start=start,
          climateModel=input$climateModel, seasonalK=isTRUE(input$seasonalK),
          temperatureGroups=if (is.null(input$temperatureGroups)) character(0) else input$temperatureGroups,
          burnInYears=input$burnInYears
        )
        "ok"
      }, error=function(e){ cat("Run failed:", conditionMessage(e), "\n"); conditionMessage(e) })
      if (!identical(result,"ok")) {
        output$modelStatus <- renderText({paste("The model could not be run:", result)})
        failure <- paste("The model could not be run:", result, "(the plots show the previous run)")
        output$runSummary <- renderText({failure})
        return()
      }
      print(TRANSMISSION_PARAMETERS)
      renderModelPlots()
      print(INTERVENTION_PARAMETERS)
      runSummary <- RUN_SUMMARY
      output$runSummary <- renderText({runSummary})
      #--------------------------------------------------------------------------
      shinyjs::enable("downloadCSVTrace"); shinyjs::enable("downloadPlotTrace")
      shinyjs::enable("downloadCSVEIR"); shinyjs::enable("downloadPlotEIR")
      shinyjs::enable("downloadPlotDemographics"); shinyjs::enable("downloadPlotVC")
      shinyjs::enable("downloadPlotHuman"); shinyjs::enable("downloadPlotR0")
      output$modelStatus <- renderText({"Done!"})
    })

    #############################################################################
    # DOWNLOADS HANDLERS ########################################################
    output$downloadParameters <- downloadHandler(
      filename = function(){paste("VCOM_Parameters","csv",sep=".")},
      content = function(file){
        df=data.frame(theta)
        df$Description=0
        names(df)=c("Value","Description")
        write.csv(df,file)
      }
    )
    output$downloadCSVTrace <- downloadHandler(
      filename = function(){paste("VCOM_Trace","csv",sep=".")},
      content = function(file){
        write.csv(IVM_traj,file)
      }
    )

    # PNG download of a ggplot built from the current trajectory (twice as wide as high)
    pngDownloadHandler <- function(fileSuffix, makePlot){
      downloadHandler(
        filename = function(){paste(input$dataset,fileSuffix,sep='')},
        content = function(file){
          device <- function(...,width,height){grDevices::png(...,width=2*width,height=height,res=300,units="in")}
          ggsave(file,plot=makePlot(IVM_traj),device=device)
        }
      )
    }
    output$downloadPlotTrace <- pngDownloadHandler('TrajectoryPlot.png', plotTrajectoryShiny)
    output$downloadPlotEIR <- pngDownloadHandler('EIRPlot.png', plotEIRShiny)
    output$downloadPlotVC <- pngDownloadHandler('VC.png', plotVC)
    output$downloadPlotHuman <- pngDownloadHandler('Human.png', plotTrajectoryHumans)
    output$downloadPlotR0 <- pngDownloadHandler('R0.png', plotR0)
    output$downloadPlotDemographics <- pngDownloadHandler('Demographics.png', barChartMosquitoDemographicsShiny)

    # PDF report with the main plots
    output$downloadPlots <- downloadHandler(
      filename = function(){"VCOMReport.pdf"},
      content = function(file) {
        pdf(file, width=12, height=7)
        print(plotR0(IVM_traj))
        print(plotVC(IVM_traj))
        print(plotTrajectoryShiny(IVM_traj))
        # Other plots that can be added to the report:
        # print(barChartMosquitoDemographicsShiny(IVM_traj))
        # print(plotTrajectoryHumans(IVM_traj))
        # print(plotEIRShiny(IVM_traj))
        dev.off()
      }
    )
    #############################################################################
  })
#---#############################################################################
