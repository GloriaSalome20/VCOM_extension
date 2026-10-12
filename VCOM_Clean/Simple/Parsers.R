########################################################################
# Parsers.R
# Import of a parameter file (CSV / XLS) into a theta object.
# (Single copy of these functions; the older, out-of-date copy that lived in
#  ODEAuxiliaryFunctions.R has been removed.)
########################################################################

# Names of all the parameters a parameter file must define (same names and order as theta)
REQUIRED_PARAMETERS_LIST_GLOBAL = c(
    "beta","muEL","muLL","muPL","durEL","durLL","durPL","durEV",
    "gamma","tau1","tau2","ITNcov","IRScov","ECScov","ECTcov","HOUcov",
    "SPRcov","PPMcov","OBTcov","SREcov","BIOcov","LARcov","ATSBcov","SSPcov",
    "OVIcov","time_ITN_on","time_IRS_on","time_HOU_on","time_ECS_on","time_ECT_on","time_SPR_on","time_PPM_on",
    "time_OBT_on","time_BIO_on","time_SRE_on","time_LAR_on","time_OVI_on","time_ATSB_on","time_SSP_on","muV",
    "Q0","phiB","phiI","rITN","sITN","rIRS","sIRS","rHOU",
    "sHOU","rSPR","sSPR","rPPM","sPPM","rECS","sECS","rECT",
    "sECT","dIRS","dHOU","aOBT","eSRE","fBIO","fLAR","fATSB",
    "fSSP","fOVI","eOVI","f0","epsilon0","recRate","NH_eq","bV","bh"
  )#,"lambdaV","NV_eq")

parseImportedCSVParameters = function(inputDataFrame){
  #. parseImportedCSVParameters: Once a csv file has been imported this function converts the data into a theta object
  values <- inputDataFrame[REQUIRED_PARAMETERS_LIST_GLOBAL, 1]
  names(values) <- REQUIRED_PARAMETERS_LIST_GLOBAL
  values
}
importCSVParametersFromDirectory = function(fileName){
  #. importCSVParametersFromDirectory: Loads the verified parameters of the simulation from a CSV file into the workflow
  inputDataFrame=read.csv(file=fileName,head=FALSE,sep=",",row.names=1,col.names=c("ID","Value","Description"),fill=FALSE,stringsAsFactors=FALSE)
  inputDataFrame[,1]=as.numeric(inputDataFrame[,1])
  return(validateCSVParameters(inputDataFrame))
}
validateCSVParameters = function(inputDataFrame){
  #. validateCSVParameters: Verifies that the user has imported a valid CSV file with all the parameters defined
  containsNA=anyNA(inputDataFrame[REQUIRED_PARAMETERS_LIST_GLOBAL,1])
  if(containsNA==FALSE){
    print("FILE IS VALID")
    return(inputDataFrame)
  }else{
    print("ERROR IN FILE")
    sub=subset(inputDataFrame,is.na(inputDataFrame))
    return(sub)
  }
}
importXLSParametersFromDirectory = function(fileName){
  #. importXLSParametersFromDirectory: Loads the verified parameters of the simulation from an XLS file into the workflow
  #  (needs the XLConnect package)
  xlsRaw = loadWorkbook(fileName)
  inputDataFrame = readWorksheet(xlsRaw,sheet=1,head=FALSE,startCol=0,rownames=1)
  inputDataFrame[,1]=as.numeric(inputDataFrame[,1])
  return(validateCSVParameters(inputDataFrame))
}
importCSVXLSParametersFromDirectoryShiny = function(fileName,typeString){
  #. importCSVXLSParametersFromDirectoryShiny: Wrapper function that works only within shiny to select the proper import function to use on a given file
  print(typeString)
  if(typeString == "text/csv"){return(importCSVParametersFromDirectory(fileName))}
  if(typeString == "application/vnd.ms-excel"){return(importXLSParametersFromDirectory(fileName))}
}
