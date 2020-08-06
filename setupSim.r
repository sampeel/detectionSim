# Test of gearGLM

setwd("/perm_storage/home/sampeel/chap2/sim2")

library(codetools)
library(raster)
library(spatstat)
library(multispeciesPP)                         # Github: sampeel version (not wfithian)
library(sp)
library(rgeos)
library(parallel)
library(vioplot)
library(RColorBrewer)
library(fields, quietly = TRUE)

source("bglm.r")
source("cellsObj.r")
source("dataObj.r")
source("domainObj.r")
source("estimateCoeffsNew.R")
source("plotLayers.r")
source("plottingObj.r")
source("resultsObject.r")
source("scenariosObject.r")
source("settingsObj.r")
source("simFuncs.r")
source("surveysObj.r")
source("utils.r")

#---------------------------------------------------------
# Settings (probably set once and leave type settings) ...
#---------------------------------------------------------

timeStartSetup <- Sys.time()

# Working directory (as tilde expansion may not work for me!)
homeDir <- "/perm_storage/home/sampeel/chap2/"

# Simulation space (for a laea projection)
simUnits <- "km"
simLonOrigin <- 180
simLatOrigin <- -72
simProj <- proj4str.laea(simUnits, simLonOrigin, simLatOrigin)
#simExt <- extent(-900, 620, -780, 620)                                     # on laea projection scale!

# True intensity functions stuff ...
numSpecies <- 21
namesSpecies <- NULL                                                        # set below.
alpha <- NULL                                                               # ln(lambda_ck) = alpha_k + beta_ik x_ci, estimated below.
beta <- NULL                                                                # estimated below.
gamma <- NULL                                                               # estimated below.
delta <- NULL                                                               # ln(bias_ck) = gamma_k + delta_j z_cj, estimated below.
lambda.formula <- "z ~ bath + I(log(chl)) + seaice"                         # log link assumed, area offset added if needed
bias.formula <- "z ~ bath"                                                  # log link assumed, area offset added if needed
isCentred <- TRUE                                                           # Centre the covariate data (e.g. xCentred = x_ci - mean(x_.i))

# Survey stuff ...
numSamples <- 2500                                                          # (scenario parameter!) Number of sample locations to simulate 
minSampleArea <- 0.001                                                      # in simUnits squared
maxSampleArea <- minSampleArea                                              # in simUnits squared
numSurveys <- 0                                                             # (scenario parameter!) Number of clusters in the sample locations, =0 gives random locations across whole domain
widthSurveys <- 5                                                           # Width of "square" each cluster potentially occupies, in number of cells (not used when numSurveys=0).

# Gear stuff ...
numGears <- 3
namesGears <- paste0("gear",1:numGears)
zeta <- NULL                                                                # ln(prDet) = zeta, estimated below.
gear.formula <- "~ factor(gear)"  
zetaMultiplier <- 1.0                                                       # (scenario parameter!) Alter effect of gears to be none (i.e. all gears have same prDet for a species) through to big effect (i.e. noticable difference between gears in the detection of a species) 
#   zeta[ ,k]' = zeta[ ,k]*zetaMult - mean(zeta[ ,k])(zetaMult-1)
detLogMean <- log(0.1)                                                      # logMean value for rlnorm function used to generate prDet values (i.e. d_kg ~ logNorm(mu,sd)) 
detLogSD <- 0.3                                                             # natural deviation from mean (this will be )
detSeed <- 1                                                                # RNG seed used to create probablity of detection starting values (NULL is no seed has been set).
skewedGearEffect <- FALSE
gearUseStrategyPA <- "rand"                                                 # How gear is assigned to sample sites ("rand" - uniform random, other options available)
gearUseStrategyPO <- NULL 

# Presence-only points ...
gammaMultiplier <- 1.0                                                      # (scenario parameter!) Alter "intensity" of the probability of observation per species to see more or less PO points.
deltaMultiplier <- 1.0                                                      # (scenario parameter!) Alter anthropogenic coefficient for probability of observation (between [0,1])

# Background points ...
numBGPoints <- 10000

# Minimum species location prevalence ...
minPrevalence <- 5                                                          # Minimum number of locations with evidence of species occupation (i.e. presences or non-zero abundances) before glms will be performed.

# Environment and mask RasterLayers (laea projection, assumes same as specified above).
# NB: assumes these cover the extent of the simulation (even if some of the values are NA).
# NB: x and y are assumed to be covariates, these are additional environment covariates.
envirDir <- paste0(homeDir, "Data/Environmental Rasters/")
envirFiles <- c("sim-rst_bathymetry.grd",
                "sim-rst_sst_summer_climatology.grd",
                "sim-rst_chl_summer_climatology.grd",
                "sim-rst_seaice_summer_variability.grd")
envirNames <- c("bath", "sst", "chl", "seaice")                             # use same names as in lambda.formula!
maskFile <- paste0(envirDir, "sim-Mask.grd")
maskValue <- NA

# Sample bias info files.
biasDir <-  envirDir
biasFiles <-  c("sim-rst_bathymetry.grd")
biasNames <- c("bath")                                                      # Use same names as in bias.formula

# Coastline data
# NB: data is SpatialPolygon in longlat projection (will be converted to simulation projection)
coastFile <-  paste0(homeDir, "Data/Coastlines/southOceanPoly.RData")

# Research base data, csv formatted table data with long, lat and species columns.
baseFile <-  paste0(homeDir,"Data/Research Stations/AntarcticBaseUniqueNames.csv")
baseXcolName <- "LONGITUDE"
baseYcolName <- "LATITUDE"
baseNamesColName <- "PLACE_NAME"
baseNamesShortColName <- "SHORT_NAME"

# Coefficient estimation data. csv formatted table data with long, lat and species columns.
coeffEstFilePO <-  paste0(homeDir, "Data/Ross Sea/RossSeaMaskPO.RData")
coeffEstFilePA <-  paste0(homeDir, "Data/Ross Sea/RossSeaMaskPA.RData")
coeffEstXcolName <- "x"
coeffEstYcolName <- "y"
coeffEstSpeciesColName <- "speciesID"
coeffEstAreaColName <- NULL                                                 # Give this or an area val to be used for all records!
coeffEstAreaVal <- 0.001                                                    # Give this or an area column name with data for all records!

# Input directory (where previous runs have saved processed data that can be reused).
inputDir <- paste0(getwd(), "/", "Input")
useSavedData <- TRUE                                                       # Use data from a previous run that is in the Input directory.

# Report progress and errors to this file, and all output (inc. plots) to this directory.
isOutputToFile <- TRUE                                                      # TRUE for file output, FALSE for console output
outDir <- paste0(getwd(), "/", "Output")                                    # Only active if isOutputToFile is TRUE
outFileName <- "setupMsg.txt"                                               # Name of file where messages go when not going to console.
outWarnOption <- 1                                                          # Option that sets the warning behaviour, see options("warn")
dataDumpDir <- "DataDumps"                                                  # Saves simulated data in the case of an error or warning from the SDM calls (assumes it is inside outDir)
resultsDir <- "Results"                                                     # Where the results are saved (assumes it is inside the outDir)

# Where to plot to (NB: isOutputToFile=TRUE should override if set to screen plotting.)
#plotToDevice <- "RStudioGD"                                                # "RStudioGD" for screen plotting in RStudio
plotToDevice <- "png"                                                       # "png" for plotting to a file of type png.
plotDir <- paste0(outDir, "/Plots")                                         # Only active if plotToDevice != "RStudio"
fileHeight <- 700
fileWidth <- 900

# Plotting - which plots are to be produced ...
doCheckPlots <- FALSE                                                       # Automatically sets all check plots to TRUE!  Helps with memory as plot layers aren't loaded until after sim is run!
checkPlotsAllRuns <- FALSE                                                  # FALSE gives plots on first run ONLY! TRUE gives plots for all runs.
envirData <- FALSE                                                          # Plot the environment covariate data (x in Fithian et al.)
maskData <- FALSE                                                           # plot the domain mask
biasData <- FALSE                                                           # plot the bias covariate data (z in Fithian et al.)
trueCoeffEstimates <- FALSE                                                 # plots for generation of 'true' coefficients
simDomain <- FALSE                                                          # called "domain" in plot object but can't be called that here!
intensity <- FALSE                                                          # lambda(x,y,alpha,beta,envirCovars)
probObs <- FALSE                                                            # bias(x,y,gamma,delta,biasCovars)
biasIntensity <- FALSE                                                      # lambda(x,y) * bias(x,y) is the intensity function for the PO simulation.
countPerCell <- FALSE                                                       # Number of simulated individuals per cell for each species
samplesPerCell <-  FALSE                                                    # Number of samples located in each cell (same for each species!)
countData <- FALSE                                                          # Number of species counted (i.e. abundunce) in each survey (well cell that the survey occurred in!)
PAData <- FALSE                                                             # Whether a species has been present (=0 in cell) or absent (=1 in cell) in a survey
POData <- FALSE                                                             # Number of observed individuals per cell for each species
BGPoints <- FALSE                                                           # Plot of the location of the generated background points (same for all runs and all species)
doResultPlots <- FALSE                                                      # Automatically do all result plots if this is TRUE.
estCoeffs <- FALSE
mseCompare <- FALSE
avgEstIntensity <- FALSE

#------------------------------------------------
# End of Settings
#------------------------------------------------

# Load last sim settings.
savedLastSimSettings <- makeFileName("LastSimSettings", inputDir, "RData")
if ( file.exists(savedLastSimSettings) ) {
  # Load data from previously saved R data file.
  load(savedLastSimSettings)
  
  # Probably an older version of simObj.
  if (  ! is.simObj(simObj.last) ) simObj.last <- initSimulation()
  
  
} else {
  # Simulation has not been previously run.
  simObj.last <- initSimulation()
}


#------------------------------------------------
# Initialise settings/simulation object for use in this simulation ...
#------------------------------------------------

simObj <- initSimulation(numRuns = 1, numCoresUse = 1)
simObj <- setSimOutput(simObj, inputDir, outFileName, outDir, isOutputToFile, dataDumpDir,
                       resultsDir, outWarnOption)


#------------------------------------------------
# Initialise plotting object for use in this simulation ...
#------------------------------------------------

# Initialise plotting ...
plotObj <- initPlotting(doCheckPlots, checkPlotsAllRuns)
plotObj <- setPlotsCheck(plotObj, envirData, maskData, biasData, simDomain, intensity,
                         biasIntensity, countPerCell, samplesPerCell, countData,
                         PAData, POData, BGPoints, probObs, trueCoeffEstimates)
plotObj <- setPlotsResults(plotObj, estCoeffs, mseCompare, avgEstIntensity)
plotObj <- setPlotsToFile(plotObj, device = plotToDevice, fileWidth = fileWidth,
                          fileHeight = fileHeight, plotDir = plotDir)

# Do we need to read in the plotting layers yet? (They use lots of memory!)
if ( any(c(doCheckPlots, envirData, maskData, biasData, simDomain, intensity, biasIntensity,
           countPerCell, samplesPerCell, countData, PAData, POData, BGPoints,
           probObs, trueCoeffEstimates)) ) {
  # Get the coastline SpatialPolygons.
  load(coastFile)
  coastLine.longlat <- cm1
  rm(cm1)
  
  # Get the research base names.
  SavedBaseNames <- makeFileName("ResearchBaseNames", inputDir, "RData")
  if ( file.exists(SavedBaseNames) && useSavedData ) {
    # Load data from previously saved R data file.
    load(SavedBaseNames)
    
    # Check the data is from the same raw data file.
    if ( basesObj$files != baseFile ) {
      stop("Research base names file has changed. Remove saved data from input directory.")
    }
    
    message("Base names data has been loaded.")
    
  } else {
    # Read in and process raw data.
    basesObj <- readDataTable(baseFile)
    basesObj$data <- postProcessBaseData(x=basesObj$data[ ,baseXcolName],
                                         y=basesObj$data[ ,baseYcolName],
                                         base=basesObj$data[ ,baseNamesColName],
                                         baseShort = basesObj$data[ ,baseNamesShortColName])
    basesObj <- as.DataPoints(basesObj, xColName = "x", yColName="y")
    
    # Save processed data for use in future.
    save(basesObj, file=SavedBaseNames)
  }
  
  # Set in the plotting object.
  plotObj <- setPlotsLayers(plotObj, simProj, coastLine.longlat, basesObj$data)
  
  # Remove objects that aren't needed anymore.
  rm(coastLine.longlat, basesObj)
}

#------------------------------------------------
# Get covariate data for use in this simulation ...
#------------------------------------------------

# Get the environmental covariates and mask.
SavedEnvirCovariates <- makeFileName("EnvirCovariates", inputDir, "RData")
if ( file.exists(SavedEnvirCovariates) && useSavedData && isCentred == simObj.last$isCentred) {
  # Load data from previously saved R data file.
  load(SavedEnvirCovariates)
  
  # Check the data is from the same raw data file.
  if ( ! all(envirObj$files == envirFiles) ) {
    stop("Environmental covariate files have changed. Remove saved data from input directory.")
  }
  
  # Check the data has the same centring setting.
  if ( envirObj$isCentred != isCentred ) {
    stop("Environmental covariates need to be changed. Remove saved data from input directory.")
  }
  
  message("Environmental data has been loaded.")
  
} else {
  # Read in the environmental covariates from the given files.
  envirObj <- readDataStack(envirFiles, envirDir, envirNames, TRUE, TRUE)
  
  # Check projection.
  if ( ! compareCRS(crs(envirObj$data), simProj) ) {
    stop("Simulation projection is not the same as environment data projection.")
  }
  
  # Plot environment data.
  if ( plotObj$envirData ) {
    fileName <- makeFileName("CheckEnvirData", plotObj$plotDir, plotObj$device)   #paste("CheckEnvirData.", plotObj$device, sep="")
    plot.plotLayers(envirObj$data, asp=1, device=plotObj$device, fileName=fileName,
                    fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                    main="Environmental covariate data")
  }
  
  # Read in the domain mask (assumes same projection and extent)
  rlMask <- raster(maskFile)
  envirObj <- setDataMask(envirObj, rlMask, maskValue)
  
  # Check projection, extent, etc.
  if ( ! compareRaster(envirObj$data[[1]], rlMask, extent=TRUE, rowcol=TRUE, crs=TRUE,
                       stopiffalse = FALSE) ) {
    stop("Environmental rasters and domain mask do not have the same raster structure.")
  }
  
  # Plot mask.
  if ( plotObj$maskData ) {
    fileName <- makeFileName("CheckMaskData", plotObj$plotDir, plotObj$device)   #paste("CheckMaskData.", plotObj$device, sep="")
    plot.plotLayers(envirObj$mask, asp=1, device=plotObj$device, fileName=fileName,
                    fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                    main="Mask for the simulation", sub="white cells are not included")
  }
  
  # Check for NA values in environment covariate data within domain.
  maskVals <- rlMask[1:ncell(rlMask)]
  validCellNums <- which(!is.na(maskVals))
  for ( layer in names(envirObj$data) ) {
    envirVals <- envirObj$data[[layer]][validCellNums]
    if ( any(is.na(envirVals)) ) {
      stop(paste0("NA values found within domain in environmental covariate ", layer))
    }
  }
  
  # Centre data if necessary and overwrite existing data.
  envirObj$isCentred <- isCentred
  if ( isCentred ) {
    # Save current data.
    envirObj$dataUncentred <- envirObj$data
    
    # Get data as it is used by the model.
    envirCovars <- getValuesMask(envirObj$data, rlMask, maskValue)
    sdm.formula <- delete.response.formula(lambda.formula)
    modelCovars <- model.frame(sdm.formula, as.data.frame(envirCovars))
    
    # Centre the data and save over current data (uncentred version).
    modelMeans <- apply(modelCovars,2,mean)
    centredModelCovars <- modelCovars
    numData <- length(modelMeans)
    for ( i in 1:numData ) {
      # Centre model data.
      centredModelCovars[ ,i] <- centredModelCovars[ ,i] - modelMeans[i]
      
      # Save the centred data.  
      rlSave <- setValuesMask(centredModelCovars[ ,i], rlMask, maskValue)
      if ( i == 1 ) {
        envirObj$data <- stack(rlSave)
      } else {
        envirObj$data <- addLayer(envirObj$data, rlSave)
      }
    }
    names(envirObj$data) <- paste0("data",1:numData)
    
    # Reset formula to reflect that data has already been processed.
    lambda.formula.old <- lambda.formula
    lambda.formula <- paste(names(envirObj$data), collapse = " + ")
    lambda.formula <- paste0("z ~ ", lambda.formula)
    
    # Save the centring values.
    envirObj$modelDataMeans <- modelMeans
  } else {
    envirObj$modelDataMeans <- NULL
  }
  
  # Save processed data for use in future.
  save(envirObj, lambda.formula, file=SavedEnvirCovariates)
  message("Environmental data has been made and saved.")
}

# Reset the simulation extent to equal that produced by the projection and cropping of data.
simExtReal <- extent(envirObj$data)

# Create sample bias covariates for this scenario (making sure extent covers that of simulation, and no NAs!).
SavedBiasCovariates <- makeFileName("BiasCovariates", inputDir, "RData")
if ( file.exists(SavedBiasCovariates) && useSavedData && isCentred == simObj.last$isCentred) {
  # Load data from previously saved R data file.
  load(SavedBiasCovariates)
  
  # Check the data is from the same raw data file.
  if ( ! all(biasObj$files == biasFiles) ) {
    stop("Bias covariate data has changed. Remove saved data from input directory.")
  }
  
  # Check the data has the same centring setting.
  if ( biasObj$isCentred != isCentred ) {
    stop("Bias covariates need to be changed. Remove saved data from input directory.")
  }
  
  message("Sample bias data has been loaded.")
  
} else {
  # Read in the environmental covariates from the given files.
  # NB: x and y are also part of the sample bias covariates but not stored in stack (so that
  #     all values of x and y are possible, not just the centroids of the raster cells).
  biasObj <- readDataStack(biasFiles, biasDir, biasNames, TRUE, TRUE)
  
  # Check raster structure same as environment data (and by default, mask)
  if ( ! compareRaster(envirObj$data[[1]], biasObj$data[[1]], extent=TRUE, rowcol=TRUE, crs=TRUE,
                       stopiffalse = FALSE) ) {
    stop("Environmental rasters and sample bias info do not have the same raster structure.")
  }
  
  # Plot sample bias info.
  if ( plotObj$biasData ) {
    fileName <- makeFileName("CheckBiasData", plotObj$plotDir, plotObj$device)
    plot.plotLayers(biasObj$data, asp=1, device=plotObj$device, fileName=fileName,
                    fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                    main="Sample bias covariate info")
  }
  
  # Read in the domain mask (assumes same projection and extent)
  rlMask <- raster(maskFile)
  biasObj <- setDataMask(biasObj, rlMask, maskValue)
  
  # Check projection, extent, etc.
  if ( ! compareRaster(biasObj$data[[1]], rlMask, extent=TRUE, rowcol=TRUE, crs=TRUE,
                       stopiffalse = FALSE) ) {
    stop("Sample bias info and domain mask do not have the same raster structure.")
  }
  
  # Plot mask.
  if ( plotObj$maskData ) {
    fileName <- makeFileName("CheckMaskData", plotObj$plotDir, plotObj$device)   #paste("CheckMaskData.", plotObj$device, sep="")
    plot.plotLayers(biasObj$mask, asp=1, device=plotObj$device, fileName=fileName,
                    fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                    main="Mask for the simulation", sub="white cells are not included")
  }
  
  # Check for NA values in bias covariate data within domain.
  maskVals <- rlMask[1:ncell(rlMask)]
  validCellNums <- which(!is.na(maskVals))
  for ( layer in names(biasObj$data) ) {
    biasVals <- biasObj$data[[layer]][validCellNums]
    if ( any(is.na(biasVals)) ) {
      stop(paste0("NA values found within domain in sample bias info ", layer))
    }
  }
  
  # Centre data if necessary and overwrite existing data.
  biasObj$isCentred <- isCentred
  if ( isCentred ) {
    # Save current data.
    biasObj$dataUncentred <- biasObj$data
    
    # Get data as it is used by the model.
    biasCovars <- getValuesMask(biasObj$data, rlMask, maskValue)
    sdm.formula <- delete.response.formula(bias.formula)
    modelCovars <- model.frame(sdm.formula, as.data.frame(biasCovars))
    
    # Centre the data and save over current data (uncentred version).
    modelMeans <- apply(modelCovars,2,mean)
    centredModelCovars <- modelCovars
    numData <- length(modelMeans)
    for ( i in 1:numData ) {
      # Centre
      centredModelCovars[ ,i] <- centredModelCovars[ ,i] - modelMeans[i]
      
      # Save the centred data.  
      rlSave <- setValuesMask(centredModelCovars[ ,i], rlMask, maskValue)
      if ( i == 1 ) {
        biasObj$data <- stack(rlSave)
      } else {
        biasObj$data <- addLayer(biasObj$data, rlSave)
      }
      
    }
    names(biasObj$data) <- paste0("bias",1:numData)
    
    # Reset formula to reflect that data has already been processed.
    bias.formula.old <- bias.formula
    bias.formula <- paste(names(biasObj$data), collapse = " + ")
    bias.formula <- paste0("z ~ ", bias.formula)
    
    # Save the centring values.
    biasObj$modelDataMeans <- modelMeans
  } else {
    biasObj$modelDataMeans <- NULL
  }
  
  # Save processed data for use in future.
  save(biasObj, bias.formula, file=SavedBiasCovariates)
  message("Sample bias data has been made and saved.")
  
}


#------------------------------------------------
# Make "true" coefficients for use in this simulation ...
#------------------------------------------------

# Can we load the data?
SavedTrueCoefficients <- makeFileName("TrueCoefficients", inputDir, "RData")
if ( file.exists(SavedTrueCoefficients) && useSavedData && isCentred == simObj.last$isCentred &&
     detLogMean == simObj.last$detLogMean && detLogSD == simObj.last$detLogSD && 
     detSeed == simObj.last$detSeed  ) {
  # Load data from previously saved R data file.
  load(SavedTrueCoefficients)
  
  # Check the data is from a simulation with the same settings. 
  if ( ( length(alpha) != numSpecies ) || ( length(gamma) != numSpecies) ||
       ( dim(beta)[1] != (numCoefficients(as.formula(lambda.formula)) - 1) ) ||
       ( length(delta) != (numCoefficients(as.formula(bias.formula)) - 1) )  ||
       ( dim(zeta)[1] != numGears ) ) {
    stop("True coefficients have changed. Remove saved data from input directory or set useSavedData = FALSE.")
  }
  
  message("True coefficients have been loaded.")
  
} else {
  # Make the domain, necessary for estimation of coefficients!
  domainObj <- makeDomain(envirObj$mask, envirObj$maskValue)
  
  # Get the PO and PA data (i.e. read actual data in from files)
  coeffEstData <- getPOPAData(domainObj, coeffEstFilePA, coeffEstFilePO, coeffEstXcolName,
                              coeffEstYcolName, coeffEstSpeciesColName, coeffEstAreaColName,
                              coeffEstAreaVal)
  
  # Check species number.
  if ( (dim(coeffEstData$PA)[2]-3 != numSpecies) ||
       (length(unique(coeffEstData$PO$species)) != numSpecies ) ) {
    stop("Specified number of species is not the same as in true PA and PO data.")
  }
  
  # Check the prevalence of each species.
  namesSpecies <- names(coeffEstData$PA)[c(-1,-2,-3)]
  prevLst <- checkPrevalence(namesSpecies, coeffEstData$PA[c(-1,-2,-3)], coeffEstData$PO, 
                             minPrevalence)
  if ( length(prevLst$namesSpeciesPA) != length(namesSpecies) )
    stop("All species in actual PA data must reach minimum prevalence value.")
  if ( dim(coeffEstData$PO)[1] > length(prevLst$whichPORows) )
    stop("All species in actual PO data must reach minimum prevalence value.")
  
  # Randomly generate probability of detection.
  set.seed(detSeed)
  prDet <- matrix(rlnorm(numGears*numSpecies, detLogMean, detLogSD), 
                  nrow=numGears, 
                  ncol=numSpecies,
                  dimnames = list(namesGears, namesSpecies))
  if ( skewedGearEffect ) {
    bigLogMean <- log(exp(detLogMean)*2)
    prDet[1, ] <- rlnorm(numSpecies, bigLogMean, detLogSD)
  }
  set.seed(NULL)
  #  if ( any(apply(prDet, 2, sum) > 1.0 ) ) 
  #    stop("One or more 'sum of gear probability of detections' are greater than one.")
  zeta <- log(prDet)                              # i.e. ln(d_kg) <- zeta_kg
  
  # Make the species intensity coefficients.
  coeffLst <- makeSpeciesIntensityCoeffs(coeffEstData$PA, lambda.formula,
                                         envirObj$data, domainObj$mask )
  
  # Re-jig the species intensity intercept to account for probility of detection < 100% and save coefficients.
  alpha.old <- coeffLst$alpha
  alpha <- alpha.old - log(apply(prDet,2,mean))     # mean(prDet[ ,k]) = mean(exp(zeta[ ,k])) 
  beta <- as.data.frame(coeffLst$beta)
  namesSpecies <- names(beta)
  rownames(beta) <- paste0("beta",1:dim(beta)[1])
  
  # Plot checks.
  coeffs <- as.data.frame(rbind(alpha, beta))
  if ( plotObj$trueCoeffEstimates ) {
    for ( species in namesSpecies) {
      # This species estimated true coefficients.
      thisSpeciesCoeffs <- coeffs[ ,species]
      
      # This species calculated true intensity
      rlIntensity <- rasterFromFunc(domainObj$ext, res(domainObj$mask), domainObj$mask,
                                    domainObj$maskValue, lambda.xy, as.formula(lambda.formula),
                                    thisSpeciesCoeffs, envirObj$data)
      
      
      # Make the plot layers (including the actual data points).
      indPresence <- which(coeffEstData$PA[ ,species] == 1)
      intensityLayers <- makeIntensityLayers(plotObj, rlIntensity,
                                             dataPoints = coeffEstData$PA[indPresence, c("x","y")],
                                             dataPointsPch = "+", dataPointsCol = "red")
      intensityLayers <- addPlotData(intensityLayers, coeffEstData$PA[-indPresence,c("x","y")],
                                     plotfunc = "points", pch=".", col="black")
      
      # Plot the plot layers
      fileName <- makeFileName(paste0("CheckTrueIntensity-",species),
                               plotObj$plotDir, plotObj$device)
      plot.plotLayers(intensityLayers, asp=1, device=plotObj$device, fileName=fileName,
                      fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                      main=paste0("Calculated true intensity for ", species))
    }
  }
  
  # Make the species probability of observation coefficients.
  coeffLst <- makeProbabilityObsCoeffs(coeffEstData$PO, bias.formula, biasObj$data,
                                       lambda.formula, envirObj$data, alpha.old, beta,
                                       domainObj, numBGPoints)
  gamma <- coeffLst$gamma
  delta <- coeffLst$delta

  # Save processed data for use in future.
  save(alpha, beta, gamma, delta, zeta, namesSpecies, namesGears, file=SavedTrueCoefficients)
  message("True coefficients have been estimated and saved.")
  
  # Plot check.
  if ( plotObj$trueCoeffEstimates ) {
    speciesNames <- unique(coeffEstData$PO$species)
    for ( species in speciesNames ) {
      # Calculate the probability of observation from this species coefficients.
      thisSpeciesCoeffs <- c(gamma[species], delta)
      rlProbObs <- rasterFromFunc(domainObj$ext, res(domainObj$mask), domainObj$mask,
                                  domainObj$maskValue, lambda.xy, as.formula(bias.formula),
                                  thisSpeciesCoeffs, biasObj$data)
      
      # Make the plot layers (including the actual data points).
      intensityLayers <- makeIntensityLayers(plotObj, rlProbObs,
                                             dataPoints = coeffEstData$PO[ , c("x","y")],
                                             dataPointsPch = ".", dataPointsCol = "black")
      indPresence <- which(coeffEstData$PO$species == species)
      intensityLayers <- addPlotData(intensityLayers, coeffEstData$PO[indPresence,c("x","y")],
                                     plotfunc = "points", pch="+", col="red")
      
      # Plot the plot layers
      fileName <- makeFileName(paste0("CheckTrueProbObs-",species),
                               plotObj$plotDir, plotObj$device)
      plot.plotLayers(intensityLayers, asp=1, device=plotObj$device, fileName=fileName,
                      fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                      main=paste0("Calculated probability of observation for ", species))
      
      # Recalculate species intensity.
      thisSpeciesCoeffs <- coeffs[ ,species]
      rlIntensity <- rasterFromFunc(domainObj$ext, res(domainObj$mask), domainObj$mask, domainObj$maskValue,
                                    lambda.xy, as.formula(lambda.formula), thisSpeciesCoeffs, envirObj$data)
      rlSampBias <- rlIntensity * rlProbObs
      
      # Make the plot layers (including the actual data points).
      intensityLayers <- makeIntensityLayers(plotObj, rlSampBias,
                                             dataPoints = coeffEstData$PO[indPresence, c("x","y")],
                                             dataPointsPch = "+", dataPointsCol = "red")
      
      # Plot the plot layers
      fileName <- makeFileName(paste0("CheckTrueBiasedIntensity-",species),
                               plotObj$plotDir, plotObj$device)
      plot.plotLayers(intensityLayers, asp=1, device=plotObj$device, fileName=fileName,
                      fileHeight = plotObj$fileHeight, fileWidth = plotObj$fileWidth,
                      main=paste0("Calculated sample biased intensity for ", species))
    }
  }
}


#------------------------------------------------
# Initialisation of simulation and results objects ...
#------------------------------------------------

# Further setting of simulation object items.
simObj <- setSimGeography(simObj, simExtReal, simProj)
simObj <- setSimSpecies(simObj, namesSpecies)
simObj <- initSimCoeffs(simObj, beta, alpha, gamma, delta, zeta)
simObj <- setSimModel(simObj, lambda.formula, bias.formula, gear.formula, isCentred)
simObj <- setSimSurveyInfo(simObj, numSamples, minSampleArea, maxSampleArea,
                           numSurveys, widthSurveys)
simObj <- setSimGearInfo(simObj, zetaMultiplier = zetaMultiplier, 
                         namesGears = namesGears, detLogMean = detLogMean, 
                         detLogSD = detLogSD, detSeed = detSeed, 
                         skewedGearEffect = skewedGearEffect, 
                         gearUseStrategyPA = gearUseStrategyPA,
                         gearUseStrategyPO = gearUseStrategyPO)
simObj <- setSimPOPoints(simObj, gammaMultiplier, deltaMultiplier)
simObj <- setSimBGPoints(simObj, numBGPoints)
simObj <- setSimUseSDMs(simObj, useSDM = c("PA","MsPP","Gear"), minPrevalence)

# Save the settings for comparison next time.
simObj.last <- simObj
save(simObj.last, file=makeFileName("LastSimSettings", inputDir, "RData"))
message("Simulation settings object created and saved.")
rm(simObj.last)

#----------------------------------
# Create scenario independent stuff
#----------------------------------


# Make the domain.
domainObj <- makeDomain(envirObj$mask, envirObj$maskValue)
if ( plotObj$domain ) plotDomain(plotObj, domainObj, withOwin=TRUE)
message("Domain object made.")

# Make the cells in the domain.
cellsObj <- makeCells(domainObj$mask, domainObj$maskValue)
cellsObj <- makeCovarData(cellsObj, envirObj$data, biasObj$data)
cellsObj <- makeIntensity(cellsObj, simObj$lambda.formula, simObj$initBeta, simObj$namesSpecies)
if ( plotObj$intensity )
  plotCellValues(plotObj, domainObj$mask, cellsObj$cells, cellsObj$trueLambda[ ,"sp351"], 
                 fileNameStart = "CheckCellIntensity",  titleTxts = "True intensity for ")
message("Cell object partially made.")

# Make the background data.
BG <- makeSimBG(simObj$numBGPoints, domainObj$owin, domainObj$mask)
BG[ ,names(cellsObj$covars)] <- getCellVals(cellsObj, BG$cell)
if ( BGPoints || doCheckPlots ) plotBGPoints(plotObj, BG, domainObj, cellsObj$cells)
message("Background points made.")

# Turn off the output to file for the setup phase.
timeEndSetup <- Sys.time()
doSimCleanUp(simObj)


# Garbage collection.
gc()
