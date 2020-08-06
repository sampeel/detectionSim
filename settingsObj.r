initSimulation <- function(numRuns=0, numCoresUse=0) {
  
  # Initialise or reset the simulation object.
  # Contains the bits and pieces that specify what sort of simulation is to be performed.
  
  # Create the simulation list object.
  simObj <- list(ext=NULL,             # Extent, domain or observation window of the simulation.
                 proj=NULL,            # Projection that the simulation is to be performed in.
                 numRuns=0,            # Number of simulation runs/repeats to perform.
                 numSpecies=0,         # Number of species to be simulated.
                 namesSpecies=NULL,    # Vector of character strings (numSpecies x 1)
                 numGears = 0,         # Number of gear types to be simulated.
                 namesGears = 0,       # Vector of character strings (numGears x 1)
                 lambda.formula=NULL,  # String that is species intensity formula (without covariates and offset) see "formula" help in R
                                       #     e.g. lambda ~ x1 + x2 + x3  
                 bias.formula=NULL,    # String that is sample bias formula (without covariates and offset) see "formula" help in R
                                       #     e.g. bias ~ z1 + z2 
                 gear.formula=NULL,    # String that is probability of detection formula
                                       #     e.g. ~ factor(gear) such that PA ~ x1 + x2 + x3 + factor(gear)
                 isCentred=FALSE,      # Whether or not to centre the model data (i.e. how the covariates are used in the formulae).
                 #coeffs=NULL,          # Log intensity function coefficients (numCoeffs x numSpecies)
                 initBeta=NULL,        # Log intensity function coefficients (numCoeffs x numSpecies)
                 initGamma=NULL,       # Log sample bias function coefficients (numSpecies x 1)
                 initDelta=NULL,       # Log sample bias function coefficients (numBiases x 1)
                 initZeta=NULL,        # Log probability detection function coefficients (numGear x numSpecies)
                 deltaMultiplier=1.0,  # Multiplier for delta (alters only the amount of anthropogenic sample bias)
                 gammaMultiplier=1.0,  # Multiplier for gamma (alters only the amount of species specific sample bias)
                 zetaMultiplier=1.0,   # Multiplier for zeta (new.zeta_kg = zeta_kg*zetaMult - mu_k(zetaMult - 1) where mu_k = mean(zeta_kg)).
                 beta=NULL,            # Log intensity function coefficients (numCoeffs x numSpecies)
                 gamma=NULL,           # Log sample bias function coefficients altered by multiplier (numSpecies x 1)
                 delta=NULL,           # Scenario log sample bias function coefficients altered by multiplier (numBiases x 1)
                 zeta=NULL,            # Log probability detection function coefficients altered by multiplier (numGear x numSpecies)
                 detLogMean = 0.0,     # Probability of detection created by log normal random number generator.
                 detLogSD = 0.0,       # as above, these are arguments for the rlnorm function.
                 detSeed = NULL,       # seed to set RNG for creation of probability of detection "true" values.
                 skewedGearEffect=FALSE,
                 gearUseStrategyPA="rand",
                 gearUseStrategyPO=NULL,
                 numSamples=NULL,      # Number of presence-absence surveys to simulate
                 minSampleArea=NULL,   # Minimum area for a survey (in square sim units, e.g. km^2)
                 maxSampleArea=NULL,   # Maximum area for a survey (set same as minimum for equal areas)
                 numSurveys=0,         # Number of clusters to group the surveys into, numSurveys=0 gives random sampling.
                 widthSurveys=NULL,   # Width of area that each cluster potentially occupies (in number of cells)
                 numBGPoints=NULL,     # Number of background points to generate for multispeciesPP
                 minPrevalence=0,      # Minimum number of locations with evidence of species occupation (i.e. presences or non-zero abundances)
                 useSDM.Count=FALSE,   # Use a glm with count data as the SDM
                 useSDM.Cloglog=FALSE, # Use a glm with presence-absence data and a cloglog as the SDM
                 useSDM.PPM=FALSE,     # Use the PPM function with presence-only data as the SDM
                 useSDM.MsPP=FALSE,    # Use the MultispeciesPP function with presence-absence and presence-only data as the SDM
                 useSDM.Gear=FALSE,    # Use the GearGLM function with presence-absence and presence-only data as the SDM
                 isOutputToFile=FALSE, # Whether or not to send messages to stdout and stderr, or, a file.
                 inputDir="",          # Input can come from this directory (for processed data that is used every time a scenario is tested)
                 outDir="",            # Output goes to this directory.
                 outFile="",           # Text output (eg. from functions message, error, warning, cat, print) goes to this file.
                 outWarnOption=0,      # The way warnings will be presented, 0 = in a group at end, 1 = immediately.
                 dataDumpDir="",       # The directory that will contain the saved data files if a warning or error occurs during an SDM run.
                 resultsDir="",        # The directory that will contain the saved results (in "RData" format)
                 numCoresUse=0,        # Number of cores to use in parallel processing (i.e. mclapply function).
                 #
                 # The following values are worked out by the code, no need to specifically set.
                 #
                 #domain=NULL,         # Replaced by domainObj!! valid points within the extent (i.e. owin object with boundary and mask)
                 numCoeffs=0,          # Saves number of coefficients (= number of rows in coeff)
                 numBiases=0,          # Saves number of sample bias covariates (from delta)
                 origWarnOption=NULL,  # Value of options("warn") when simulation starts, reset to this at end.
                 isError=FALSE         # Whether or not there has been an error within a function
  )
  
  # Set the number of simulations and number of cores to use.
  simObj$numRuns <- numRuns
  simObj$numCoresUse <- numCoresUse
  
  # Return value.
  return(simObj)
  
}  

#-----------------------------------------------------------------------------------------

is.simObj <- function(obj) {
  
  # Test whether the argument is a valid simulation object.  Returns true if it is, false otherwise.
  # NB: only tests names of items at this stage, not classes of items!
  
  # Check argument is the right class (as far as we can!)
  if ( !is(obj, "list") ) {
    # The argument is not even a list.  It is not a simulation object.
    return(FALSE)
  }
  
  # Get the expected names of the items for a simulation object.
  objectItemNames <- names(initSimulation())
  
  # Check the object has the same items.
  if ( all(names(obj) %in% objectItemNames) ) {
    # The same item names, hence, a valid simulation object.
    return(TRUE)
    
  } else {
    # Not the same item names, hence, an invalid simulation object.
    return(FALSE)
  }
  
}

#-----------------------------------------------------------------------------------------

setSimGeography <- function(simObj, extSim, projSim) {

  # Set the extent and projection of the simulation.
  simObj$ext <- extSim
  simObj$proj <- projSim
  
  # Return value.
  return(simObj)

}

#-----------------------------------------------------------------------------------------

setSimSpecies <- function(simObj, namesSpecies) {
  
  # Set the species info.
  
  numSpecies <- length(namesSpecies)
  if ( numSpecies > 0 ) {
    simObj$namesSpecies <- namesSpecies
    simObj$numSpecies <- numSpecies
  } else {
    simObj$isError <- TRUE
    stop("Names of species not set as none provided (or not provided in a vector).")
  }
  

  # Return value.
  return(simObj)
  
}
#-----------------------------------------------------------------------------------------

initSimCoeffs <- function(simObj, beta, alpha = rep(0,dim(beta)[[2]]), 
                          gamma = NULL, delta=NULL, zeta=NULL) {
  
  # Initialise the coefficients and the numbers of things related to each.  Models are:
  #     log(lambda[c,k]) = alpha[k] + sum(beta[ ,k] * X[c, ])
  #       log(bias[c,k]) = gamma[k] + sum(delta[ ] * Z[c, ])
  #        log(det[g,k]) = zeta[g,k]
  # where lambda is the species intensity, bias is the thinning of the species intensity
  # due to biased sampling, and det is the thinning of the species intensity due to the 
  # inefficiencies of the gear used to collect a sample and/or observation 
  # (c = cell, g = gear, k = species, X and Z are covariate matrices).  
  
  # Check species names and number have been set.
  numSpecies <- simObj$numSpecies
  namesSpecies <- simObj$namesSpecies
  if ( numSpecies == 0 || is.null(namesSpecies) ) {
    simObj$isError <- TRUE
    stop("Please set species names before continuing.")
  }
  
  # Check dimensions relating to species ...
  if ( ! match.names(names(alpha), namesSpecies) ) 
    stop("Names of species in 'alpha' does not match given species names.")
  if ( ! match.names(dimnames(beta)[[2]], namesSpecies) )
    stop("Names of species in 'beta' does not match given species names.")
  if ( ! match.names(names(gamma), namesSpecies) )
    stop("Names of species in 'gamma' does not match given species names.")
  if ( ! match.names(dimnames(zeta)[[2]], namesSpecies) )
    stop("Names of species in 'zeta' does not match given species names.")
  
  # Set other dimensions ...
  numBetas <- dim(beta)[1]
  simObj$numCoeffs <- numBetas + 1                   # Includes alpha!
  simObj$numGears <- dim(zeta)[1] 
  simObj$numBiases <- length(delta)
  
  # Initialise coefficients values ...
  simObj$initBeta <- rbind(alpha[namesSpecies], beta[ ,namesSpecies])
  simObj$initZeta <- zeta[ ,namesSpecies]
  simObj$initGamma <- gamma[namesSpecies]
  simObj$initDelta <- delta
    
  # Set coefficient names ...
  rownames(simObj$initBeta) <- c("alpha", paste0("beta", 1:numBetas))
  rownames(simObj$initZeta) <- paste0("zeta", 1:simObj$numGears)
  names(simObj$initDelta) <- paste0("delta", 1:simObj$numBiases)  
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimCoeffs <- function(simObj, beta, gamma, delta, zeta) {

  # Set the coefficient value for the simulation (NB: beta includes alpha in the first row.)
  # Assumes coefficients are in the same order as simObj$namesSpecies (for species dims).
  
  # Check species names and number have been set.
  namesSpecies <- simObj$namesSpecies
  if ( simObj$numSpecies == 0 || is.null(namesSpecies) ) {
    simObj$isError <- TRUE
    stop("Please set species names before continuing.")
  }
   
  # Check dimensions relating to species ...
  if ( ! match.names(dimnames(beta)[[2]], namesSpecies) )
    stop("Names of species in 'beta' does not match given species names.")
  if ( ! match.names(names(gamma), namesSpecies) )
    stop("Names of species in 'gamma' does not match given species names.")
  if ( ! match.names(dimnames(zeta)[[2]], namesSpecies) )
    stop("Names of species in 'zeta' does not match given species names.")
  
  # Check other dimensions ...
  if ( simObj$numCoeffs != dim(beta)[1] ) 
    stop("Number of 'beta' rows does not match given number of coefficients.")
  if ( simObj$numBiases != length(delta) ) 
    stop("Number of 'delta' elements does not match given number of biases.")
  if ( simObj$numGears != dim(zeta)[1] )
    stop("Number of 'zeta' elements does not match given number of gears.")

  # Set coefficients values ...
  simObj$beta <- beta[ ,namesSpecies]
  simObj$zeta <- zeta[ ,namesSpecies]
  simObj$gamma <- gamma[namesSpecies]
  simObj$delta <- delta
  
  # Set coefficient names ...
  numBetas <- simObj$numCoeffs - 1
  rownames(simObj$beta) <- c("alpha", paste0("beta", 1:numBetas))
  rownames(simObj$zeta) <- paste0("zeta", 1:simObj$numGears)
  names(simObj$delta) <- paste0("delta", 1:simObj$numBiases)  

  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------
  
setSimModel <- function(simObj, lambda.formula, bias.formula, gear.formula, isCentred=TRUE) {
  
  # Sets these values in the settings object.  

  # Check coefficients have already been set.
  if ( simObj$numCoeffs == 0 )
    stop("Please initialise coefficients before continuing.")

  # Set formulae.
  simObj$lambda.formula <- as.formula(lambda.formula)
  simObj$bias.formula <- as.formula(bias.formula)
  simObj$gear.formula <- as.formula(gear.formula)

  # Check consistency of intensity formula with previously set coefficients.
  terms.attributes <- attributes(terms(simObj$lambda.formula))
  if ( terms.attributes$intercept == 1) {
    numBetas <- simObj$numCoeffs - 1
  } else {
    simObj$isEffor <- TRUE
    stop("An intercept is currently assumed to be present in this code.")
  }
  if ( numBetas != length(terms.attributes$term.labels) ) {
    simObj$isError <- TRUE
    msg <- paste0("Argument 'lambda.formula' has the wrong number of covariate terms, should be = ",
                  numBetas, ".")
    stop(msg)
  }

  # Check consistency of sample bias formula with previously set coefficients.
  terms.attributes <- attributes(terms(simObj$bias.formula))
  if ( simObj$numBiases != length(terms.attributes$term.labels) ) {
    simObj$isError <- TRUE
    msg <- paste0("Argument 'bias.formula' has the wrong number of covariate terms, should be = ",
                  simObj$numBiases, ".", sep="")
    stop(msg)
  }
  isIntercept <- terms.attributes$intercept
  if ( !is.null(simObj$initGamma) == isIntercept ) {
    # All good.
  } else {
    simObj$isError <- TRUE
    stop("The presence or absence of a bias intercept is not consistent with previously set gamma.")
  }

  # Check consistency of gear formula with previously set coefficients.
  terms.attributes <- attributes(terms(simObj$bias.formula))
  if ( terms.attributes$intercept != 1 || length(terms.attributes$term.labels) != 1 )
    stop("Form of argument 'gear.formula' is not accepted in code at this time.")
  
  # Set whether or not to centre the model data.
  simObj$isCentred = isCentred
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimGearInfo <- function(simObj, detLogMean = 0.1, detLogSD = 0.3, zetaMultiplier = 1.0, 
                           namesGears=paste0("gear", 1:simObj$numGears), detSeed = NULL,
                           skewedGearEffect = FALSE, gearUseStrategyPA = "rand", 
                           gearUseStrategyPO = NULL) {
  
  # Sets these values in the settings object.  
  
  # Check and set the number of gears.
  if ( length(namesGears) != simObj$numGears ) {
    stop("Names of gears is not consistent with number of gears set previously.")
  }
  simObj$namesGears <- namesGears

  # Set formula and other bits and pieces.
  simObj$zetaMultiplier <- zetaMultiplier
  simObj$detLogMean <- detLogMean
  simObj$detLogSD <- detLogSD
  simObj$detSeed <- detSeed
  simObj$skewedGearEffect <- skewedGearEffect
  simObj$gearUseStrategyPA <- gearUseStrategyPA
  simObj$gearUseStrategyPA <- gearUseStrategyPO

  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimSurveyInfo <- function(simObj, numSamples, minSampleArea=1, maxSampleArea=minSampleArea, 
                             numSurveys=0, widthSurveys=Inf){
  
  # Sets the simulation list with these values.
  
  # Set survey information.
  simObj$numSamples <- numSamples
  simObj$minSampleArea <- minSampleArea
  simObj$maxSampleArea <- maxSampleArea
  simObj$numSurveys <- numSurveys
  simObj$widthSurveys <- widthSurveys
  
  # Check min and max sample area are in right order.
  if ( minSampleArea > maxSampleArea ) {
    simObj$isError <- TRUE
    stop("Minimum sample area is greater than maximum sample area, check settings.")
  }
  
  # Check that the number of surveys is positive.
  if ( numSurveys < 0 ) {
    simObj$isError <- TRUE
    stop("Unable to make a negative number of surveys for sampling locations, check settings.")
  }
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimBGPoints <- function(simObj, numBGPoints=10000) {
  
  # Sets the simulation list with these values.
  
  # Set background points information.
  simObj$numBGPoints <- numBGPoints

  # Check number of points is positive.
  if ( numBGPoints < 1 ) {
    simObj$isError <- TRUE
    stop("The number of background points must be a positive integer.")
  }
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimPOPoints <- function(simObj, gammaMultiplier=1.0, deltaMultiplier=1.0){
  
  # Sets the simulation list with these values.
  
  simObj$gammaMultiplier <- gammaMultiplier
  simObj$deltaMultiplier <- deltaMultiplier
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimUseSDMs <- function(simObj, useSDM = c("PA", "MsPP", "Gear"), minPrevalence=10) {
  
  # Sets which SDMs to use in the simulation list object.
  
  # Check there are only valid SDMs specified.
  useValidSDM <- match.arg(useSDM, several.ok = TRUE)   # Produces error for no matches!
  
  # Set use SDM items in list.
  simObj$useSDM.Count <- "AB" %in% useValidSDM
  simObj$useSDM.Cloglog <- "PA" %in% useValidSDM
  simObj$useSDM.PPM <- "PO" %in% useValidSDM
  simObj$useSDM.MsPP <- "MsPP" %in% useValidSDM
  simObj$useSDM.Gear <- "Gear" %in% useValidSDM
  
  # Check for any invalid SDM.
  invalidSDMs <- setdiff(useSDM, useValidSDM)
  if ( length(invalidSDMs) > 0 ) {
    warning("Sim was asked to use one or more invalid SDMs.  These were discarded.")
  }

  # Set minimum prevalence.
  simObj$minPrevalence <- minPrevalence
  
  # Return value.
  return(simObj)
  
}

#-----------------------------------------------------------------------------------------

setSimOutput <- function(simObj, inputDir=getwd(), outFileName="simOut.txt", outDir=getwd(), 
                         isOutputToFile=FALSE, dataDumpDir="DataDumps", resultsDir="Results",
                         outWarnOption=getOption("warn"), dirSep="/") {
  
  # Sets the location of output for the simulation runs.  Also, options for simulations runs.
  # Creates directory, if necessary.  Opens connection to file, if necessary.
  
  # Set the values that control the output location.
  simObj$inputDir <- inputDir
  simObj$outDir <- outDir
  simObj$outFile <- paste0(outDir, dirSep, outFileName)
  simObj$isOutputToFile <- isOutputToFile
  simObj$dataDumpDir <- paste0(outDir, dirSep, dataDumpDir)
  simObj$resultsDir <- paste0(outDir, dirSep, resultsDir)
  simObj$outWarnOption <- outWarnOption
  #simObj$hugeRndObs <- hugeRndObs

  # Get the current warn option and then reset to the one specified for the simulation.
  simObj$origWarnOption <- getOption("warn")    
  options(warn=simObj$outWarnOption)
  
  # Get the current warning level for huge random observations and then reset.
  #simObj$origHugeRndObs <- spatstat.options("huge.npoints")
  #spatstat.options(huge.npoints=simObj$hugeRndObs)
  
  # If necessary, open file ready for diversion of output from stdout and stderr.
  if ( isOutputToFile ) {
    # Create directory, if necessary.
    if ( ! dir.exists(outDir) ) {
      dir.create(outDir, recursive = TRUE)
    }
    
    # Open file for output.
    # NB: using file automatically creates new file or overwrites/clears existing file.
    outConn <- file(simObj$outFile, open="wt")
    
    # Start diversion of output from stdout and stderror to a file.
    sink(file = outConn, type="output")
    sink(file = outConn, type="message")
    
    # File opened message.
    message(paste("Output file opened at", Sys.time()))
  } else {
    message(paste("Output started at", Sys.time()))
  }

  # Create the input directory for processed data, if necessary (NB: same across scenarios).
  if ( ! dir.exists(inputDir) ) {
    dir.create(inputDir, recursive = TRUE)
  }
  
  # Create data dump directory, if necessary.
  if ( ! dir.exists(simObj$dataDumpDir) ) {
    dir.create(simObj$dataDumpDir, recursive = TRUE)
  }
  
  # Create results directory, if necessary.
  if ( ! dir.exists(simObj$resultsDir) ) {
    dir.create(simObj$resultsDir, recursive = TRUE)
  }
  
  # Return value.
  return(simObj)

}

#-----------------------------------------------------------------------------------------

doSimCleanUp <- function(simObj) {
  
  # Turn off diversion of output to file, if necessary.
  if ( simObj$isOutputToFile ) {
    # File closed message.
    message(paste("File closed at", Sys.time()))
    
    # Turn off diversion.
    sink(type="message")
    sink(type="output")
  } else {
    message(paste("Finished at", Sys.time()))
  }
  
  # Change the warning option back to whatever it was before the simulation started.
  options(warn=simObj$origWarnOption)
  
  # Reset the spatstat options back to whatever it was before the simulation started.
  #spatstat.options(huge.npoints=simObj$origHugeRndObs)                                           

  # Garbage collection.
  gc()
  
  # Return value.
  invisible(NULL)
  
}

#-----------------------------------------------------------------------------------------

setSimScenarioInfo <- function(simObj, numRuns, numCoresUse, numSamples, numSurveys, 
                               gammaMultiplier, deltaMultiplier, zetaMultiplier,
                               gearUseStrategyPA, gearUseStrategyPO, useSDM) {
  
  # Set the simulation parameters that are also scenario parameters.  
  
  # Execution settings ...
  simObj$numRuns <- numRuns
  simObj$numCoresUse <- numCoresUse
  
  # Scenario parameters ...
  simObj$numSamples <- numSamples
  simObj$numSurveys <- numSurveys
  simObj$gammaMultiplier <- gammaMultiplier
  simObj$deltaMultiplier <- deltaMultiplier
  simObj$zetaMultiplier <- zetaMultiplier
  
  # Create this sceanrio's coefficients from these parameters.
  simObj <- setSimCoeffs(simObj, 
                         simObj$initBeta, 
                         simObj$initGamma * gammaMultiplier,
                         simObj$initDelta * deltaMultiplier,
                         makeScenarioZeta(simObj$initZeta, zetaMultiplier) )

  # Other settings ...
  simObj$gearUseStrategyPA <- gearUseStrategyPA
  simObj$gearUseStrategyPO <- gearUseStrategyPO
  simObj <- setSimUseSDMs(simObj, useSDM, simObj$minPrevalence)

  # Return simulation object.
  return(simObj)
  
}
  
#-----------------------------------------------------------------------------------------

deSink <- function() {
  
  # Turn off diversion of output to file (usually lost ones caused by stopping the program).
  
  # Turn off message diversion.
  sink(type="message")

  # Turn off any other diversions.
  while ( sink.number() > 0) {
    sink()
    message("sink reversed")
  } 
  
}

#-----------------------------------------------------------------------------------------

