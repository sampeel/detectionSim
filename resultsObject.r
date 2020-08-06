#-----------------------------------------------------------------------------------------

initResults <- function(numRuns, numCoeffs, numBiases, namesSpecies, namesGears=NULL,
                        numSpecies=length(namesSpecies), numGears=length(namesGears)){
  
  # Initialise or reset the results list object.
  #
  # Arguments ...
  # numRuns:    the number of runs in the simulation
  # numSpecies: the number of species in the data
  # numCoeffs:  the number of coefficients in the intensity function including intercept 
  # numBiases:  the number of sample bias covariates (gives number of deltas)
  # numSpecies: the number of species in the data
  # numGears:   the number of gear types in the data (zero if this feature is not to be used)

  # Dummy names for dimensions.
  namesCoeffsLambda <- c("alpha", paste0("beta", 1:(numCoeffs-1)))
  namesCoeffsBias <- paste0("delta", 1:numBiases)
  namesCoeffsDet <- paste0("zeta", 2:numGears)                # NB: This is a factor!
  namesRuns <- paste0("run", 1:numRuns)
  
  # Create a skeleton list to store the results from each SDM.
  matEmpty <- matrix(data=NA, nrow=numSpecies, ncol=numRuns, 
                     dimnames = list(namesSpecies, namesRuns))
  beta <- array(data=NA, dim=c(numCoeffs, numSpecies, numRuns), 
                dimnames=list(namesCoeffsLambda, namesSpecies, namesRuns))
  zeta <- array(data=NA, dim=c(numGears-1, numSpecies, numRuns), 
                dimnames=list(namesCoeffsDet, namesSpecies, namesRuns))
  gamma <- matEmpty
  delta <- matrix(data=NA, nrow=numBiases, ncol=numRuns, 
                  dimnames=list(namesCoeffsBias, namesRuns))
  coeffs <- list(beta=beta, zeta=zeta, gamma=gamma, delta=delta)
  lstEmpty <- list(tryMe=FALSE, fit=NA, coeffs=coeffs, SE=coeffs)
  
  # Initialise warning and error recording.
  errorsEmpty <- as.data.frame(matrix(nrow=0,ncol=4), stringsAsFactors=FALSE)
  colnames(errorsEmpty) <- c("species", "SDM", "run", "msg")

  # Create the results list object.
  resLst <- list(numRuns=numRuns,            # Number of runs
                 numCoeffs=numCoeffs,        # Number of coefficients
                 numBiases=numBiases,        # Number of sample bias covariates
                 numSpecies=numSpecies,      # Number of species
                 namesSpecies=namesSpecies,  # a vector (numSpecies x 1) of species names
                 numGears=numGears,            
                 namesGears=namesGears,
                 numPresencesPA=matEmpty,    # a matrix (numSpecies x numRuns) of the number of presences per species per run in the simulated PA data
                 numPresencesPO=matEmpty,    # a matrix (numSpecies x numRuns) of the number of presences per species per run in the simulated PO data                
                 Count=lstEmpty,             # SDM of count/survey/abundance data 
                 Cloglog=lstEmpty,           # SDM of presence-absense data using complementary-log-log
                 PPM=lstEmpty,               # SDM of presence-only data using PPM function (Poisson point process model)
                 MsPP=lstEmpty,              # SDM of presence-absence and presence-only data using MultispeciesPP function (combined Poisson point process)
                 Gear=lstEmpty,              # SDM of presence-absence and presence-only data using gearGLM function (combined Poisson point process with gear)
                 errors=errorsEmpty,         # Errors: info about place of occurrence and error message.
                 warnings=errorsEmpty,       # Warnings: info about place of occurrence and warning message.
                 timeSetup = NULL,           # Time taken to run the setup part of the simulation (getting of data, creation of results storage, etc)
                 timeSim = NULL,             # Time taken to run the simulation for the given scenario
                 timeResPlots = NULL,        # Time taken to create the result plots for the given scenario (total time = timeSim + timeResPlots)
                 #
                 # The following values are worked out by the object's functions, no need to set.
                 #
                 validSDM=NULL,              # Vector of the SDM names from above items (not user set)
                 namesValidSDM=NULL,         # Vector of the SDM names for pretty plotting
                 numSDMs=0                   # Number of SDM to be used.
                )
  
  # Return value.
  return(resLst)
}

#-----------------------------------------------------------------------------------------

setResTryMe <- function(simLst, resLst) {
  
  # Sets the tryMe value for each SDM type from the appropriate value in the sim list object.
  # Also sets the validSDM vector (names must be the same as SDM item names in result object list).
  
  # Reset this value to nothing.
  resLst$validSDM <- NULL
  
  # Set tryMe values using useSDM values in the simulation object list.
  resLst$Count$tryMe <- simLst$useSDM.Count
  if ( resLst$Count$tryMe ) {
    resLst$validSDM <- c(resLst$validSDM, "Count")
    resLst$namesValidSDM <- c(resLst$namesValidSDM, "AB")
  } else {
    resLst$Count <- NULL
  }
  
  resLst$Cloglog$tryMe <- simLst$useSDM.Cloglog
  if ( resLst$Cloglog$tryMe ) {
    resLst$validSDM <- c(resLst$validSDM, "Cloglog")
    resLst$namesValidSDM <- c(resLst$namesValidSDM, "PA")
  } else {
    resLst$Cloglog <- NULL
  }
  
  resLst$PPM$tryMe <- simLst$useSDM.PPM
  if ( resLst$PPM$tryMe ) {
    resLst$validSDM <- c(resLst$validSDM, "PPM")
    resLst$namesValidSDM <- c(resLst$namesValidSDM, "PO")
  } else {
    resLst$PPM <- NULL
  }
  
  resLst$MsPP$tryMe <- simLst$useSDM.MsPP
  if ( resLst$MsPP$tryMe ) {
    resLst$validSDM <- c(resLst$validSDM, "MsPP")
    resLst$namesValidSDM <- c(resLst$namesValidSDM, "MsPP")
  } else {
    resLst$MsPP <- NULL
  }
  
  resLst$Gear$tryMe <- simLst$useSDM.Gear
  if ( resLst$Gear$tryMe ) {
    resLst$validSDM <- c(resLst$validSDM, "Gear")
    resLst$namesValidSDM <- c(resLst$namesValidSDM, "Gear")
  } else {
    resLst$Gear <- NULL
  }

  # Number of SDMs to be used.
  resLst$numSDMs <- length(resLst$validSDM)
  
  # Return value.
  return(resLst)
  
}
  
#-----------------------------------------------------------------------------------------

setResults <- function(resLstSDM, thisRun, betas, isCoeffs=TRUE, 
                       gammas=NULL, deltas=NULL, zetas=NULL) {
  
  # Sets the estimated coefficients OR the standard errors for the given SDM for thisRun 
  # and all species.
  #
  # Arguments ...
  # resLstSDM: the SDM part of the result list (e.g. resLst$Count or resLst$PPM)
  # thisRun:   the run number under which these estimated coefficients are to be stored.
  # betas:     a matrix with a column for each species and rows that are alpha, beta1, ..., 
  #            betaN, where N is the number of covariates.
  # isCoeffs:  store values as coefficients if TRUE and standard errors if FALSE.
  # gammas:    a matrix with numSpecies elements that are the estimates of the bias intercept.
  # deltas:    a matrix with numBiases elements that are the estimates of the sample bias 
  #            covariate coefficients.
  # zetas:     an array with numGears-1 rows and numSpecies columns that are the estimates of 
  #            the gear effect.

  # Coefficient estimated values or standard errors?
  if ( isCoeffs ) {
    thisItem <- "coeffs"
  } else {
    thisItem <- "SE"
  }
  
  # Get the names of the species that have coefficient estimates OR standard errors
  namesSpecies <- colnames(betas)
  
  # Store the values in the given result list.
  resLstSDM[[thisItem]]$beta[ ,namesSpecies,thisRun] <- betas
  if ( ! is.null(gammas) ) resLstSDM[[thisItem]]$gamma[namesSpecies,thisRun] <- gammas
  if ( ! is.null(deltas) ) resLstSDM[[thisItem]]$delta[ ,thisRun] <- deltas
  if ( ! is.null(zetas) ) resLstSDM[[thisItem]]$zeta[ ,namesSpecies,thisRun] <- zetas

  # Return value.
  return(resLstSDM)
  
}

#-----------------------------------------------------------------------------------------

setResError <- function(resLst, run, SDM=NA, species=NA, msg="", isWarning=FALSE) {

  # Sets the error stuff for the given run (and SDM and species, if given).
  # 
  # Arguments ...
  # resLst:  result object
  # run:     scalar integer giving the run that has had an error.
  # SDM:     character scalar or vector giving the method that has had an error (NA if error  
  #          has occurred for all methods).  If a vector, must have same length as errMsg.
  # species: integer scalar or vector giving the species that has had an error (NA if error
  #          has occurred for all species).  If a vector, must have same length as errMsg.
  # msg:     character scalar or vector that contains the error message for each error.
  # isWarning: if true, set the warning stuff instead.

  # How many errors are to be recorded?
  numErrors <- length(msg)
  
  # Make sure run is a scalar (and then make a vector of same value, if necessary).
  if ( length(run) > 1) {
    stop("A vector of run numbers is not accepted in function setResError.")
  } else {
    run <- rep(run, numErrors)
  }
  
  # Check species is the same size as errMsg.
  if ( length(species) != numErrors ) {
    stop("Error messages and species identifier are not the same length in error setting.")
  }
  
  # If SDM is a vector, check it is also the same size.
  if ( length(SDM) > 1 ) {
    if ( length(SDM) != numErrors ) {
      stop("Error messages and SDM identifier are not the same length in error setting.")
    }
  } else {
    SDM <- rep(SDM, numErrors)
  }
    
  # Add this error to the error list.
  newRow <- data.frame(species=species, SDM=SDM, run=run, msg=msg, stringsAsFactors = FALSE)
  if ( isWarning ) {
    resLst$warnings <- rbind(resLst$warnings, newRow)
  } else {
    resLst$errors <- rbind(resLst$errors, newRow)
  }
  
  # Return value.
  return(resLst)
  
}

#-----------------------------------------------------------------------------------------

isSimError <- function(errors, run, SDM, species) {
  
  # Did the given run, SDM and species have an error during the simulation?  Search the 
  # list of errors for a matching identification (combination of run, SDM and species).
  # Remember that "NA" is recorded for errors occurring in all SDMs and/or all species.
  # Thus, any value of SDM will "match" any NA, ditto for species.
  
  # If there are no errors.
  if ( length(errors$run) == 0 ) return(FALSE)
  
  # Search for the given run in the list of errors.
  indRun <- which(errors$run == run)
  if ( length(indRun) == 0 ) return(FALSE)
  
  # Cut to possibles and search for the given SDM.
  possErrors <- errors[indRun,]
  indSDM <- which((possErrors$SDM == SDM) || ( is.na(possErrors$SDM)))
  if ( length(indSDM) == 0 ) return(FALSE)
  
  # Cut to possibles and search for the given species.
  possErrors <- possErrors[indSDM]
  indSpecies <- which((possErrors$species == species) || ( is.na(possErrors$species)))
  if ( length(indSpecies) == 0 ) return(FALSE)
  
  # There was an error during the simulation.
  return(TRUE)
  
}

#-----------------------------------------------------------------------------------------

setResOneToMany <- function(resLst, res1, whichRun) {
  
  # Assigns the results from one simulation (res1) to the full result list (resLst) at the 
  # whichRun position.
  # Returns the result list with the added results.
  
  # Error stuff.
  numErrors <- dim(res1$errors)[1] 
  if ( numErrors > 0 ) {
    resLst <- setResError(resLst, whichRun, res1$errors$SDM, 
                          res1$errors$species, res1$errors$msg)    
  }
  
  # Warning stuff.
  numWarnings <- dim(res1$warnings)[1]
  if ( numWarnings > 0 ) {
    resLst <- setResError(resLst, whichRun, res1$warnings$SDM, 
                          res1$warnings$species, res1$warnings$msg, TRUE)    
  }

  # Numbers of presences.
  resLst$numPresencesPA[ ,whichRun] <- res1$numPresencesPA[ ,1]
  resLst$numPresencesPO[ ,whichRun] <- res1$numPresencesPO[ ,1]
  
  # How many SDMs where used?
  #numValidSDMs <- length(resLst$validSDM)
  
  # Which SDMs where used?
  for ( nameSDM in resLst$validSDM ) {
    # Which SDM is this?
    #nameSDM <- resLst$validSDM[j]
    
    # Coefficients ...
    resLst[[nameSDM]] <- setResults(resLst[[nameSDM]], 
                                    whichRun, 
                                    res1[[nameSDM]]$coeffs$beta,
                                    TRUE, 
                                    res1[[nameSDM]]$coeffs$gamma, 
                                    res1[[nameSDM]]$coeffs$delta,
                                    res1[[nameSDM]]$coeffs$zeta)

    # Standard Errors ...
    resLst[[nameSDM]] <- setResults(resLst[[nameSDM]], 
                                    whichRun, 
                                    res1[[nameSDM]]$SE$beta,
                                    FALSE, 
                                    res1[[nameSDM]]$SE$gamma, 
                                    res1[[nameSDM]]$SE$delta,
                                    res1[[nameSDM]]$SE$zeta)
  }

  # Return value.
  return(resLst)

}

#-----------------------------------------------------------------------------------------

setResTimes <- function(resLst, startTime, endTime, whichPart=0) {
  
  # NB: whichPart = 1 for setup, 
  #               = 2 for sim timing, or 
  #               = 3 for result plots timing.
  
  # Calculate time taken to run the part of the sim for this scenario.
  timeTaken <- format(endTime - startTime)
  
  # Which part of the whole code is this time for?
  if ( whichPart == 1 ) {
    resLst$timeSetup <- timeTaken
  } else if ( whichPart == 2 ) {
    resLst$timeSim <- timeTaken
  } else if ( whichPart == 3) {
    resLst$timeResPlots <- timeTaken
  } else {
    # Do nothing?
    warning("Unrecognised value given in 'whichPart' argument.")
  }
  
  # Return value.
  return(resLst)
  
}

#-----------------------------------------------------------------------------------------

setResNumPresences <- function(resLst, PAsp, POsp, thisRun=1){

  # Sets the number of presences for thisRun and each species.
  # Note that the number of absences is the number of surveys minus the number of presences.
  # for each run.
  #
  # Arguments ...
  # resLst:  result object
  # PAsp:    just the species columns of a PA data.frame.
  # POsp:    just the species name column of a PO data.frame.
  # thisRun: the run number that applies to this data.
  
  # Calculate the number of things.
  numSurveys <- dim(PAsp)[1]
  numSpecies <- dim(PAsp)[2]
  if ( numSpecies != dim(resLst$numPresencesPA)[1] ) {
    stop("Result dimensions don't match the number of species given in PA data.")  
  }
  
  # Calculate the number of presences for each species from the PA data.
  resLst$numPresencesPA[ ,thisRun] <- apply(PAsp, 2, sum)
  
  # Calculate the number of presences for each species from the PO data.
  numPresencesPerSpecies <- table(POsp)
  speciesNames <- names(numPresencesPerSpecies) 
  resLst$numPresencesPO[speciesNames,thisRun] <- numPresencesPerSpecies
  
  # Return value.
  return(resLst)
  
}

#-----------------------------------------------------------------------------------------

makeCoeffBiasVar <- function(resLstSDM, coeffEstimates, coeffs) {
  
  # Calculate the bias and variance of the estimated coefficients where 
  #                mse = variance + bias^2
  # (see wiki "Mean Square Error")
  #
  # Arguments ...
  # resLstSDM:      the SDM part of the result list (e.g. resLst$Count or resLst$PPM)
  # coeffEstimates: the estimated value of the coefficients (numCoeffs x numSpecies x numRuns)
  # coeffs:         the true value of the coefficients (numCoeffs x numSpecies)
  
  # Numbers and names of things.
  numCoeffs <- dim(coeffEstimates)[1]
  namesCoeffs <- dimnames(coeffEstimates)[[1]]
  numSpecies <- dim(coeffEstimates)[2]
  namesSpecies <- dimnames(coeffEstimates)[[2]]
  numRuns <- dim(coeffEstimates)[3]
  
  # Check dimensions are the same.
  if ( (dim(coeffs)[1] != numCoeffs) || (dim(coeffs)[2] != numSpecies) ) 
    stop("Dimensions of arguments don't match.")
  
  # Initialise return values.
  bias <- matrix(nrow=numCoeffs, ncol=numSpecies, dimnames=list(namesCoeffs, namesSpecies))
  bias <- as.data.frame(bias, stringsAsFactors=FALSE)
  var <- bias

  # Calculate bias for each coeff x species combination.
  for ( coeff in namesCoeffs ) {
    for ( species in namesSpecies ) {
      # Which runs were successful for this species and this coefficient?
      indSuccessfulRuns <- which(! is.na(coeffEstimates[coeff,species, ]))
      resLstSDM$numCoeffsEst[coeff,species] <- length(indSuccessfulRuns)
      
      if ( length(indSuccessfulRuns) > 0 ) {
        # Mean of the estimated values.
        meanEstimates <- mean(coeffEstimates[coeff,species, indSuccessfulRuns])
        
        # Calculate bias.
        bias[coeff,species] <- meanEstimates - coeffs[coeff,species,indSuccessfulRuns] 
        
        # Calculate variance.
        var[coeff,species] <- mean(((coeffEstimates[coeff,species,indSuccessfulRuns] - meanEstimates)^2))
      }
    }
  }

  # Return value.
  resLstSDM$bias <- bias
  resLstSDM$var <- var
  return(resLstSDM)
}  
  
#-----------------------------------------------------------------------------------------

makeNumSuccessfulRuns <- function(resLstSDM, betaEstimates) {
  
  # The number of successful runs (i.e. estimates of coefficients) for this SDM.
  # Only uses the first coefficient (assumes that if it works for one then will work for all)!
  # Returned as part of the resLstSDM list.
  # 
  # Arguments ...
  # resLstSDM:      the SDM part of the result list (e.g. resLst$Count or resLst$PPM)
  # betaEstimates:  the estimated value of the coefficients (numCoeffs x numSpecies x numRuns)

  # Numbers and names of things.
  numSpecies <- dim(betaEstimates)[2]
  namesSpecies <- dimnames(betaEstimates)[[2]]
  numRuns <- dim(betaEstimates)[3]
  
  # Initialise return values.
  resLstSDM$numSuccessfulRuns <- matrix(nrow=1, ncol=numSpecies, 
                                   dimnames=list(NULL, namesSpecies))

  # Calculate bias for each coeff x species combination.
  for ( species in namesSpecies ) {
    # Which runs were successful for this species and the first coefficient.
    indSuccessfulRuns <- which(! is.na(betaEstimates[1,species, ]))
    resLstSDM$numSuccessfulRuns[1,species] <- length(indSuccessfulRuns)
  }

  # Return value.
  return(resLstSDM)
}  

#-----------------------------------------------------------------------------------------
