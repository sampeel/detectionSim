runSim <- function(i, simObj, domainObj, cellObj, plotObj, BG) {
  
  # Run one simulation with the given objects. "i" is the number of this run.
  # NB: results are organised so that this function can be used in parallel.
  
  # Catch if an error occurs in this run and report (rather than stopping everything!)
  tryCatchOut <- tryCatch(
    {
      message(paste(i, "Starting sim run", Sys.time()))
      
      # Initialise data dump values.
      fileDump <- paste0("DataDump-run", i)
      fileDump <- makeFileName(fileDump, simObj$dataDumpDir, "RData")
      cellObj <- resetSimDataCells(cellObj)
      surveyObj <- resetSimDataSurveys(surveyObj)

      # Initialise the result list object for one run.
      resLst <- initResults(1, simObj$numCoeffs, simObj$numBiases, cellObj$namesSpecies,
                            cellObj$namesGears)
      resLst <- setResTryMe(simObj, resLst)
      
      # Make the number of individuals per cell for each species.
      message(paste(i, "Begin makeNumIndivids"))
      cellObj <- makeNumIndivids(cellObj, i)
      if ( plotObj$countPerCell && ( i == 1 || plotObj$checkPlotsAllRuns ) )
        plotCellValues(plotObj, domainObj$mask, cellObj$cells, cellObj$N, i, "CheckCellNumIndivids",
                       titleTxts = "Count per cell of simulated individuals for ")
      message(paste(i, "End makeNumIndivids"))
      
      
      # Make the presence-only data (no need for mask as individuals are already in the domain!)
      message(paste(i, "Begin makePOObservations"))
      cellObj <- makePOObservations(cellObj, i, simObj$numGears, simObj$namesGears,
                                    simObj$gearUseStrategyPO, cellObj$covars[ ,1])
      if ( plotObj$POData  && ( i == 1 || plotObj$checkPlotsAllRuns ) )
        plotPOPoints(plotObj, domainObj$mask, cellObj, cellObj$PO, i, "CheckPOPoints")
      message(paste(i, "End makePOObservations"))
      
      
      # Make the survey locations.
      message(paste(i, "Begin makeSampleLocations"))
      surveyObj <- makeSampleLocations(domainObj$mask, cellObj, simObj$numSamples, 
                                       simObj$minSampleArea, simObj$maxSampleArea, simObj$numSurveys,
                                       simObj$bias.formula, c(0,simObj$initDelta), cellObj$biases,
                                       simObj$widthSurveys)
      if ( plotObj$samplesPerCell && ( i == 1 || plotObj$checkPlotsAllRuns ))
        plotSampleAreas(plotObj, domainObj$mask, cellObj, surveyObj, "CheckSampleCells")

      # Assign gear type used at each sample location.
      surveyObj <- assignGearUsedSamples(surveyObj,
                                         simObj$numGears, 
                                         simObj$namesGears,
                                         simObj$gearUseStrategyPA, 
                                         cellObj$covars[surveyObj$rowsInCells,"data1"],
                                         range(cellObj$covars[ ,"data1"]))
      
      if ( plotObj$samplesPerCell && ( i == 1 || plotObj$checkPlotsAllRuns ))
        plotSampleGears(plotObj, domainObj$mask, cellObj, surveyObj, "CheckSampleGears")
      message(paste(i, "End makeSampleLocations"))
      
      # Make the count data and the presence-absence data
      message(paste(i, "Begin makeNSampled"))
      surveyObj <- makeNSampled(surveyObj, cellObj$N, cellObj$trueD, cellObj$areaCell)
      if ( plotObj$countData && ( i == 1 || plotObj$checkPlotsAllRuns ) )
        plotCellCount(domainObj, plotObj, surveyObj$cells, surveyObj$N, 
                      cellObj$namesSpecies, i, "CheckCountData")
      if ( plotObj$PAData && ( i == 1 || plotObj$checkPlotsAllRuns ) )
        plotSurveyPA(domainObj, plotObj, surveyObj$xy, surveyObj$Y, 
                     cellObj$namesSpecies, i, "CheckPAData")
      message(paste(i, "End makeNSampled"))

      # Record the number of presences per species (from both types of data)      
      resLst <- setResNumPresences(resLst, surveyObj$Y[ ,cellObj$namesSpecies], 
                                   cellObj$PO$species)

      # Restrict SDM usage to just those species with enough data.
      prevLst <- checkPrevalence(cellObj$namesSpecies, surveyObj$Y, cellObj$PO, 
                                 simObj$minPrevalence, cellObj$P)
      if ( prevLst$isErrors ) {
        resLst <- setResError(resLst, 1, species = prevLst$errors$species,
                              msg = prevLst$errors$msg)
      }
      if ( prevLst$isWarnings ) {
        resLst <- setResError(resLst, 1, species = prevLst$warnings$species, 
                              msg = prevLst$warnings$msg, isWarning=TRUE)
      }
      
      # Run the SDMs to estimate the true intensity (and the MSE of the estimated coefficients).
      # if ( simObj$useSDM.Count ) {
      #   # glm with count data.
      #   message(paste(i, "Begin AB SDM"))
      #   tryRes <- try(runSDM.Count(simObj$lambda.formula, surveyObj$cells, 
      #                              cellObj$covars[surveyObj$rowsInCells, ],
      #                              surveyObj$N[ ,prevLst$namesSpeciesPA], surveyObj$areas, 
      #                              simObj$gear.formula, surveyObj$gears))
      #                             
      #   # tryRes <- try(runSDM.Count(simObj$lambda.formula, Count$cell,
      #   #                            cellObj$covars[surveyObj$rowsInCells, ],
      #   #                            Count[ ,prevLst$namesSpeciesPA], Count$area))
      #   if ( inherits(tryRes, "try-error") ) {
      #     # An error occured that stopped the function.
      #     resLst <- setResError(resLst, 1, "AB", msg = tryRes[1])
      #   } else {
      #     # Some errors may have occured inside the function that did not require it to stop.
      #     if ( tryRes$isError ) {
      #       resLst <- setResError(resLst, 1, "AB", tryRes$errors$species, tryRes$errors$msg)
      #     }
      #     
      #     # Set the estimated coefficients.
      #     resLst$Count <- setResEstCoeffs(resLst$Count, tryRes$estCoeffs, 1, tryRes$stdErrors)
      #   }
      #   message(paste(i, "End AB SDM"))
      # }
      
      if ( simObj$useSDM.Cloglog ) {
        # glm with PA data.
        message(paste(i, "Begin PA SDM"))
        tryRes <- try(runSDM.Cloglog(simObj$lambda.formula, surveyObj$cells, 
                                   cellObj$covars[surveyObj$rowsInCells, ],
                                   surveyObj$Y[ ,prevLst$namesSpeciesPA], surveyObj$areas, 
                                   cellObj$namesSpecies, simObj$gear.formula, surveyObj$gears))
        # tryRes <- try(runSDM.Cloglog(simObj$lambda.formula, PA$cell,
        #                              cellObj$covars[surveyObj$rowsInCells, ],
        #                              PA[ ,prevLst$namesSpeciesPA], PA$area))
        if ( inherits(tryRes, "try-error") ) {
          # An error occured that stopped the function.
          resLst <- setResError(resLst, 1, "PA", msg = tryRes[1])
        } else {
          # Some errors may have occured inside the function that did not require it to stop.
          if ( tryRes$isError ) {
            resLst <- setResError(resLst, 1, "PA", tryRes$errors$species,
                                  tryRes$errors$msg)
          }
          
          # Set the estimated coefficients and standard errors.
          resLst$Cloglog <- setResults(resLst$Cloglog, 1, tryRes$coeffs$beta, TRUE, 
                                       tryRes$coeffs$gamma, tryRes$coeffs$delta, 
                                       tryRes$coeffs$zeta)
          resLst$Cloglog <- setResults(resLst$Cloglog, 1, tryRes$SE$beta, FALSE, 
                                       tryRes$SE$gamma, tryRes$SE$delta, 
                                       tryRes$SE$zeta)
        }
        message(paste(i, "End PA SDM"))
      }
      
      # if ( simObj$useSDM.PPM ) {
      #   # glm with PO and BG data.
      #   message(paste(i, "Begin PO SDM"))
      #   tryRes <- try(runSDM.PPM(simObj$lambda.formula, cellObj$PO[prevLst$whichPORows, ], 
      #                            BG, domainObj, cellObj))
      #   
      #   if ( inherits(tryRes, "try-error") ) {
      #     # An error occured that stopped the function.
      #     resLst <- setResError(resLst, 1, "PO", msg = tryRes[1])
      #   } else {
      #     # Some errors may have occured inside the function that did not require it to stop.
      #     if ( tryRes$isError ) {
      #       resLst <- setResError(resLst, 1, "PO", tryRes$errors$species,
      #                             tryRes$errors$msg)
      #     }
      #     
      #     # Set the estimated coefficients.
      #     resLst$PPM <- setResEstCoeffs(resLst$PPM, tryRes$estCoeffs, 1, tryRes$stdErrors)
      #   }
      #   message(paste(i, "End PO SDM"))
      # }
      
      if ( simObj$useSDM.MsPP ) {
        # glm with PA, PO and BG data.
        message(paste(i, "Begin MsPP SDM"))
        tryRes <- try(runSDM.MsPP(simObj$lambda.formula, simObj$bias.formula, surveyObj$cells, 
                                  surveyObj$areas, surveyObj$Y[ ,prevLst$namesSpeciesPA], 
                                  cellObj$PO[prevLst$whichPORows, ], BG, cellObj))
        
        if ( inherits(tryRes, "try-error") ) {
          # An error occured that stopped the function.
          resLst <- setResError(resLst, 1, "MsPP", msg = tryRes[1])
          
        } else {
          # Some errors may have occured inside the function that did not require it to stop.
          if ( tryRes$isError ) {
            resLst <- setResError(resLst, 1, "MsPP", species = tryRes$errors$species,
                                  msg = tryRes$errors$msg)
          }
          
          # Set the estimated coefficients and standard errors.
          resLst$MsPP <- setResults(resLst$MsPP, 1, tryRes$coeffs$beta, TRUE, 
                                       tryRes$coeffs$gamma, tryRes$coeffs$delta, 
                                       tryRes$coeffs$zeta)
          resLst$MsPP <- setResults(resLst$MsPP, 1, tryRes$SE$beta, FALSE, 
                                       tryRes$SE$gamma, tryRes$SE$delta, 
                                       tryRes$SE$zeta)
        }
        message(paste(i, "End MsPP SDM"))
      }
      
      if ( simObj$useSDM.Gear ) {  
        # glm with PA, PO and gear data.
        message(paste(i, "Begin Gear SDM"))
        tryRes <- try(runSDM.Gear(simObj$lambda.formula, 
                                  simObj$bias.formula,
                                  simObj$gear.formula,
                                  surveyObj$rowsInCells,
                                  surveyObj$areas, 
                                  surveyObj$Y[ ,prevLst$namesSpeciesPA],
                                  cellObj$P[ ,prevLst$namesSpeciesPO],  
                                  surveyObj$gears,
                                  cellObj))
          
        if ( inherits(tryRes, "try-error") ) {
          # An error occured that stopped the function.
          resLst <- setResError(resLst, 1, "Gear", msg = tryRes[1])

        } else {
          # Some errors may have occured inside the function that did not require it to stop.
          if ( tryRes$isError ) {
            resLst <- setResError(resLst, 1, "Gear", species = tryRes$errors$species,
                                  msg = tryRes$errors$msg)
          }

          # Set the estimated coefficients and standard errors.
          resLst$Gear <- setResults(resLst$Gear, 1, tryRes$coeffs$beta, TRUE, 
                                       tryRes$coeffs$gamma, tryRes$coeffs$delta, 
                                       tryRes$coeffs$zeta)
          resLst$Gear <- setResults(resLst$Gear, 1, tryRes$SE$beta, FALSE, 
                                       tryRes$SE$gamma, tryRes$SE$delta, 
                                       tryRes$SE$zeta)
        }

        message(paste(i, "End Gear SDM"))
      }
      
      # Check if there has been an error in any of the methods.
      if ( dim(resLst$errors)[1] > 0 ) {
        # Dump data necessary for an SDM run (for analysis if there is an error).
        save(i, BG, cellObj, surveyObj, file = fileDump)
        
      } else if ( i == 1 || i == simObj$numRuns) {
        # This ensures that there are at least two data dumps to be used for problem analysis.
        # Used a different file name but it is the same data as the error dump.
        # Different file name so that they are not mistaken as an error having occured!
        fileExamp <- paste0("DataExample-run", i)
        fileExamp <- makeFileName(fileExamp, simObj$dataDumpDir, "RData")
        save(i, BG, cellObj, surveyObj, file = fileExamp)
      }
      
      # Return results.
      message(paste(i,"Finished sim run", Sys.time()))
      return(resLst)
    },
    error = function(cond) {
      # What to do if there was an error in other bits of code.
      msg <- conditionMessage(cond)
      msg <- paste0("ERROR on run ", i, ": ", msg)
      resLst <- setResError(resLst, 1, msg = msg)
      
      # Dump data necessary for an SDM run (for analysis if there is an error).
      save(i, BG, cellObj, surveyObj, file = fileDump)
      
      # Return results.
      return(resLst)
    }
  )
  
  # Any uncaught errors?
  if ( inherits(tryCatchOut, "try-error")) {
    resLst <- setResError(resLst, 1, msg = tryRes$errors$msg)
    return(resLst)
  }
  
  # Return value.
  return(tryCatchOut)
  
}

#-----------------------------------------------------------------------------------------

getMaxLambda <- function(domain, res, myFormula, coeffs, covars) {

  # Get the maximum value of lambda (the intensity function) for each set of coefficients.
  # Uses rasters to calculate the value of lambda at the specified resolution.  The value
  # of the resolution will dictate how "accurate" the result is.  Assumes there is a
  # function called "lambda" that has the arguments (x,y,coeffs,covars) in the workspace.
  #
  # Arguments ...
  # domain:    a domain object on which the lambda values are to be calculated.
  # res:       the resolution at which the lambda values are to be calculated.
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    the vector of the coefficients of the lambda function.  If coeffs is a
  #            matrix, the coefficients of a lambda function are assumed to be in each
  #            column and the number of columns is the number of lambda functions whose
  #            maximum is to be found.
  # covars:    the RasterStack of covariates for the lambda function.  Must have values at
  #            all points in the domain (combination of extent and mask!)

  # The number of intensity function coefficients provided.
  numFuncs <- dim(coeffs)[2]

  # Initialise the return value.
  maxLambda <- drop(matrix(nrow=numFuncs, ncol=1))

  # Find maximum for each intensity function ...
  for (k in 1:numFuncs) {
    # This intensity function's coefficients.
    thisCoeffs <- coeffs[ ,k]

    # Make a RasterLayer of the values of lambda.
    rsLambda <- rasterFromFunc(domain$ext, res, domain$mask, domain$maskValue, lambda,
                               myFormula, thisCoeffs, covars)

    # The maximum intensity function value is ...
    maxLambda[k] <- maxValue(rsLambda)
  }

  # Return value.
  return(maxLambda)

}

#-----------------------------------------------------------------------------------------

getMaxLambda2 <- function(domain, res, myFormula, coeffs, covars,
                          numCoresUse=max(1, detectCores() - 1)) {

  # Get the maximum value of lambda (the intensity function) for each set of coefficients.
  # Uses rasters to calculate the value of lambda at the specified resolution.  The value
  # of the resolution will dictate how "accurate" the result is.  Assumes there is a
  # function called "lambda" that has the arguments (x,y,coeffs,covars) in the workspace.
  #
  # Arguments ...
  # domain:    a domain object on which the lambda values are to be calculated.
  # res:       the resolution at which the lambda values are to be calculated.
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    the vector of the coefficients of the lambda function.  If coeffs is a
  #            matrix, the coefficients of a lambda function are assumed to be in each
  #            column and the number of columns is the number of lambda functions whose
  #            maximum is to be found.
  # covars:    the RasterStack of covariates for the lambda function.  Must have values at
  #            all points in the domain (combination of extent and mask!)
  # numCoresUse: number of cores to use in the mclapply function that paralellises the code.

  # The number of intensity function coefficients provided.
  numFuncs <- dim(coeffs)[2]

  ## How many cores are there (don't use all of them)
  #numCoresUse <- max(1, detectCores() - 1)

  # Find maximum for each intensity function ...
  maxLambda <- mclapply(1:numFuncs, getOneMaxLambda, domain, res, myFormula, coeffs, covars,
                        mc.cores=numCoresUse)
  maxLambda <- unlist(maxLambda)

  # Return value.
  return(maxLambda)

}

#-----------------------------------------------------------------------------------------

getOneMaxLambda <- function(k, domain, res, myFormula, coeffs, covars){

  # This intensity function's coefficients.
  thisCoeffs <- coeffs[ ,k]

  # Make a RasterLayer of the values of lambda.
  rsLambda <- rasterFromFunc(domain$ext, res, domain$mask, domain$maskValue, lambda,
                             myFormula, thisCoeffs, covars)

  # The maximum intensity function value is ...
  maxLambda <- maxValue(rsLambda)

  # Return value.
  return(maxLambda)

}

#-----------------------------------------------------------------------------------------

lambda.xy <- function (x, y, myFormula, coeffs, covars, lnLambda=FALSE) {

  # Intensity function.  For inhomogeneous Poisson point distribution and model of this as
  #
  #        ln(lambda(x,y)) = formula with coeffs and covars(x,y)
  #
  # Returns a vector containing the value of lambda(x,y) for all points (x,y).
  # Assumes (x,y) will have valid values in covars (i.e. only use (x,y) that do!!!)
  #
  # ARGUMENTS ...
  # x:         a vector of longitude values of the points
  # y:         a vector of latitude values of the points
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    a vector containing the coefficients for the intensity formula.
  # covars:    a raster stack containing numEnvirs layers and values for all x and y.
  # lnLambda:  return ln(lambda(x,y)) for TRUE or lambda(x,y) for FALSE

  # Prepare the formula by removing any response variable.
  myTerms <- terms(myFormula, keep.order = TRUE)
  myTerms <- delete.response(myTerms)

  # Get the relevant data.
  xy <- cbind(x,y)
  vals <- davesExtract.v3(covars, xy)   # NB: one column each for each raster layer in stack
  vals <- as.data.frame(cbind(xy, vals), stringsAsFactors=FALSE)
  mat <- as.matrix(model.frame(myTerms, vals))

  # Is there an intercept term?
  if ( attr(myTerms, "intercept") ) {
    ones <- rep(1, length.out=length(x))
    mat <- cbind(ones, mat)
  }

  # Evaluate the function.
  lnRes <- mat %*% coeffs
  if ( lnLambda ) {
    return(drop(lnRes))
  } else {
    return(drop(exp(lnRes)))
  }

}

#-----------------------------------------------------------------------------------------

lambda.cell <- function (myFormula, coeffs, covars, lnLambda=FALSE) {

  # Intensity function.  For inhomogeneous Poisson point distribution and model of this as
  #
  #        ln(lambda(cell)) = formula with coeffs and covars(cell)
  #
  # Returns a vector containing the value of lambda(cell) for all cells with covariate 
  # values (i.e. all rows in covars data.frame).
  #
  # ARGUMENTS ...
  # myFormula: a formula object that contains how to combine the information for an intensity function.
  # coeffs:    a vector containing the coefficients for the intensity formula.
  # covars:    a data.frame containing columns of covariate data with the values in each
  #            cell forming the rows.  Only contains information for cells in the domain.
  #            NB: covariate columns need to be named the same as in the formula!
  # lnLambda:  return ln(lambda(cell)) for TRUE or lambda(cell) for FALSE

  # Prepare the formula by removing any response variable.
  myTerms <- terms(myFormula, keep.order = TRUE)
  myTerms <- delete.response(myTerms)

  # Get the relevant data.
  mat <- as.matrix(model.frame(myTerms, covars))

  # Is there an intercept term?
  if ( attr(myTerms, "intercept") ) {
    ones <- rep(1, length.out=nrow(covars))
    mat <- cbind(ones, mat)
  }

  # Evaluate the function.
  lnRes <- mat %*% coeffs
  if ( lnLambda ) {
    return(drop(lnRes))
  } else {
    return(drop(exp(lnRes)))
  }

}

#-----------------------------------------------------------------------------------------

makeSimBG <- function(numPoints, domainOwin, domainMask, numReps=5) {

  # Make the quasi-random background points for use in the multispeciesPP method.
  # Returns a data frame with three columns, called x, y and cell.
  #
  # Arguments ...
  # numPoints: the number of background points required.
  # domainOwin: a domain observation window (see spatstat).
  # domainMask: a raster layer that indicates which cells are contained within the domain.
  # numReps:   the number of repeats of the number of the number of points.
  #            This combats the fact that rQuasi returns the same points each time by
  #            choosing a different set of numPoints from the numReps*numPoints.
  #            Set numReps = 1 if you do want the same points!

  # Check there has been a number specified.
  if ( numPoints < 1 ) {
    stop("Set the number of required background points before making background data.")
  }

  # Generate quasi-random points (lots of points as owin affects how many are made)
  numReps <- 5
  BG <- as.data.frame(rQuasi(numPoints*numReps, domainOwin, type="Halton"))

  # Choose a random starting point (as the same sequence of points every time).
  numQuasiPoints <- dim(BG)[1]
  indStart <- sample(1:(numQuasiPoints-numPoints+1), size=1)
  indEnd <- indStart + numPoints - 1
  BG <- BG[indStart:indEnd, ]

  # Extract cell numbers from points?
  BG$cell <- cellFromXY(domainMask, BG)

  # Return value.
  return(BG)

}

#-----------------------------------------------------------------------------------------

runSDM.Count <- function(lambda.formula, sampleCells, covars, counts, sampleAreas,
                         gear.formula = as.formula(~ factor(gear)), sampleGears = NULL) {
  
  # Run a species distribution model for count/abundance data.  Use a Poisson glm with log link.
  # Returns the estimated coefficients for each species' intensity function (lambda.formula).
  #
  # ARGUMENTS ...
  # lambda.formula: a formula object that contains how to combine the information for an
  #                 intensity function.
  # sampleCells:    a vector of cell numbers indicating where the PA data was collected 
  #                 within the domain (i.e. which cells contain the sample areas)
  # covars:         a data.frame with covariate values in these cells for each covariate
  #                 (numSamples x numCovars)
  # counts:         matrix of the number of individuals counted in each survey (same order 
  #                 as cells) rows are for each survey site and columns are for each 
  #                 species (numSamples x numSpecies).
  # sampleAreas:    vector of the area of the samples defined by cell number (numSamples)
  # gear.formula:   a formula object that contains the gear specific component of the 
  #                 formula (use "gear" in this formula).  Combined with intensity
  #                 formula if sampleGears data is given.  
  # sampleGears:    a vector containing an indicator for which gear type was used in each
  #                 sample (numSamples).  If this is NULL, then this data is not available.

  # Add the gear factor and area offset to the formula and remove the response variable (if present)
  if ( ! is.null(sampleGears) ) {
    # Add gear formula.
    lambda.terms <- attr(terms(lambda.formula),"term.labels")
    gear.terms <- attr(terms(gear.formula),"term.labels")
    all.terms <- c(lambda.terms, gear.terms)
    sdm.formula <- reformulate(all.terms, response = NULL, intercept=TRUE)
  } else {
    sdm.formula <- lambda.formula
  }
  sdm.formula <- update.formula(sdm.formula, "N ~ . + offset(log(area))")

  # Put the data together into one structure (except for species PA)
  dfData <- data.frame(cell = sampleCells,
                       gear = sampleGears,
                       area = sampleAreas,
                       covars,
                       N = counts[ ,1])
  
  # Numbers of things.
  namesSpecies <- names(counts)
  numSpecies <- length(namesSpecies)
  modMat <- model.matrix(sdm.formula, dfData)
  numCoeffs <- length(attr(modMat,"assign"))
  if ( ! is.null(sampleGears) ) {
    numGears <- length(unique(sampleGears))
  } else {
    numGears <- 0
  }
  
  # Initialise the return value.
  estCoeffs <- matrix(nrow=numCoeffs, ncol=numSpecies)
  colnames(estCoeffs) <- namesSpecies
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retVal <- list(estCoeffs=estCoeffs, stdErrors=estCoeffs, errors=errors, isErrors=FALSE,
                 glmFit=NULL)

  # Estimate the coefficients for each species' intensity function.
  for ( species in namesSpecies) {
    # Get this species' count data.
    dfData$N <- counts[ ,species]

    # Estimate the coefficients for the model.
    glmFit <- glm(sdm.formula, family=poisson(), data=dfData, control=list(maxit=100))

    # Did the glm converge to an answer?
    if ( ! glmFit$converged ) {
      retVal$isErrors <- TRUE
      newError <- data.frame(species=species, msg="AB SDM has not converged.", stringsAsFactors = FALSE)
      retVal$errors <- rbind(retVal$errors, newError)
      return(retVal)
    } else {
      # Save estimated coefficients
      retVal$estCoeffs[ ,species] <- glmFit$coefficients
      retVal$stdErrors[ ,species] <- summary(glmFit)$coefficients[ ,2]
      retVal$glmFit <- glmFit
    }
  }

  # Return results.
  return(retVal)

}

#-----------------------------------------------------------------------------------------

runSDM.Cloglog <- function(lambda.formula, sampleCells, covars, PA, sampleAreas, 
                           namesSpecies = colnames(PA),
                           gear.formula = as.formula(~ factor(gear)), sampleGears = NULL) {
  
  # Run a species distribution model for presence-absence data.  Use a binomial glm with
  # a complimentary log log link.  Returns the estimated coefficients from the glm of the
  # centred data (using centreCovars).  Assumes there is enough data (both presences and
  # absences; each species) for the glm to work.  Adds areas as on offset in the glm.
  #
  # ARGUMENTS ...
  # lambda.formula: a formula object that contains how to combine the information for an
  #                 intensity function.
  # sampleCells:    a vector of cell numbers indicating where the PA data was collected 
  #                 within the domain (i.e. which cells contain the sample areas)
  # covars:         a data.frame with covariate values in these cells for each covariate
  #                 (numSamples x numCovars)
  # PA:             matrix of whether or not individuals are present (1) or absent (0) in
  #                 each sample (same order as cells) rows are for each sample site and
  #                 columns are for each species (numSamples x numSpecies)
  # sampleAreas:    vector of the area of the samples defined by cell number (numSamples)
  # namesSpecies:   a vector of ALL the species names (as PA data could be cropped to just
  #                 those that reach a minimum prevalence).
  # gear.formula:   a formula object that contains the gear specific component of the 
  #                 formula (not a valid formula on its own???).  Combined with intensity
  #                 formula only if sampleGears data is given.
  # sampleGears:    a vector containing an integer indicator for which gear type was used in each
  #                 sample (numSamples).  If this is NULL, then this data is not available.

  # Initialise return value.
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retVal <- list(errors=errors, 
                 isErrors=FALSE,
                 glmFit=NULL)

  # Numbers of things.
  numSpecies <- length(namesSpecies)
  numCoeffs <- numCoefficients(lambda.formula)
  if ( ! is.null(sampleGears) ) {
    numGears <- length(unique(sampleGears))
    namesZetas <- paste0("zeta", 1:numGears)
  } else {
    numGears <- 0
    namesZetas <- NULL
  }
  numSamples <- length(sampleCells)
  
  # Check species names in PA data match that in namesSpecies argument.
  namesSpeciesPA <- colnames(PA)
  if ( ! all(namesSpeciesPA %in% namesSpecies)) {
    retVal$isErrors <- TRUE
    newError <- data.frame(species=species, 
                           msg="Some species names within the PA data are different to the given list of species names.", 
                           stringsAsFactors = FALSE)
    retVal$errors <- rbind(retVal$errors, newError)
    return(retVal)
  }
  
  # Add the gear factor and area offset to the formula and remove the response variable (if present)
  if ( ! is.null(sampleGears) ) {
    # Add gear formula.
    lambda.terms <- attr(terms(lambda.formula),"term.labels")
    gear.terms <- attr(terms(gear.formula),"term.labels")
    all.terms <- c(lambda.terms, gear.terms)
    sdm.formula <- reformulate(all.terms, response = NULL, intercept=TRUE)
  } else {
    sdm.formula <- lambda.formula
  }
  sdm.formula <- update.formula(sdm.formula, "Y ~ . + offset(log(area))")
  
  # Create a skeleton list to store the results from each SDM.
  namesCoeffs <- c("alpha",paste0("beta", 1:(numCoeffs-1)))
  beta <- array(data=NA, dim=c(numCoeffs, numSpecies), 
                dimnames=list(namesCoeffs, namesSpecies))
  if ( numGears > 1 ) {
    zeta <- array(data=NA, dim=c(numGears-1, numSpecies), 
                  dimnames=list(namesZetas[-1], namesSpecies))
  } else {
    zeta <- NULL
  }
  coeffs <- list(beta=beta, zeta=zeta, gamma=NULL, delta=NULL)
  retVal$coeffs <- coeffs
  retVal$SE <- coeffs

  # Put the data together into one structure (except for species PA)
  dfData <- data.frame(cell = sampleCells,
                       gear = switch(as.character(is.null(sampleGears)),
                                     "TRUE" = rep(NA,numSamples),
                                     "FALSE" = sampleGears),
                       area = sampleAreas, covars)
  
  # Estimate the coefficients for each species' intensity function.
  for ( species in namesSpeciesPA) {
    # Get this species' PA data.
    dfData$Y <- PA[ ,species]

    # Estimate the coefficients for the model.
    glmFit <- glm(sdm.formula, family=binomial(link=cloglog), data=dfData, control=list(maxit=100))

    # Did the glm converge to an answer?
    if ( ! glmFit$converged ) {
      retVal$isErrors <- TRUE
      newError <- data.frame(species=species, msg="PA SDM has not converged.", stringsAsFactors = FALSE)
      retVal$errors <- rbind(retVal$errors, newError)
      return(retVal)
    } else {
      # Save estimated coefficient values and standard errors
      allSE <- summary(glmFit)$coefficients[ ,2]
      retVal$coeffs$beta[ ,species] <- glmFit$coefficients[1:numCoeffs]
      retVal$SE$beta[ ,species] <- allSE[1:numCoeffs]
      if ( numGears > 1 ) {
        retVal$coeffs$zeta[ ,species] <- glmFit$coefficients[(numCoeffs+1):(numCoeffs+numGears-1)]
        retVal$SE$zeta[ ,species] <- allSE[(numCoeffs+1):(numCoeffs+numGears-1)]
      }
      # retVal$estCoeffs[ ,species] <- glmFit$coefficients
      # retVal$stdErrors[ ,species] <- summary(glmFit)$coefficients[ ,2]
      # if ( ! is.null(sampleGears) )
      #   retVal$zeta[ ,species] <- glmFit$coefficients[(numCoeffs+1):(numCoeffs+numGears-1)]
      retVal$glmFit <- glmFit
    }
  }

  # Return results.
  return(retVal)

}

#-----------------------------------------------------------------------------------------

runSDM.PPM <- function(myFormula, PO, BG, domainObj, cellObj) {

  # GLM for presence-only data using PPM function.
  # NB: I think it might work out area itself! No, doesn't use area as they are points!!!!
  # NB: ONE SPECIES per run.
  # NB: need x and y as input for the PPM function even if they are not used in the formula.
  #
  # ARGUMENTS ...
  # myFormula: a formula object that contains how to combine the information for an
  #            intensity function.
  # PO:        the presence-only data points in a data frame (two columns: cell, species).
  # BG:        the quadrature or background data points in a data frame (columns: x, y, covar1, ...).
  # domainObj: a domain object that contains the mask and an owin object giving the region
  #            within the extent that is available.
  # cellObj:   a cell object that contains information at a cell level (see "cellsObj.r")

  # Initialise return value.
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retVal <- list(errors=errors, 
                 isErrors=FALSE,
                 glmFit=NULL)

  # Numbers of things.
  namesSpecies <- cellObj$namesSpecies
  POSpecies <- sort(unique(PO$species))
  if ( ! setequal(POSpecies, intersect(POSpecies, namesSpecies)) ) {
    retVal$isErrors <- TRUE
    newError <- data.frame(species=species, 
                           msg="Problem with species identification in PO data.", 
                           stringsAsFactors = FALSE)
    retVal$errors <- rbind(retVal$errors, newError)
    return(retVal)
  }
  numSpecies <- length(POSpecies)
  numCoeffs <- numCoefficients(myFormula)

  # Initialise the return value.
  estCoeffs <- matrix(nrow=numCoeffs, ncol=numSpecies)
  colnames(estCoeffs) <- POSpecies
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retVal$estCoeffs <- estCoeffs
  retVal$stdErrors <- estCoeffs

  # Get background (BG) data in right format for PPM function.  Same for all species.
  pppBG <- ppp(BG$x, BG$y, window=domainObj$owin)

  # Get the environmental covariate data at the presence-only cells
  # (background points now done elsewhere as same no matter what run!).
  namesCovars <- names(cellObj$covars)
  PO[ ,namesCovars] <- getCellVals(cellObj, PO$cell, "covars")

  # Make xy points for the PO data (uniformly random location within given cell)
  numPOPoints <- dim(PO)[1]
  resCell <- res(domainObj$mask)
  cellxrange <- c(0, resCell[1]) - (resCell[1] / 2.0)  # cell x boundaries shifted to have zero at centre!
  cellyrange <- c(0, resCell[2]) - (resCell[2] / 2.0)  # cell y boundaries shifted to have zero at centre!
  PO[ ,c("x","y")] <- getCellVals(cellObj, PO$cell, "xy")  # cell centres!
  x <- runif(numPOPoints, min=cellxrange[1], max=cellxrange[2])
  y <- runif(numPOPoints, min=cellyrange[1], max=cellyrange[2])
  PO$x <- x + PO$x
  PO$y <- y + PO$y

  # Remove response and replace x and y as PPM has these reserved.
  # No need for an area offset as these are point observations, not area observations!
  sdm.formula <- delete.response.formula(myFormula)
  sdm.formula <- replaceVarName(sdm.formula, "x", "x1")  # This works for formula without x!
  sdm.formula <- replaceVarName(sdm.formula, "y", "y1")

  # Estimate coefficients but only for the species in PO data.
  for ( species in POSpecies) {
    # Get this species' PO data.
    indWhichSpecies <- which(PO$species == species)
    thisPO <- PO[indWhichSpecies, c("x","y", namesCovars)]

    # Get presence-only (PO) data in right format for PPM function.
    pppPO <- ppp(thisPO$x, thisPO$y, window=domainObj$owin)

    # Create a quadrature scheme that contains the presence-only points and the
    # quadrature (or background) points.
    qsData <- quadscheme(data = pppPO, dummy = pppBG, method = "dirichlet")

    # Put both data together (in the same columns, not extra columns!).
    dfPPMData <- rbind(thisPO, BG[ , names(thisPO)])
    colnames(dfPPMData)[1] <- "x1"      # Rename column x as x1
    colnames(dfPPMData)[2] <- "y1"      # Rename column y as y1

    # Fit the PPM (non-stationary Poisson point process model).
    glmFit <- ppm(qsData, trend=sdm.formula, interaction=Poisson(), covariates=dfPPMData,
                  gcontrol=list(maxit=100))

    # Did the glm converge to an answer?
    if ( ! glmFit$internal$glmfit$converged ) {
      # Create error message ...
      retVal$isErrors <- TRUE
      newError <- data.frame(species=species, msg="PPM SDM has not converged.", stringsAsFactors = FALSE)
      retVal$errors <- rbind(retVal$errors, newError)
      return(retVal)
      
    } else {
      # Save estimates ...
      retVal$estCoeffs[ ,species] <- glmFit$coef
      retVal$stdErrors[ ,species] <- summary.ppm(glmFit)$coefs.SE.CI[ ,2]
      retVal$glmFit <- glmFit
    }
  }

  # Return results.
  return(retVal)

}

#-----------------------------------------------------------------------------------------

runSDM.MsPP <- function(lambda.formula, bias.formula, surveyCells, areaSurveys, PA, PO, BG,
                        cellObj) {

  # Run the MultispeciesPP species distribution model for PO + PA data.
  #
  # Arguments ...
  # lambda.formula: a formula object that contains how to combine the information for a
  #                 species intensity function.
  # bias.formula:   a formula object that contains how to combine the information for a
  #                 sample bias intensity function.
  # surveyCells:    the cell locations of the surveys (rows of PA data).
  # areaSurveys:    the areas of the surveys (in the same order as surveyCells and PA rows).
  # PA:             the presence-absence data in a data frame (numSurveys x numSpecies)
  # PO:             the presence-only data in a data frame (two columns, cell,species) 
  #                 with a row for each observation.
  # BG:             the background data in a data frame (columns: cell + covars)
  # cellObj:        a cell object that contains the environmental and sample bias covariates
  #                 for each cell in the domain

  # Initialise return value.
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retLst <- list(errors = errors, 
                 warnings = errors,
                 isErrors = FALSE)
  
  # Numbers and names of things.
  namesSpeciesAll <- cellObj$namesSpecies
  numSpeciesAll <- cellObj$numSpecies
  numEnvirs <- dim(cellObj$covars)[2]
  namesEnvirs <- names(cellObj$covars)
  numBiases <- dim(cellObj$biases)[2]
  namesBiases <- names(cellObj$biases)
  numCoeffs <- numCoefficients(lambda.formula)
  numDeltas <- numCoefficients(bias.formula, includeIntercept = FALSE)

  # Create a skeleton list to store the results from each SDM.
  beta <- array(data=NA, dim=c(numCoeffs, numSpeciesAll), 
                dimnames=list(NULL, namesSpeciesAll))
  gamma <- array(NA, dim=c(numSpeciesAll), dimnames=list(namesSpeciesAll))
  delta <- rep(NA, times=numDeltas)
  coeffs <- list(beta=beta, zeta=NULL, gamma=gamma, delta=delta)
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retLst$coeffs <- coeffs
  retLst$SE <- coeffs

  # Check survey cells and area surveys are the right length.
  numSurveys <- dim(PA)[1]
  if ( length(surveyCells) != numSurveys ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Argument surveyCells is not the correct length.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  }
  if ( length(areaSurveys) != numSurveys ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Argument areaSurveys is not the correct length.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  }
  
  # Set the number and name of the species to use.
  namesSpeciesPA <- names(PA)
  numSpeciesPA <- length(namesSpeciesPA)
  namesSpeciesPO <- sort(unique(PO$species))
  numSpeciesPO <- length(namesSpeciesPO)
  namesSpeciesInBoth <- intersect(namesSpeciesPA, namesSpeciesPO)
  numSpeciesInBoth <- length(namesSpeciesInBoth)
  namesSpeciesInEither <- unique(c(namesSpeciesPA,namesSpeciesPO))

  # Write an error message for those species that we lose because they aren't in both data sets.
  # This is an error that only MsPP needs (i.e. not Count, Cloglog or PPM).
  if ( numSpeciesInBoth <  length(namesSpeciesInEither) ) {
    namesSpeciesInOnlyOne <- namesSpeciesInEither[which(! namesSpeciesInEither %in% namesSpeciesInBoth)]
    for ( species in namesSpeciesInOnlyOne ) {
      # retLst$isErrors <- TRUE
      # newError <- data.frame(species=species,
      #                        msg="This species is available in only one of the data sets (PA or PO).",
      #                        stringsAsFactors = FALSE)
      # retLst$errors <- rbind(retLst$errors, newError)
      newWarning <- data.frame(species=species,
                               msg="This species is available in only one of the data sets (PA or PO).",
                               stringsAsFactors = FALSE)
      retLst$warnings <- rbind(retLst$warnings, newWarning)
    }
  }

  # Initialise temporary data frames for each type of data.
  dfPA <- data.frame(cell=surveyCells, stringsAsFactors = FALSE)
  dfPO <- PO
  dfBG <- BG
  #dfBG <- BG[ ,c("cell",namesEnvirs)]

  # Add the environmental covariate data at each of the data points.  BG version done elsewhere now.
  dfPA[ ,namesEnvirs] <- getCellVals(cellObj, surveyCells, item="covars")
  dfPO[ ,namesEnvirs] <- getCellVals(cellObj, PO$cell, item="covars")
  dfBG[ ,namesEnvirs] <- getCellVals(cellObj, BG$cell, item="covars")

  # Add the species presence-absence indicator columns (assume it is necessary for these
  # to be after the covariate values, otherwise, they could be added with cell above.)
  dfPA[ ,namesSpeciesInBoth] <- PA[ ,namesSpeciesInBoth]

  # Add the bias covariate values at each of the presence-only points and background points.
  # Exclude bias covariates that have already been added by environment covariates!
  indNewBiasNames <- ! is.element(namesBiases, namesEnvirs)
  newBiasNames <- namesBiases[indNewBiasNames]
  if ( length(newBiasNames) > 0 ) {
    tmp <- getCellVals(cellObj, PO$cell, item="biases")
    dfPO[ ,newBiasNames] <- tmp[ ,newBiasNames]
    tmp <- getCellVals(cellObj, BG$cell, item="biases")
    dfBG[ ,newBiasNames] <- tmp[ ,newBiasNames]
  }

  # Re-arrange the presence-only data into a list as required by multispeciesPP.
  lstPO <- list()
  for ( species in namesSpeciesInBoth ) {
    # Which rows of the PO data are this species?
    indWhichSpecies <- which(dfPO$species == species)

    # Add covariates (inc. x and y).
    lstPO[[species]] <- dfPO[indWhichSpecies, c("cell", namesEnvirs, newBiasNames)] # Everything but the species indicator column.
  }

  # Make the model specifications.
  sdmLambda.formula <- delete.response.formula(lambda.formula)
  sdmBias.formula <- delete.response.formula(bias.formula)

  # Call the function.
  areaDomain <- cellObj$areaCell * cellObj$numCells
  glmFit <- multispeciesPP( sdmLambda.formula, sdmBias.formula, dfPA, lstPO, dfBG,
                            quadrat.size=areaSurveys, region.size=areaDomain, control=list(maxit=100))

  # Did the glm converge to an answer?
  if ( ! glmFit$converged ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="MsPP SDM has not converged.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  } else {
    # Store the estimated coefficients for this fit.
    retLst$coeffs$beta[ ,namesSpeciesInBoth] <- glmFit$species.coef[1:numCoeffs, namesSpeciesInBoth] # alphas and betas
    retLst$coeffs$gamma[namesSpeciesInBoth] <- glmFit$species.coef[numCoeffs+1, namesSpeciesInBoth]  # gammas
    retLst$coeffs$delta <- glmFit$bias.coef                                                          # deltas
    
    # Store the standard errors for this fit.
    allSE <- readStdErrors.MsPP(glmFit$std.errs, numCoeffs, namesSpeciesInBoth,
                              length(glmFit$bias.coef), rownames(glmFit$species.coef)[1:numCoeffs])
    retLst$SE$beta[ ,namesSpeciesInBoth] <- allSE$beta[ ,namesSpeciesInBoth]
    retLst$SE$gamma[namesSpeciesInBoth] <- allSE$gamma[namesSpeciesInBoth]
    retLst$SE$delta <- allSE$delta
  }

  # Return results.
  return(retLst)

}

#-----------------------------------------------------------------------------------------

runSDM.Gear <- function(lambda.formula, bias.formula, gear.formula, rowsInCells, sampleAreas, 
                        PA, cellPO, sampleGears, cellsObj) {

  # Run the gearGLM species distribution model for PO + PA data.  PA data has gear 
  # type information for each sample but PO observations don't have gear type info.
  # 
  # Arguments ...
  # lambda.formula: a formula object that contains how to combine the information for an
  #                 intensity function.
  # bias.formula:   a formula object that contains how to combine the information for a
  #                 sample bias intensity function.
  # gear.formula:   a formula object that contains the gear specific component of the 
  #                 formula (not a valid formula on its own???).  Combined with intensity
  #                 formula if sampleGears data is given.
  # rowsInCells:    rows in cell object that contain the cells where the samples were taken.
  #                 (see surveysObj.r)
  # sampleCells:    a vector of cell numbers indicating where the PA data was collected 
  #                 within the domain (i.e. which cells contain the sample areas).  In reality
  #                 this will need to be the row in cellsObj where the sample cells are 
  #                 located (i.e. use 'rowsInCells' in the surveys object).
  # sampleAreas:    vector of the area of the samples defined by cell number (numSamples)
  # PA:             matrix of PA results from each sample (same order as cells) rows are
  #                 for each sample site and columns are for each species (numSamples x numSpecies)
  # cellPO:         the number of presence-only data per cell in a data frame (numCell*numSpecies).
  # sampleGears:    a vector containing an indicator for which gear type was used in each
  #                 sample (numSamples).  If this is NULL, then this data is not available.
  # cellsObj:       a cell object that contains the environmental and sample bias covariates
  #                 for each cell in the domain.
  #
  # NB: to trick this into doing exactly the same as MsPP, set gear.formula = NULL.

  # Initialise return value.
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retLst <- list(errors = errors, 
                 warnings = errors,
                 isErrors = FALSE)
  
  # Numbers and names of things.
  namesSpeciesAll <- cellsObj$namesSpecies
  numSpeciesAll <- cellsObj$numSpecies
  numEnvirs <- dim(cellsObj$covars)[2]
  namesEnvirs <- names(cellsObj$covars)
  numBiases <- dim(cellsObj$biases)[2]
  namesBiases <- names(cellsObj$biases)
  if ( is.null(gear.formula) ) {
    numGears <- 1
    namesGears <- cellsObj$namesGears[1]
  } else {
    numGears <- cellsObj$numGears
    namesGears <- cellsObj$namesGears
  }
  numCoeffs <- numCoefficients(lambda.formula)
  numDeltas <- numCoefficients(bias.formula, includeIntercept = FALSE)
  numSamples <- dim(PA)[1]
  
  # Create a skeleton list to store the results from each SDM.
  beta <- array(data=NA, dim=c(numCoeffs, numSpeciesAll), 
                dimnames=list(NULL, namesSpeciesAll))
  gamma <- array(NA, dim=c(numSpeciesAll), dimnames=list(namesSpeciesAll))
  if ( numDeltas > 0 ) {
    delta <- rep(NA, times=numDeltas)
  } else {
    delta <- NULL
  }
  if ( numGears > 1 ) {
    zeta <- array(data=NA, dim=c(numGears-1, numSpeciesAll), 
                  dimnames=list(namesGears[-1], namesSpeciesAll))
  } else {
    zeta <- NULL
  }
  coeffs <- list(beta=beta, zeta=zeta, gamma=gamma, delta=delta)
  errors <- as.data.frame(matrix(nrow=0, ncol=2), stringsAsFactors=FALSE)
  colnames(errors) <- c("species", "msg")
  retLst$coeffs <- coeffs
  retLst$SE <- coeffs

  # Check argument vectors relating to samples are the right length.
  if ( length(rowsInCells) != numSamples ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Argument rowsInCells is not the correct length.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  }
  if ( length(sampleAreas) != numSamples ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Argument sampleAreas is not the correct length.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  }
  if ( length(sampleGears) != numSamples ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Argument sampleGears is not the correct length.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
  }
  
  # Set the number and name of the species to use.
  namesSpeciesPA <- names(PA)
  numSpeciesPA <- length(namesSpeciesPA)
  namesSpeciesPO <- names(cellPO)
  numSpeciesPO <- length(namesSpeciesPO)
  namesSpeciesInBoth <- intersect(namesSpeciesPA, namesSpeciesPO)
  numSpeciesInBoth <- length(namesSpeciesInBoth)
  namesSpeciesInEither <- unique(c(namesSpeciesPA,namesSpeciesPO))
  
  # Write an error message for those species that we lose because they aren't in both data sets.
  # This is an error that only joint SDMs need (i.e. not Count, Cloglog or PPM).
  if ( numSpeciesInBoth <  length(namesSpeciesInEither) ) {
    namesSpeciesInOnlyOne <- namesSpeciesInEither[which(! namesSpeciesInEither %in% namesSpeciesInBoth)]
    for ( species in namesSpeciesInOnlyOne ) {
      newWarning <- data.frame(species=species,
                               msg="This species is available in only one of the data sets (PA/AB or PO).",
                               stringsAsFactors = FALSE)
      retLst$warnings <- rbind(retLst$warnings, newWarning)
    }
  }

  # Rejig formulae to work with gearglm function.
  lambda.response <- paste0("cbind(",paste0(namesSpeciesInBoth, collapse=","),")")
  lambda.terms <- paste0(attr(terms(lambda.formula),"term.labels"), collapse = "+")
  lambda.formula <- as.formula(paste0(lambda.response, "~", lambda.terms))
  bias.formula <- update(bias.formula, isP ~ .)
  if ( ! is.null(gear.formula) ) gear.formula <- update(gear.formula, ~ .)
  
  # Setup data
  d.po <- cbind.data.frame(cellPO[ ,namesSpeciesInBoth],  
                           cellsObj$covars[ ,namesEnvirs,drop=FALSE],
                           cellsObj$biases[ ,namesBiases,drop=FALSE],
                           gear = NA,
                           area = cellsObj$areaCell,
                           isP = TRUE)
  d.pa <- cbind.data.frame(PA[ ,namesSpeciesInBoth],
                           cellsObj$covars[rowsInCells,namesEnvirs,drop=FALSE],
                           matrix(0, nrow=numSamples, ncol=numBiases, dimnames = list(NULL, namesBiases)),
                           gear = sampleGears,
                           area = sampleAreas,
                           isP = FALSE)
  d <- rbind(d.pa,d.po)
  
  # Create A and B matrices and estimate coefficients using block glm.
  glmFit <- gearglm(lambda.formula, bias.formula, data = d, gear.formula, area = d$area,
                      method = "bglmQR", control = glm.control(maxit = 100)) 
    
  # Did the glm converge to an answer?
  if ( ! glmFit$converged ) {
    retLst$isErrors = TRUE
    newError <- data.frame(species="NA",
                           msg="Gear SDM has not converged.",
                           stringsAsFactors = FALSE)
    retLst$errors <- rbind(retLst$errors, newError)
    return(retLst)
    
  } else {
    # Unpack and store the estimated coefficients and SE for this fit.
    fitCoeffs <- glmFit$coefficients
    fitSE <- sqrt(diag(glmFit$cov))
    rowSpeciesEnds <- 0
    for ( sp in 1:numSpeciesInBoth ) {
      # This species name
      species <- namesSpeciesInBoth[sp]
      
      # This species coefficients and SEs.
      rowSpeciesStarts <- rowSpeciesEnds + 1
      rowSpeciesEnds <- sp * (numCoeffs + numGears)      # alpha + betas + zetas - 1 + gamma
      thisSpeciesCoeffs <- fitCoeffs[rowSpeciesStarts:rowSpeciesEnds]
      thisSpeciesSE <- fitSE[rowSpeciesStarts:rowSpeciesEnds]
      
      # Alpha and Betas ...
      rowsBeta <- 1:numCoeffs
      retLst$coeffs$beta[ ,species] <- thisSpeciesCoeffs[rowsBeta]
      retLst$SE$beta[ ,species] <- thisSpeciesSE[rowsBeta]  
      
      # Gammas ...  
      rowGamma <- numCoeffs + 1
      retLst$coeffs$gamma[species] <- thisSpeciesCoeffs[rowGamma]
      retLst$SE$gamma[species] <- thisSpeciesSE[rowGamma]  
      
      # Zetas ...
      if ( numGears > 1 ) {
        rowsZeta <- (rowGamma+1):(rowGamma+numGears-1)
        retLst$coeffs$zeta[ ,species] <- thisSpeciesCoeffs[rowsZeta]
        retLst$SE$zeta[ ,species] <- thisSpeciesSE[rowsZeta]
      }
    }
    
    # Deltas ...
    if ( numDeltas > 0 ) {
      rowsDelta <- (rowSpeciesEnds + 1):(rowSpeciesEnds + numDeltas)
      retLst$coeffs$delta <- fitCoeffs[rowsDelta]
      retLst$SE$delta <- fitSE[rowsDelta]
    }
  }
  retLst$fit <- glmFit

  # Return results.
  return(retLst)
  
}

#-----------------------------------------------------------------------------------------

readStdErrors.MsPP <- function(stdErrors, numCoeffs, namesSpecies, numDeltas,
                               coeffNames=1:numCoeffs){

  # Rearrange the stdErrors, given as a vector by MsPP, into the required format.

  # Species stuff.
  namesSpeciesSorted <- sort(namesSpecies)   # output from MsPP has species in alphanumeric order!
  numSpecies <- length(namesSpecies)

  # Initialise return value.
  retMat <- matrix(data=NA, nrow = numCoeffs, ncol = numSpecies)
  colnames(retMat) <- namesSpeciesSorted
  rownames(retMat) <- coeffNames
  retVal <- list(beta=retMat, gamma=matrix(NA,numSpecies,1), delta=NULL)
  rownames(retVal$gamma) <- namesSpeciesSorted

  # Unpack vector.
  for ( k in 1:numSpecies ) {
      ind <- (k-1)*(numCoeffs+1)
      retVal$beta[ ,k] <- stdErrors[(ind+1):(ind+numCoeffs)]
      retVal$gamma[k,1] <- stdErrors[ind+numCoeffs+1]
  }
  deltaStart <- ind + numCoeffs + 2
  deltaEnd <- deltaStart + numDeltas - 1
  retVal$delta <- stdErrors[deltaStart:deltaEnd]

  # Change order back to order in namesSpecies
  retVal$beta <- retVal$beta[ , namesSpecies]
  retVal$gamma <- retVal$gamma[namesSpecies,1]

  # Return Value
  return(retVal)
}

#-----------------------------------------------------------------------------------------

checkPrevalence <- function(namesSpecies, PA, PO, minPrevalence, cellPO=NULL){

  # Checks which species have enough data (using minPrevalence as the limit of acceptable).
  # Returns those species that have enough for each data set and for both data sets.
  # NB: Assumes that the AB (i.e. count) data is acceptable for the same species as the PA data!
  #
  # Arguments ...
  # namesSpecies:  the names of the species to be included in the simulation.
  # PA:            the presence-absence data points in a data frame (numSamples x numSpecies).
  # PO:            the presence-only observation points in a data frame (two columns: cell, species).
  # minPrevalence: minimum number of oberservations per species that are required to run glm.
  # cellPO:        the number of presence-only observations per species per cell in a 
  #                data frame (numCells x numSpecies).

  # Intialise return value.
  retVal <- list(namesSpeciesPA=namesSpecies, namesSpeciesPO=NA, whichPORows=1:dim(PO)[1],
                 isErrors=FALSE, errors=NULL, isWarnings=FALSE, warnings)

  # Numbers and names of things ...
  numSurveys <- dim(PA)[1]
  PASpecies <- colnames(PA)
  POSpecies <- unique(PO$species)
  notThesePORows <- NULL
  if ( ! is.null(cellPO) ) {
    cellPOSpecies <- colnames(cellPO)  # POSpecies should = cellPOSpecies but no check is made!
    retVal$namesSpeciesPO <- namesSpecies
  }

  # Check for species that are in data but not in names.
  if ( ! setequal(POSpecies, intersect(POSpecies, namesSpecies)) ) {
    retVal$isErrors <- TRUE
    newError <- data.frame(species=NA,
                           msg="Problem with species identification in PO data.",
                           stringsAsFactors = FALSE)
    retVal$errors <- rbind(retVal$errors, newError)
    return(retVal)
  }
  if ( ! setequal(PASpecies, intersect(PASpecies, namesSpecies)) ) {
    retVal$isErrors <- TRUE
    newError <- data.frame(species=NA,
                           msg="Problem with species identification in PA data.",
                           stringsAsFactors = FALSE)
    retVal$errors <- rbind(retVal$errors, newError)
    return(retVal)
  }

  for ( species in namesSpecies ) {
    # For PA data ...
    if ( is.na(match(species,PASpecies)) ) {
      # Missing species column in PA data?
      retVal$isWarnings <- TRUE
      newWarning <- data.frame(species=species,
                               msg="This species is missing from the AB/PA data columns.",
                               stringsAsFactors = FALSE)
      retVal$warnings <- rbind(retVal$warnings, newWarning)
      indWhich <- which(retVal$namesSpeciesPA == species)
      retVal$namesSpeciesPA <- retVal$namesSpeciesPA[-indWhich]

    } else if ( sum(PA[ ,species]) < minPrevalence ) {
      # Number of presences is less than the minimum prevalence required (including no presences!).
      retVal$isWarnings <- TRUE
      PAmsg <- paste("This species has less than", minPrevalence, "presences in the AB or PA data.")
      newWarning <- data.frame(species=species,
                               msg=PAmsg,
                               stringsAsFactors = FALSE)
      retVal$warnings <- rbind(retVal$warnings, newWarning)
      indWhich <- which(retVal$namesSpeciesPA == species)
      retVal$namesSpeciesPA <- retVal$namesSpeciesPA[-indWhich]

    } else if ( (numSurveys - sum(PA[ ,species])) < minPrevalence ) {
      # Number of absences is less than the minimum prevalence required (including no absences!).
      retVal$isWarnings <- TRUE
      PAmsg <- paste("This species has less than", minPrevalence, "absences in the AB or PA data.")
      newWarning <- data.frame(species=species,
                               msg=PAmsg,
                               stringsAsFactors = FALSE)
      retVal$warnings <- rbind(retVal$warnings, newWarning)
      indWhich <- which(retVal$namesSpeciesPA == species)
      retVal$namesSpeciesPA <- retVal$namesSpeciesPA[-indWhich]

    } else {
      # All good.
    }

    # For cell PO data ...
    if ( ! is.null(cellPO) ) {
      if ( is.na(match(species,cellPOSpecies)) ) {
        # Missing species column in PO data?
        retVal$isWarnings <- TRUE
        newWarning <- data.frame(species=species,
                                 msg="This species is missing from the cell PO data columns.",
                                 stringsAsFactors = FALSE)
        retVal$warnings <- rbind(retVal$warnings, newWarning)
        indWhich <- which(retVal$namesSpeciesPO == species)
        retVal$namesSpeciesPO <- retVal$namesSpeciesPO[-indWhich]
        
      } else if ( sum(cellPO[ ,species]) < minPrevalence ) {
        # Number of presences is less than the minimum prevalence required (including no presences!).
        retVal$isWarnings <- TRUE
        Pmsg <- paste("This species has less than", minPrevalence, "presences in the cell PO data.")
        newWarning <- data.frame(species=species,
                                 msg=Pmsg,
                                 stringsAsFactors = FALSE)
        retVal$warnings <- rbind(retVal$warnings, newWarning)
        indWhich <- which(retVal$namesSpeciesPO == species)
        retVal$namesSpeciesPO <- retVal$namesSpeciesPO[-indWhich]
        
      } else {
        # All good.
      }
    }
    
    # For PO data ...
    indWhichRows <- which(PO$species == species)
    if ( length(indWhichRows) == 0 ) {
      # No observations for this species.
      retVal$isWarnings <- TRUE
      newWarning <- data.frame(species=species,
                               msg="This species is missing from the PO observations.",
                             stringsAsFactors = FALSE)
      retVal$warnings <- rbind(retVal$warnings, newWarning)
      
    } else if ( length(indWhichRows) < minPrevalence ) {
      # Number of presences is less than the minimum prevalence required.
      retVal$isWarnings <- TRUE
      POmsg <- paste("This species has less than", minPrevalence, "presences in the PO observations.")
      newWarning <- data.frame(species=species,
                               msg=POmsg,
                               stringsAsFactors = FALSE)
      retVal$warnings <- rbind(retVal$warnings, newWarning)
      notThesePORows <- c(notThesePORows, indWhichRows)

    } else {
      # Enough data for this species.
    }

  }

  # Reset acceptable rows for PO data, if necessary.
  if ( ! is.null(notThesePORows) ) retVal$whichPORows <- retVal$whichPORows[-notThesePORows]

  # Return value.
  return(retVal)

}

#-----------------------------------------------------------------------------------------
