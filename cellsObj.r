#-----------------------------------------------------------------------------------------

initCellsObj <- function(){
  
  # Initialise or reset the cells object as a list with the required items but no values.
  
  # Create the object.
  obj <- list(numCells = 0,               # The number of cells in the domain.
              numRowsExt = 0,             # The number of rows in the extent of the domain (so can recover lost cell numbers).
              numColsExt = 0,             # The number of cols in the extent of the domain (ditto).   
              numSpecies = 0,             # The number of species populations.
              namesSpecies = NULL,        # A vector of the names of the species (used as column headings and/or identifiers).
              numGears = 0,               # The number of gear types used to detect/sample/collect/observe individuals in each population.
              namesGears = NULL,          # A vector of the names of the gear types (used as column headings and/or identifiers).
              resCell = NULL,             # A vector of length two that contains the width and length of a cell (width is x-axis and length is y-axis)
              areaCell = 0,               # The area of each cell in the domain (assumes area is the same for all cells).
              cells = NULL,               # A vector of length numCells that contains the cell number of the cells within the domain (cells of a raster).
              xy = NULL,                  # A two column data.frame with numCells rows where each row contains the centrepoint (x,y) of each cell.
              covars = NULL,              # A data.frame (numCells x numCovars) containing the environmental covariate values for each cell in the domain.
              biases = NULL,              # A data.frame (numCells x numBiases) containing the bias covariate values for each cell in the domain.
              trueLambda=NULL,            # The intensity (value of lambda) per cell per species (data.frame: numCells x numSpecies)
              trueD=NULL,                 # The thinning rate for each gear x species combination (data.frame: numGears x numSpecies)
              trueB=NULL,                 # The sample bias (value of b in Fithian, et al.) per cell per species (data.frame: numCells x numSpecies)
              N = NULL,                   # The number of individuals per cell per species for a single run (data.frame: numCells x numSpecies)
              P = NULL,                   # The number of individuals observed per cell per species for a single run (data.frame: numCells x numSpecies)
              PO = NULL,                  # The presence-only observations, one observation per row for a single run (data.frame: cell, species, gear)
              runNum = 0,                 # The run number that has created the values in N.
              isError = FALSE             # Error indicator
  )
  
  # Return value.
  return(obj)
}

#-----------------------------------------------------------------------------------------

is.cellsObj <- function(obj) {
  
  # Test whether the argument is a valid cells object.  Returns true if it is, false otherwise.
  # NB: only tests names of items at this stage, not classes of items!
  
  # Check argument is the right class (as far as we can!)
  if ( !is(obj, "list") ) {
    # The argument is not even a list.  It is not a cells object.
    return(FALSE)
  }
  
  # Get the expected names of the items for a cells object.
  objectItemNames <- names(initCellsObj())
  
  # Check the object has the same items.
  if ( all(names(obj) %in% objectItemNames) ) {
    # The same item names, hence, a valid cells object.
    return(TRUE)
    
  } else {
    # Not the same item names, hence, an invalid cells object.
    return(FALSE)
  }
  
}

#-----------------------------------------------------------------------------------------

makeCells <- function(mask, maskValue){
  
  # Make the valid cell numbers from the raster mask (assumes the raster structure of mask
  # is the same as the covariate rasters used elsewhere in the simulation).
  
  # Initialise the return value.
  obj <- initCellsObj()
  
  # Get the valid cell numbers from the mask raster.
  obj$cells <- cellFromMask(mask, mask, maskValue)
  
  # Get the centre point for each of these cells.
  obj$xy <- as.data.frame(xyFromCell(mask, cell=obj$cells), stringsAsFactors=FALSE)
  
  
  # Set the number of valid cells.
  obj$numCells <- length(obj$cells)[]
  
  # Set the number of columns and rows in the extent of the domain.
  obj$numRowsExt <- mask@nrows
  obj$numColsExt <- mask@ncols
  
  # Set the area of each cell.
  obj$resCell <- res(mask)
  obj$areaCell <- prod(obj$resCell)
  
  # Return value.
  return(obj)

}

#-----------------------------------------------------------------------------------------

makeCovarData <- function(cellsObj, covars, biases=NULL) {
  
  # Makes a data.frame version of the covars from the raster data.  
  # Returns the cell objected with the covariate data.
  #
  # Arguments ...
  # cellsObj: a cells object that contains the valid cells of the domain.
  # covars:   is a raster stack with numEnvirs layers that returns a valid value for all 
  #           points in the domain.
  # biases:   a raster stack with numBiases layers that returns a valid value for all 
  #           points in the domain.

  # The number of species.
  if ( is.null(cellsObj) || is.null(cellsObj$cells) || cellsObj$numCells == 0 ) {
    cellsObj$isError <- TRUE
    stop("Please make cells of the domain before extracting the covariate data of these cells.")
  }
  
  # The covariate values for the given cells.
  cellsObj$covars <- davesExtract.v3(covars, cellsObj$xy)
  
  # Convert to data.frame with column names.
  cellsObj$covars <- as.data.frame(cellsObj$covars, stringsAsFactors=FALSE)
  names(cellsObj$covars) <- names(covars)
  
  # Also for the sample bias covariates, if given.
  if ( ! is.null(biases) ) {
    cellsObj$biases <- davesExtract.v3(biases, cellsObj$xy)
    cellsObj$biases <- as.data.frame(cellsObj$biases, stringsAsFactors=FALSE)
    names(cellsObj$biases) <- names(biases)
  }
  
  # Return value.
  return(cellsObj)
  
}
  
#-----------------------------------------------------------------------------------------

makeIntensity <- function(cellsObj, myFormula, coeffs, namesSpecies=names(coeffs)) {
  
  # Calculates the value of the true intensity using the given arguments. 
  # Assumes there is a function 'lambda.cell' that calculates the intensity argument.
  # Assumes that the cellsObj already has values for those covariates in myFormula.
  # Results are returned in the given cells object.
  #
  # Assumes that there is an intensity function called "lambda.cell" (cell version of lambda).  
  #
  # Arguments ...
  # cellsObj:     a cells object that contains the valid cells of the domain.
  # myFormula:    a formula that specifies the form of the lambda function to be used.
  # coeffs:       is a matrix of coefficients for the intensity function lambda (with 
  #               nrows=numCoeffs and ncols=numSpecies)
  # namesSpecies: a vector of character strings containing the names of the species.
  
  
  # Check cells have been made.
  if ( is.null(cellsObj) || is.null(cellsObj$cells) || cellsObj$numCells == 0 ) {
    cellsObj$isError <- TRUE
    stop("Please make cells of the domain before calculating the true species intensities.")
  }
  
  # Check covars have been extracted.
  if ( is.null(cellsObj$covars) ) {
    cellsObj$isError <- TRUE
    stop("Please extract covariate data before calculating the true species intensities.")
  }

  # Convert formula to formula class if it is a character.
  if ( inherits(myFormula, "character") ) myFormula <- as.formula(myFormula)
  
  # Get the number of species.
  cellsObj$namesSpecies <- namesSpecies
  cellsObj$numSpecies <- length(namesSpecies)
  
  # Check this matches species columns in coeffs.
  if ( is.vector(coeffs) ) coeffs <- as.matrix(coeffs) 
  coeffColNames <- colnames(coeffs)
  if ( ! setequal(coeffColNames,namesSpecies) ) {
    cellsObj$isError <- TRUE
    stop("Species in columns of 'coeffs' does not match given species names.")
  }

  # Initialise the return values.
  cellsObj$trueLambda <- as.data.frame(matrix(nrow = cellsObj$numCells, ncol = cellsObj$numSpecies), 
                                  stringsAsFactors=FALSE)
  names(cellsObj$trueLambda) <- namesSpecies
  
  # Create each species' true intensity values in each cell.
  for ( species in namesSpecies) {
    cellsObj$trueLambda[ ,species] <- lambda.cell(myFormula, coeffs[ ,species], cellsObj$covars)
  }
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

setProbDet <- function(cellsObj, zeta, namesGears = dimnames(zeta)[[1]]) {
  
  # Sets the probability of detection per species x gear type in the given cells object.
  # Returns the given cells object with the given probability of detection.
  
  # Numbers of stuff.
  numGears <- length(namesGears)
  numSpecies <- dim(zeta)[2]
  if ( numSpecies != cellsObj$numSpecies ) 
    stop("Argument 'zeta' has the wrong number of species columns.")
  namesSpecies <- dimnames(zeta)[[2]]
  if ( ! setequal(namesSpecies, cellsObj$namesSpecies) && ! is.null(namesSpecies) ) 
    stop("Species in 'zeta' are not the same as those in cells object.")
  if ( numGears != dim(zeta)[1] ) 
    stop("Argument 'namesGears' has the wrong number of gears.")
  
  # Set probability of detection.
  cellsObj$trueD <- as.data.frame(exp(zeta))
  cellsObj$numGears <- numGears
  cellsObj$namesGears <- namesGears
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

makeProbObs <- function(cellsObj, myFormula, gamma, delta) {
  
  # Calculates the value of the true probability of observation (b) using the given arguments. 
  # Results are returned in the given cells object.
  #
  # Assumes that there is an intensity function called "lambda.cell" (cell version of lambda).
  # Formula is for the regression relationship ...
  #                      ln(b(s)) = gamma + delta * covars(s)
  #
  # Arguments ...
  # cellsObj:  a cells object that contains the valid cells of the domain.
  # myFormula: a formula that specifies the form of the lambda function to be used.
  # gamma:     a vector containing the coefficients related to species probability of 
  #            observation (length of numSpecies)
  # delta:     a vector containing the coefficients related to the human probability of
  #            observation (length of nlayer(biases))

  # Check cells have been made.
  if ( is.null(cellsObj) || is.null(cellsObj$cells) || cellsObj$numCells == 0 ) {
    cellsObj$isError <- TRUE
    stop("Please make cells of the domain before calculating the true species intensities.")
  }

  # Check sample bias covariates have been extracted.
  if ( is.null(cellsObj$biases) ) {
    cellsObj$isError <- TRUE
    stop("Please extract sample bias covariate data before calculating the true probability of observations.")
  }
  
  # Convert formula to formula class if it is a character.
  if ( inherits(myFormula, "character") ) myFormula <- as.formula(myFormula)
  
  # Check this matches species names in gamma.
  if ( ! setequal(names(gamma), cellsObj$namesSpecies) && ! is.null(names(gamma))) {
    cellsObj$isError <- TRUE
    stop("Species names in gamma do not match given species names.")
  }
  
  # Initialise the return values.
  cellsObj$trueB <- as.data.frame(matrix(nrow = cellsObj$numCells, ncol = cellsObj$numSpecies), 
                                  stringsAsFactors=FALSE)
  names(cellsObj$trueB) <- cellsObj$namesSpecies
  
  # Create each species' true probability of observation values in each cell.
  for ( species in cellsObj$namesSpecies) {
    thisSpeciesCoeffs <- c(gamma[species], delta)
    cellsObj$trueB[ ,species] <- lambda.cell(myFormula, thisSpeciesCoeffs, cellsObj$biases)
  }
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

makeNumIndivids <- function(cellsObj, runNum=0) {
  
  # Simulates the species populations using the given arguments.  Uses the Poisson
  # distribution to randomly generate the number of each species in each cell of the domain.
  # That is, N_ck ~ Pois(lambda * areaCell).
  #
  # Returns a cells object with the simulated species numbers per cell added in the item N.
  #
  # Arguments ...
  # cellsObj: a cells object that contains the valid cells of the domain.
  # runNum:   an integer that indicates which run the numbers of individuals are from
  #           (gets overwritten for each run)

  # Check cells have been made.
  if ( is.null(cellsObj) || is.null(cellsObj$cells) || cellsObj$numCells == 0 ) {
    cellsObj$isError <- TRUE
    stop("Please make cells of the domain before making species populations.")
  }
  
  # Check the true intensity for each species has been calculated.
  if ( is.null(cellsObj$trueLambda) ) {
    cellsObj$isError <- TRUE
    stop("Please make each species' true intensity before making species populations.")
  }
  
  # Initialise the return values.
  cellsObj$N <- as.data.frame(matrix(nrow = cellsObj$numCells, 
                                     ncol = cellsObj$numSpecies), 
                              stringsAsFactors=FALSE)
  names(cellsObj$N) <- cellsObj$namesSpecies
  
  # Create each species' population.
  for ( species in cellsObj$namesSpecies ) {
    # Calculate the expected number of individuals in each cell of the domain.
    # This assumes that the intensity (or value of lambda) is homogeneous within each cell.
    muN <- cellsObj$trueLambda[ ,species] * cellsObj$areaCell

    # Generate the number of individuals in each cell (NB: different for each run!)
    cellsObj$N[ ,species] <- rpois(cellsObj$numCells, muN)
  }
  
  # The run number.
  cellsObj$runNum <- runNum 
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

getCellVals <- function(cellsObj, cells, item="covars") {
  
  # Get the specified item's values at the given cell numbers.  'item' is the name of any
  # list item in cellsObj that has a first dimension with length = numCells (e.g. "xy", "trueLambda")
  # 
  # Returns the values associated with the given item for the given cell values.
  
  # Check values are available.
  if ( is.null(cellsObj[[item]]) ) {
    cellsObj$isError <- TRUE
    stop("These cell values have not yet been specified so can not be returned.")
  }

  # Initialise return value.
  numWantedCells <- length(cells)  
  if ( numWantedCells == 0 ) {
    cellsObj$isError <- TRUE
    stop("There don't appear to be any cells for which to get cell values.")
  }  
  # numCovars <- dim(cellsObj$covars)[2]
  # retVals <- as.data.frame(matrix(nrows=numWantedCells, ncols=numCovars), stringsAsFactors=FALSE)
  # names(retVals) <- names(cellsObj$covars)
  
  # Which rows are the wanted cell numbers in the pop list of valid cell nums?
  indRows <- match(cells, cellsObj$cells)
    #which(cellsObj$cells %in% cellNums)
  
  # Get values for these cells
  if ( inherits(cellsObj[[item]], "matrix") ) {
    # Must return as a matrix even if has one column.
    retVals <- as.matrix(cellsObj[[item]][indRows, ])  
    names(retVals) <- names(cellsObj[[item]])
    
  } else if ( inherits(cellsObj[[item]], "data.frame") ) {
    # Must return as a data.frame even if has one column.
    retVals <- as.data.frame(cellsObj[[item]][indRows, ]) 
    names(retVals) <- names(cellsObj[[item]])
    
  } else if ( inherits(cellsObj[[item]], "vector") ) { 
    retVals <- cellsObj[[item]][indRows]
    
  } else {
    cellsObj$isError <- TRUE
    stop("Unable to return cell values for this class of item.")
  }
  
  # Return value.
  return(retVals)

}

#-----------------------------------------------------------------------------------------

plotCellVals <- function(cellsObj, vals="trueLambda", cols=NULL, titles=NULL, ...) {
  
  # Plot the values in all cells of the cells object.  Uses xy values within cellsObj and
  # plots the requested values as a raster heat map (or image) style of plot.
  # 
  # Arguments ...
  # cellsObj: a cells object defined as a list of items
  # vals:     name of the item within the cells object to plot against xy values of cells
  #           OR data.frame of z values to plot (with numCell rows in same order as xy)
  # cols:     vector of name (string) or index (integer) of particular columns of data to
  #           plot (value of NULL will plot all).
  # titles:   vector of strings (same as number of columns to be used).  Uses column names 
  #           if NULL.
  
  # Get requested item's values from the list object.
  if ( is.character(vals) ) {
    allValues <- getCellVals(cellsObj, cellsObj$cells, vals)
  } else {
    allValues <- vals
  }
  
  # Restrict to just the required columns.
  if ( is.null(cols) ) cols <- 1:dim(allValues)[2]
  nPlots <- length(cols)

  # Check there are xy values.
  if ( is.null(cellsObj$xy) ) stop("There are no xy values in the given cell object.")
  
  # Titles for plots?
  if ( is.null(titles) ) titles <- names(allValues[ ,cols])
  if ( length(titles) != nPlots ) {
    if ( length(titles) == 1 ) {
      titles <- rep(titles, nPlots)
    } else {
      stop("The titles vector is the wrong length for the number of plots requested.")
    }
  }
  
  # Plot.
  for ( i in 1:nPlots ) {
    thisPlot <- cols[i]
    rlPlot <- rasterFromXYZ(cbind(cellsObj$xy, allValues[ ,thisPlot]), res = cellsObj$resCell)
    plot(rlPlot, ...)
    title(titles[i])
  }
  
}


#-----------------------------------------------------------------------------------------

makeSpeciesNumbers <- function(lambda, coeffs, covars, namesSpecies=colnames(coeffs),
                               nMeanIndivids=NULL, doPlotCheck=TRUE, doEstimateCheck=TRUE){
  
  # Simulates the numbers of each species present in each cell using the lambda formula,
  # coefficients and covariates.  Returns a cell object.
  #
  # OBSOLETE 14/08/2018.
  #
  # Arguments ...
  # lambda:        a formula that specifies the form of the intensity for each cell in the 
  #                domain per species.
  # coeffs:        a matrix that contains the coefficients of the formula for each species
  #                (nrow=nCoeffs, ncol=nSpecies).
  # covars:        a raster stack that contains the values of the covariates given in the 
  #                formula for each cell in the domain.
  # namesSpecies:  a vector that contains the names of the species (length=nSpecies)
  # nMeanIndivids: a scalar or vector that contains the mean number of indviduals for each
  #                species that are to be simulated (coeffs[1, ] will be scaled accordingly)
  # doPlotCheck:   produce plots to check simulated numbers.
  # doEstimateCheck: produce glm estimates of coefficients from simulated data and compare 
  #                to true coefficients.
  
  
  # Numbers of things.
  nCoeffs <- dim(coeffs)[1]
  nSpecies <- dim(coeffs)[2]
  
  # Check arguments.
  if ( ! inherits(lambda,"formula") ) lambda <- as.formula(lambda)
  if ( length(unique(namesSpecies)) != nSpecies ) 
    stop("Number of species in 'namesSpecies' does not match number in 'coeffs'.")
  if ( ! is.null(nMeanIndivids) ) {
    # Check it matches number of species or is a scalar.
    if ( length(nMeanIndivids) == 1 ) {
      nMeanIndivids <- rep(nMeanIndivids, nSpecies)
    } else if ( length(nMeanIndivids) != nSpecies ) {
      stop("Length of 'nMeanIndivids' should be number of species.")
    } else {
      # Use as is.
    }
  }
  
  # Take out alpha, wont use if nMeanIndivids is given.    
  beta <- coeffs
  if ( ! is.null(nMeanIndivids) ) beta[1, ] <- 0.0
  
  # Create cell object.
  cellsObj <- makeCells(domainObj$mask, domainObj$maskValue)
  cellsObj <- makeCovarData(cellsObj, dataObj$data)
  cellsObj <- makeIntensity(cellsObj, lambda, beta, namesSpecies)
  
  # Adjust alpha values, if necessary
  if ( ! is.null(nMeanIndivids) ) {
    # Use given alpha values.
    beta[1, ] <- log(nMeanIndivids / apply(cellsObj$trueLambda, 2, sum))
    expAlpha <- matrix(exp(beta[1, ]), nrow=cellsObj$numCells, ncol=nSpecies, byrow = TRUE)
    cellsObj$trueLambda <- cellsObj$trueLambda * expAlpha
  }
  # Check above works by using below as comparison!
  # if ( ! is.null(nMeanIndivids) ) {
  #   tmp <- NULL
  #   for (sp in 1:nSpecies ) {
  #     beta[1,sp] <- log(nMeanIndivids[sp]/sum(cellsObj$trueLambda[ ,sp]))
  #     tmp  <- cbind(tmp, cellsObj$trueLambda[ ,sp] * exp(beta[1,sp]))
  #   }
  #   if ( any(tmp - cellsObj$trueLambda > 1e-8)) stop("Not the same.")
  # }
  
  # Make the number of each species present in each cell.
  cellsObj <- makeNumIndivids(cellsObj, 1)
  
  # Estimation check ...
  if ( doEstimateCheck ) {
    # Compare estimated coefficients to true coefficients in a plot!
    estCoeffs <- matrix(nrow=nCoeffs, ncol=nSpecies, dimnames = dimnames(coeffs))
    offsetSDM <- rep(log(cellsObj$areaCell), cellsObj$numCells)
    for ( sp in 1:nSpecies ) {
      formSDM <- update.formula(lambda, N ~ .)
      dataSDM <- data.frame(cellsObj$covars, N=cellsObj$N[ ,sp], offset=offsetSDM)
      all.fit <- glm(formSDM, family=poisson(), data=dataSDM, offset=offset)
      estCoeffs[ ,sp] <- all.fit$coefficients
    }
    diffCoeffs <- abs((estCoeffs - beta)/beta)
    plot(c(1,nCoeffs), c(0, max(diffCoeffs)), ylim=c(0, max(diffCoeffs)), type="n", 
         xaxt="n", xlab="", ylab="abs(diff/true)")
    for ( sp in 1:nSpecies ) {
      points(1:nCoeffs, diffCoeffs[ ,sp], pch="_", col="red")
    }
    abline(h = 0, lty="dotted")
    axis(1, 1:nCoeffs, labels = dimnames(coeffs)[[1]], las=2)
    title("Comparison of true and estimated coefficients using nIndivids as data", xlab = "coeffs")
  }
  
  # Plot check
  if ( doPlotCheck ) {
    opar <- par(mfrow=c(2,2), mar=c(2,4,1,2), oma=c(0,0,1,0))
    doTitle <- TRUE
    for ( sp in 1:nSpecies ) {
      rlLambda <- rasterFromXYZ(cbind(cellsObj$xy, cellsObj$trueLambda[ ,sp]), res=cellsObj$resCell)
      plot(rlLambda, asp=1, ylab=paste0("lambda ", namesSpecies[sp]))
      rlN <- rasterFromXYZ(cbind(cellsObj$xy, cellsObj$N[ ,sp]), res=cellsObj$resCell)
      plot(rlN, asp=1, ylab=paste0("sim numbers ", namesSpecies[sp]))
      if ( doTitle ) {
        title("Simulated species numbers from lambda", outer=TRUE) 
        doTitle <- FALSE
      } else {
        doTitle <- TRUE
      }
    }
    par(opar)
  }
  
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------
  
makePOObservations <- function(cellsObj, runNum, 
                               numGears, namesGears, gearUseStrategy = NULL,
                               covar = cellsObj$covars[ ,1], meanGear = NULL, 
                               sdGear = NULL) {
  
  # Make the presence-only observations (numbers per cell P and list PO).  These are 
  # simulated from the simulated numbers of individuals for each species, the sample bias 
  # information and the average probability of detection per species.  That is, thin the 
  # species populations using N * trueB * sum_g(trueD).  Note that N are the simulated 
  # number of each species per cell (where N ~ Pois(trueLambda * areaCell)) as previously 
  # calculated by function 'makeNumIndivids'.  Thus, P ~ bin(N,trueB * sum_g(trueD)) for 
  # each cell.  FYI: p_{ckg} ~ Pois(\lambda_{ck} |A_c| b_{ck} \pi_{cg} d_{kg}) when 
  # gearUseStrategy != NULL. 
  #
  # Returns the cell object with the simulated number of presence-only data points per cell 
  # in 'P', and a data.frame 'PO' containing simulated data with an observed individual in 
  # each row.
  #
  # Arguments ...
  # cellsObj:         a cells object that contains the species populations.
  # runNum:           the current run number (as PO observations are created for each run).
  #                   Used to check that N stored in cellsObj are from the current run.
  # numGears:         the number of gears available (set in cellsObj here).
  # nameGears:        the names of each gear (set in cellsObj here).
  # gearUseStrategy:  the gear use strategy to apply for the PO data simulation.  If
  #                      NULL    - don't worry about gears for PO data, info is assumed not known.
  #                      covar   - NOT IMPLEMENTED! use cell covariate value to assign \pi_cg for each g. 
  # covar:            a vector of covariate data to use to weight gear usage (only necessary 
  #                   if gearUseStrategy = "covar").  
  # meanGear:         Sets the mean for each gear's normal distribution ( min(covar) < 
  #                   meanGear[gear] < max(covar) ).  When meanGear is NULL, the range of
  #                   the covariate is divided into numGears equal sections and meanGear
  #                   is the centre of each of these sections. Only necessary if 
  #                   gearUseStrategy = "covar".
  # sdGear:           Sets the standard deviation for each gear's normal distribution.  
  #                   When sdGear is NULL, all gears have an sd = half the width of the
  #                   sections (see meanGear). Only necessary if gearUseStrategy = "covar".
  
  # Initialise values
  numCells <- cellsObj$numCells
  
  # Check there are simulated populations.
  if ( is.null(cellsObj$N) || (runNum > cellsObj$runNum) ) {
    stop("Must simulate species cell numbers before making presence-only data.")
  }
  
  # Check there are probability of observations
  if ( is.null(cellsObj$trueB) || ! setequal(cellsObj$namesSpecies, names(cellsObj$trueB)) ) {
    stop("Please check that the probability of observing a species has been set correctly.")
  }
  
  # Check there are probability of detection
  if ( is.null(cellsObj$trueD) || ! setequal(cellsObj$namesSpecies, names(cellsObj$trueD)) ) {
    stop("Please check that the probability of detecting a species has been set correctly.")
  }
  
  # Check that numGears has been set and is a positive number.
  if ( numGears < 1 || is.null(numGears) || is.na(numGears)) 
    stop("The number of gear types is not a positive integer.") 
  
  # Check arguments ...
  if ( length(namesGears) != numGears ) 
    stop("Argument 'namesGears' is the wrong length.")
  
  # Initialise return values
  cellsObj$P <- as.data.frame(matrix(nrow=numCells, ncol=cellsObj$numSpecies, 
                                    dimnames=list(NULL, cellsObj$namesSpecies)))
  cellsObj$PO <- data.frame(cell=NULL, species=NULL, gear=NULL, stringsAsFactors = FALSE)
  
  # Work out PO data for each species' population.
  if ( is.null(gearUseStrategy) ) {
    avgTrueD <- apply(cellsObj$trueD, 2, mean)   # sum over gear types.
    for ( species in cellsObj$namesSpecies ) {
      # Calculate the number of individuals that are observed per cell.
      # NB: probability of successfully observing any of the N_ck individuals is b_ck * mean(d_kg)
      cellsObj$P[ ,species] <- rbinom(numCells, 
                                     cellsObj$N[ ,species], 
                                     cellsObj$trueB[ ,species] * avgTrueD[species])
      PO.cell <- rep(cellsObj$cells, times=cellsObj$P[ ,species])  
      numPO <- sum(cellsObj$P[ ,species])
      PO.species <- rep(species, times=numPO)
      PO <- data.frame(cell = PO.cell, 
                       species = PO.species, 
                       gear = rep(0, numPO), 
                       stringsAsFactors = FALSE)
      cellsObj$PO <- rbind(cellsObj$PO, PO)
 
      # Check none of the numbers observed in a cell are actually bigger than the number in the cell.
      if ( any(cellsObj$P[ ,species] > cellsObj$N[ ,species]) )
        warning(paste("Some of the cells have observed numbers greater than the number present!"))
      
      # temporary plot check
      # rlplot[cellsObj$cells] <- cellsObj$trueB[ ,species] * cellsObj$trueLambda[ ,species]
      # plot(rlplot, asp=1, main=paste0("Number of species ", species, " observed in each cell"))
      # indPOCells <- cellPOCount > 0
      # points(cellsObj$xy[indPOCells,1], cellsObj$xy[indPOCells,2], pch=as.character(cellPOCount[indPOCells]), col=cellPOCount[indPOCells])
      # title(sub="NB: sample biased intensity included")
    }
    
  # } else if ( gearUseStrategy == "covar" ) {
  #   # Gear types are assigned based on value of covariate in observation's cell.
  #   # Uses multiple overlapping normal distributions to calculate probability for each gear
  #   # being assigned in a cell.
  #   
  #   # Create normal distributions.
  #   if ( is.null(meanGear) || is.null(sdGear) ) {
  #     # Equally spread normal curves.
  #     rangeDomainCovar <- range(covar)
  #     binBreaks <- seq(from=floor(rangeDomainCovar[1]), to=ceiling(rangeDomainCovar[2]), 
  #                      length.out=numGears+1)
  #     halfBin <- abs(binBreaks[2] - binBreaks[1])/2.0
  #     means <- binBreaks[1:numGears] + halfBin     # i.e. centre of bins.
  #     stdevs <- rep(halfBin*2.0, numGears)
  #   } else {
  #     # User specified means and sds.
  #     means <- meanGear
  #     stdevs <- sdGear
  #   }
  #   
  #   # Get probability of using each gear in each cell.
  #   gearProbs <- matrix(nrow=numCells, ncol=numGears)
  #   for ( g in 1:numGears ) {
  #     gearProbs[ ,g] <- dnorm(covar, means[g], stdevs[g])
  #   }
  #   sumGearProbs <- matrix(apply(gearProbs, 1, sum), nrow=numCells, ncol=numGears)
  #   gearProbs <- gearProbs / sumGearProbs    # matrix where rows sum to one!
  #   
  #   # Select gear used in each observation.
  #   for ( ob in 1:numObs ) {
  #     cell <- cellObj$PO$cell[ob]
  #     cellsObj$PO$gear[ob] <- sample(1:numGears, size = 1, prob = gearProbs[cell, ]) 
  #   }
  #   
  #   # Work out what was observed based on assigned gear.
  #   for ( species in cellObj$namesSpecies ) {
  #     # Number of potential individuals that could be observed (exist in cell and observer in cell).
  #     potObs <- rbinom(cellsObj$numCells, cellObj$N[ ,species], cellsObj$trueB[ ,species])
  #     
  #     # Up to here 18/09/2018.
  #     # Number of potential individuals that could be observed by each gear in each cell.
  #     for ( cell in 1:numCells ) {
  #       potObsGear <- rmultinom(potObs[cell], 1, gearProbs[cell, ])
  #     }
  #     
  #     # Number of individuals observed (potential obs thinned by probability of detection).
  #   }
    
    
  } else {
    
    cellsObj$isError <- TRUE
    stop("Gear use strategy for PO data simulation not recognised.")
  }
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------

calcProbGearUse <- function(cellsObj, numGears, namesGears,
                            gearUseStrategy=c("rand", "prDet","covar"), 
                            covarName=names(cellsObj$covars)[1], meanGear=NULL, 
                            sdGear=NULL) {
  
  # Assign the gear type used to collect each observation.  There is a choice of the strategy 
  # used to assign gears.  Returns the cell object with the gear info included.
  # FYI: p_{ckg} ~ Pois(\lambda_{ck}|A_c|b_{ck}\pi_{cg}d_{kg}) 
  # 
  # OBSOLETE 18/09/2018
  #
  # Arguments ...
  # cellsObj:         the cells object that contains the number and location of the observations
  # numGears:         the number of gears available (set in cellsObj here).
  # nameGears:        the names of each gear (set in cellsObj here).
  # gearUseStrategy:  the strategy used to assign gears to observations, one of
  #                     rand  - uniform random assignment of one of the gears to each 
  #                             observation (\pi_{cg} = \pi, same for all observations)
  #                     prDet - use the probability of detection to weight the number of each
  #                             gear type assigned (\pi_{cg} = \pi_g, same for all cells).
  #                     covar - use the given covariate value to give probability of each 
  #                             gear being used for a particular sample.
  # covarName:        the name of the covariate to use to weight gear usage (only necessary 
  #                   if gearUseStrategy = "covar").  
  # meanGear:         Sets the mean for each gear's normal distribution ( min(covar) < 
  #                   meanGear[gear] < max(covar) ).  When meanGear is NULL, the range of
  #                   the covariate is divided into numGears equal sections and meanGear
  #                   is the centre of each of these sections. Only necessary if 
  #                   gearUseStrategy = "covar".
  # sdGear:           Sets the standard deviation for each gear's normal distribution.  
  #                   When sdGear is NULL, all gears have an sd = half the width of the
  #                   sections (see meanGear). Only necessary if gearUseStrategy = "covar".
  
  # Initialise values
  numObs <- dim(cellObj$PO)[1]
  numCells <- cellObj$numCells
  
  # Check that numGears has been set and is a positive number.
  if ( numGears < 1 || is.null(numGears) || is.na(numGears)) 
    stop("The number of gear types is not a positive integer.") 
  
  # Check arguments ...
  if ( length(namesGears) != numGears ) 
    stop("Argument 'namesGears' is the wrong length.")
  gearUseStrategy <- match.arg(gearUseStrategy)
  
  # Check other arguments have been given when they are required.
  if ( gearUseStrategy == "covar" ) {
    
    if ( is.null(covarName) ) 
      stop("Argument 'covarName' is missing but is required when gearUseStrategy = covar.")
    if ( ! covarName %in% names(cellsObj$covars) )
      stop("Argument 'covarName' is not one of the covariates present in cellsObj.")
    if ( ! is.null(meanGear) && length(meanGear) != numGears ) 
      stop("Argument 'meanGear' is the wrong length.")
    if ( ! is.null(meanGear) && length(sdGear) != numGears ) 
      stop("Argument 'sdGear' is the wrong length.")
  }
  
  # Assign gear type
  if ( gearUseStrategy == "rand") {
    # Random assignment of gear type to each observation (probably breaks model assumptions!).
    cellObj$PO$gear <- sample(1:numGears, numObs, replace=TRUE)
    
  } else if ( gearUseStrategy == "prDet" ) {
    # Gear types are assigned based on the probability of detection value (same for all cells).
    # Split into gear types.
    sumTrueD <- apply(cellsObj$trueD, 2, sum)   # sum over gear types.
    for ( species in cellsObj$namesSpecies ) {
      prGear <- cellsObj$trueD[ ,species] / sumTrueD[species]
      rowStart <- 1
      for ( cell in 1:numCells ) {
        numPOCell <- cellsObj$P[cell,species]
        if ( numPOCell > 0 ) {
          rowEnd <- rowStart + numPOCell - 1
          numPOCellGear <- rmultinom(1, numPOCell, prGear)
          PO$gear[rowStart:rowEnd] <- rep(1:numGears, times=numPOCellGear)
          rowStart <- rowEnd + 1
        }
      }
      cellsObj$PO <- rbind(cellsObj$PO, PO)
    }
    
  } else if ( gearUseStrategy == "covar" ) {
    # Gear types are assigned based on value of covariate in observation's cell.
    # Uses multiple overlapping normal distributions to calculate probability for each gear
    # being assigned in a cell.
    
    # Create normal distributions.
    if ( is.null(meanGear) || is.null(sdGear) ) {
      # Equally spread normal curves.
      rangeDomainCovar <- range(cellObj$covars[ ,covarName])
      binBreaks <- seq(from=floor(rangeDomainCovar[1]), to=ceiling(rangeDomainCovar[2]), 
                       length.out=numGears+1)
      halfBin <- abs(binBreaks[2] - binBreaks[1])/2.0
      means <- binBreaks[1:numGears] + halfBin     # i.e. centre of bins.
      stdevs <- rep(halfBin*2.0, numGears)
    } else {
      # User specified means and sds.
      means <- meanGear
      stdevs <- sdGear
    }
    
    # Get probability of using each gear in each cell.
    gearProbs <- matrix(nrow=cellObj$numCells, ncol=numGears)
    for ( g in 1:numGears ) {
      gearProbs[ ,g] <- dnorm(cellObj$covars[ ,covarName], means[g], stdevs[g])
    }
    sumGearProbs <- matrix(apply(gearProbs, 1, sum), nrow=numCells, ncol=numGears)
    gearProbs <- gearProbs / sumGearProbs
    
    # Select gear used in each observation.
    for ( ob in 1:numObs ) {
      cell <- cellObj$PO$cell[ob]
      cellsObj$PO$gear[ob] <- sample(1:numGears, size = 1, prob = gearProbs[cell, ]) 
    }
  }
  
  # Return value.
  cellsObj$numGears <- numGears
  cellsObj$namesGears <- namesGears
  return(cellsObj)
  
  
}

#-----------------------------------------------------------------------------------------

assignGearUsedObs <- function(cellsObj, numGears, namesGears,
                              gearUseStrategy=c("rand", "prDet","covar"), 
                              covarName=names(cellsObj$covars)[1], meanGear=NULL, 
                              sdGear=NULL) {
    
    # Assign the gear type used to collect each observation.  There is a choice of the strategy 
    # used to assign gears.  Returns the cell object with the gear info included.
    # FYI: p_{ckg} ~ Pois(\lambda_{ck}|A_c|b_{ck}\pi_{cg}d_{kg}) 
    # 
    # OBSOLETE 18/09/2018
    #
    # Arguments ...
    # cellsObj:         the cells object that contains the number and location of the observations
    # numGears:         the number of gears available (set in cellsObj here).
    # nameGears:        the names of each gear (set in cellsObj here).
    # gearUseStrategy:  the strategy used to assign gears to observations, one of
    #                     rand  - uniform random assignment of one of the gears to each 
    #                             observation (\pi_{cg} = \pi, same for all observations)
    #                     prDet - use the probability of detection to weight the number of each
    #                             gear type assigned (\pi_{cg} = \pi_g, same for all cells).
    #                     covar - use the given covariate value to give probability of each 
    #                             gear being used for a particular sample.
    # covarName:        the name of the covariate to use to weight gear usage (only necessary 
    #                   if gearUseStrategy = "covar").  
    # meanGear:         Sets the mean for each gear's normal distribution ( min(covar) < 
    #                   meanGear[gear] < max(covar) ).  When meanGear is NULL, the range of
    #                   the covariate is divided into numGears equal sections and meanGear
    #                   is the centre of each of these sections. Only necessary if 
    #                   gearUseStrategy = "covar".
    # sdGear:           Sets the standard deviation for each gear's normal distribution.  
    #                   When sdGear is NULL, all gears have an sd = half the width of the
    #                   sections (see meanGear). Only necessary if gearUseStrategy = "covar".
    
    # Initialise values
    numObs <- dim(cellObj$PO)[1]
    numCells <- cellObj$numCells

    # Check that numGears has been set and is a positive number.
    if ( numGears < 1 || is.null(numGears) || is.na(numGears)) 
      stop("The number of gear types is not a positive integer.") 
    
    # Check arguments ...
    if ( length(namesGears) != numGears ) 
      stop("Argument 'namesGears' is the wrong length.")
    gearUseStrategy <- match.arg(gearUseStrategy)
    
    # Check other arguments have been given when they are required.
    if ( gearUseStrategy == "covar" ) {
      
      if ( is.null(covarName) ) 
        stop("Argument 'covarName' is missing but is required when gearUseStrategy = covar.")
      if ( ! covarName %in% names(cellsObj$covars) )
        stop("Argument 'covarName' is not one of the covariates present in cellsObj.")
      if ( ! is.null(meanGear) && length(meanGear) != numGears ) 
        stop("Argument 'meanGear' is the wrong length.")
      if ( ! is.null(meanGear) && length(sdGear) != numGears ) 
        stop("Argument 'sdGear' is the wrong length.")
    }
    
    # Assign gear type
    if ( gearUseStrategy == "rand") {
      # Random assignment of gear type to each observation (probably breaks model assumptions!).
      cellObj$PO$gear <- sample(1:numGears, numObs, replace=TRUE)
  
    } else if ( gearUseStrategy == "prDet" ) {
      # Gear types are assigned based on the probability of detection value (same for all cells).
      # Split into gear types.
      sumTrueD <- apply(cellsObj$trueD, 2, sum)   # sum over gear types.
      for ( species in cellsObj$namesSpecies ) {
        prGear <- cellsObj$trueD[ ,species] / sumTrueD[species]
        rowStart <- 1
        for ( cell in 1:numCells ) {
          numPOCell <- cellsObj$P[cell,species]
          if ( numPOCell > 0 ) {
            rowEnd <- rowStart + numPOCell - 1
            numPOCellGear <- rmultinom(1, numPOCell, prGear)
            PO$gear[rowStart:rowEnd] <- rep(1:numGears, times=numPOCellGear)
            rowStart <- rowEnd + 1
          }
        }
        cellsObj$PO <- rbind(cellsObj$PO, PO)
      }
      
    } else if ( gearUseStrategy == "covar" ) {
      # Gear types are assigned based on value of covariate in observation's cell.
      # Uses multiple overlapping normal distributions to calculate probability for each gear
      # being assigned in a cell.
      
      # Create normal distributions.
      if ( is.null(meanGear) || is.null(sdGear) ) {
        # Equally spread normal curves.
        rangeDomainCovar <- range(cellObj$covars[ ,covarName])
        binBreaks <- seq(from=floor(rangeDomainCovar[1]), to=ceiling(rangeDomainCovar[2]), 
                         length.out=numGears+1)
        halfBin <- abs(binBreaks[2] - binBreaks[1])/2.0
        means <- binBreaks[1:numGears] + halfBin     # i.e. centre of bins.
        stdevs <- rep(halfBin*2.0, numGears)
      } else {
        # User specified means and sds.
        means <- meanGear
        stdevs <- sdGear
      }
      
      # Get probability of using each gear in each cell.
      gearProbs <- matrix(nrow=cellObj$numCells, ncol=numGears)
      for ( g in 1:numGears ) {
        gearProbs[ ,g] <- dnorm(cellObj$covars[ ,covarName], means[g], stdevs[g])
      }
      sumGearProbs <- matrix(apply(gearProbs, 1, sum), nrow=numCells, ncol=numGears)
      gearProbs <- gearProbs / sumGearProbs
      
      # Select gear used in each observation.
      for ( ob in 1:numObs ) {
        cell <- cellObj$PO$cell[ob]
        cellsObj$PO$gear[ob] <- sample(1:numGears, size = 1, prob = gearProbs[cell, ]) 
      }
    }
    
    # Return value.
    cellsObj$numGears <- numGears
    cellsObj$namesGears <- namesGears
    return(cellsObj)
    
    
}

#-----------------------------------------------------------------------------------------

resetSimDataCells <- function(cellsObj) {
  
  # Reset the simulated data items (N, PO and P) to null at the beginning of a run.
  
  cellsObj$N <- NULL
  cellsObj$P <- NULL
  cellsObj$PO <- NULL
  cellsObj$runNum <- 0
  
  # Return value.
  return(cellsObj)
  
}

#-----------------------------------------------------------------------------------------
