initialiseScenarios <- function(scenariosDir, numPA, gammaMultiplier=1.0, deltaMultiplier=1.0, 
                                zetaMultiplier = 1.0, 
                                numClusters=0, numRuns=10, sep="", prettyNames = NULL ){
  
  # Initialise or setup the scenarios object from the given arguments.  Arguments numPA, 
  # gammaMultiplier, deltaMultiplier, zetaMultiplier and numClusters and can each be scalars or vectors.
  # Argument scenariosDir is the directory where this group of scenarios is to be located
  # The unique combinations of these arguments form the settings for each scenario.
  # Returns the total number of scenarios, settings and name for each scenario.
  
  
  # Initialise return value.
  scenObj <- list(scenariosDir=scenariosDir,
                  numPA = numPA, 
                  gammaMultiplier = gammaMultiplier, 
                  deltaMultipler = deltaMultiplier, 
                  zetaMultiplier = zetaMultiplier,
                  numClusters = numClusters, 
                  numRuns = numRuns,
                  prettyNamesScenarios = prettyNames,  
                  #
                  # Things created within this function.
                  #
                  numSettings = 0,                # Number of scenario settings for each scenario (that can be vectors).
                  numSettingVals = NULL,          # a vector containing the number of values per setting to be tried.
                  numScenarios=0,                 # number of scenarios to be performed (calculated product of numSettingVals)
                  settings = NULL,                # a data.frame that contains the expanded settings with a row for each scenario and a column for each setting.
                  namesScenarios = NULL,          # a vector of strings containing a unique name for each scenario (can be used as directory name or on plots).
                  nameResultsFile = "Results",                  
                  isError=FALSE
                  )
  
  # What is the length of each setting argument?
  # This is where it will change if the number of settings changes (as well as function line).
  settingsLst <- list(nPA = numPA, GMult = gammaMultiplier, DMult = deltaMultiplier, 
                      nClust = numClusters, ZMult = zetaMultiplier)
  namesSettings <- names(settingsLst)           # Prefix used to create unique directorys.
  scenObj$numSettings <- length(settingsLst)
  scenObj$numSettingVals <- sapply(settingsLst, length)

  # How many scenarios need to be performed?
  scenObj$numScenarios <- prod(scenObj$numSettingVals)
  
  # Create settings for each unique scenario.
  scenObj$settings <- as.data.frame(matrix(nrow=scenObj$numScenarios, ncol=scenObj$numSettings),
                                    stringsAsFactors = TRUE)
  names(scenObj$settings) <- namesSettings
  scenObj$namesScenarios <- rep("", scenObj$numScenarios) 
  firstSettingVec <- TRUE
  for ( i in 1:scenObj$numSettings ) {
    # Repeat settings to form rows of unique combinations (in whole data.frame)
    if ( i < scenObj$numSettings ) {
      numEach <- prod(scenObj$numSettingVals[(i+1):scenObj$numSettings])
    } else {
      numEach <- 1
    }
    if ( i > 1 ) {
      numTimes <- prod(scenObj$numSettingVals[1:(i-1)])
    } else {
      numTimes <- 1
    }
    scenObj$settings[ ,i] <- rep(rep(settingsLst[[i]],each=numEach), times=numTimes)
    
    # Concatenate descriptive strings to form unique directory names for each scenario.
    # Don't bother for settings that have only one value.
    if ( scenObj$numSettingVals[i] > 1) {
      # Are the values for this setting doubles, characters or something else?  
      #strVals <- format(settingsLst[[i]], trim=TRUE)
      if ( is.double(settingsLst[[i]]) ) {
        strVals <- format(settingsLst[[i]], trim=TRUE, nsmall=1)
      } else if ( is.character(settingsLst[[i]]) ) {
        strVals <- settingsLst[[i]]
      } else {
        # Integers, logical.
        strVals <- format(settingsLst[[i]], trim=TRUE, nsmall=0)
      }
      if ( firstSettingVec ) {
        firstSettingVec <- FALSE
        strSep <- ""
      } else {
        strSep <- sep
      }
      strAdd <- rep(rep(paste0(strSep, namesSettings[i], strVals), each=numEach, times=numTimes))
      scenObj$namesScenarios <- paste0(scenObj$namesScenarios, strAdd)
    }
  }
  
  # Are there pretty names (print or plot worthy) provided by the user for each scenarios?
  if ( length(prettyNames) == 1 && scenObj$numScenarios > 1) {
    # Same scalar value for each, append number of scenario so names not all the same!
    scenObj$prettyNamesScenarios <- paste0(rep(prettyNames, scenObj$numScenarios),
                                           1:(scenObj$numScenarios)) 

  } else if ( is.null(prettyNames) ) {
    # No pretty names, use names created within this function.
    scenObj$prettyNamesScenarios <- scenObj$namesScenarios

  } else if ( length(prettyNames) == scenObj$numScenarios ) {
    # All pretty names supplied.
    # Check they are unique (as they will be used as identifiers for statistics and plotting).
    if ( length(unique(prettyNames)) < scenObj$numScenarios ) {
      scenObj$isError <- TRUE
      stop("There are duplicated values in the 'prettyNames' argument.  Duplicates are not permitted.")
    } else {
      scenObj$prettyNamesScenarios <- prettyNames  
    }

  } else {
    # Possible error in number of pretty names supplied?
    scenObj$isError <- TRUE
    stop("The length of the 'prettyNames' argument does not match the calculated number of scenarios.")
  }
  
  # Add the pretty names for plotting.
  scenObj <- addPrettyNamesToScenarios(scenObj, prettyNames)
  
  # Return value
  return(scenObj)

}

#-----------------------------------------------------------------------------------------

addPrettyNamesToScenarios <- function(scenObj, prettyNames=NULL) {

  # Add pretty names to the scenarios object.
  
  # Are there pretty names (print or plot worthy) provided by the user for each scenarios?
  if ( length(prettyNames) == 1 && scenObj$numScenarios > 1) {
    # Same scalar value for each, append number of scenario so names not all the same!
    scenObj$prettyNamesScenarios <- paste0(rep(prettyNames, scenObj$numScenarios),
                                           1:(scenObj$numScenarios)) 
    
  } else if ( is.null(prettyNames) ) {
    # No pretty names, use names created for directories.
    scenObj$prettyNamesScenarios <- scenObj$namesScenarios
    
  } else if ( length(prettyNames) == scenObj$numScenarios ) {
    # All pretty names supplied.
    # Check they are unique (as they will be used as identifiers for statistics and plotting).
    if ( length(unique(prettyNames)) < scenObj$numScenarios ) {
      scenObj$isError <- TRUE
      stop("There are duplicated values in the 'prettyNames' argument.  Duplicates are not permitted.")
    } else {
      scenObj$prettyNamesScenarios <- prettyNames  
    }
    
  } else {
    # Possible error in number of pretty names supplied?
    scenObj$isError <- TRUE
    stop("The length of the 'prettyNames' argument does not match the calculated number of scenarios.")
  }

  # Return value
  return(scenObj)

}

#-----------------------------------------------------------------------------------------

addPrettyNamesToStats <- function(statsObj, prettyNames=NULL) {
  
  # Add pretty names to the statsObj object by replacing existing names on scenario dimension.  
  # Assume they have already been made or added to scenarios object, so do less error checking.
  
  # Are there pretty names (print or plot worthy) provided by the user for each scenario?
  numScenarios <- dim(statsObj$stats[[1]]$avg)[[1]]
  if ( length(prettyNames) != numScenarios ) {
    # Possible error in number of pretty names supplied?
    stop("The length of the 'prettyNames' argument does not match the scenarios dimension in stats.")
  }
  
  # Replace names in each stat.
  for ( i in 1:statsObj$numStats ) {
    # Replace names in avg.
    dimnames(statsObj$stats[[i]]$avg)[[1]] <- prettyNames
    
    # Replace names in sd.
    dimnames(statsObj$stats[[i]]$sd)[[1]] <- prettyNames
  }
  
  # Return value
  return(statsObj)
  
}

#-----------------------------------------------------------------------------------------

runScenarios <- function(scenariosObj, simObj, cellsObj, useSDM, BG, domainObj, plotObj,
                         gearUseStrategyPA, gearUseStrategyPO,
                         numCoresUse = max(1, detectCores()-1), 
                         isOutputToFile=TRUE, randomSeed=NULL) {
  
  # Run all the scenarios using the information contained within the scenarios object.

  # Start time.
  timeTotal <- Sys.time()
  
  # Set random seed.
  set.seed(randomSeed)

  # Initialise result storage.
  resLstAll <- list()
  timeStart <- rep(NA, scenariosObj$numScenarios)
  class(timeStart) <- class(timeTotal)
  
  # Run each scenario.
  for ( i in 1:scenariosObj$numScenarios ) {
    # Start time for scenario.
    timeStart[i] <- Sys.time()
    
    # Create a scenario directory and sub-directories for ouput. 
    scenarioDir <- paste0(scenariosObj$scenariosDir, "/", scenariosObj$namesScenarios[i])
    simObj <- setSimOutput(simObj, simObj$inputDir, "allmsg.txt", scenarioDir, 
                           simObj$isOutputToFile)
    
    # Set scenario specific information (i.e. change coefficients to this scenario's version!)  
    simObj <- setSimScenarioInfo(simObj, 
                                 scenariosObj$numRuns, 
                                 numCoresUse, 
                                 numSamples = scenariosObj$settings$nPA[i], 
                                 gammaMultiplier = scenariosObj$settings$GMult[i], 
                                 deltaMultiplier = scenariosObj$settings$DMult[i], 
                                 zetaMultiplier = scenariosObj$settings$ZMult[i], 
                                 numSurveys = scenariosObj$settings$nClust[i],
                                 gearUseStrategyPA = gearUseStrategyPA, 
                                 gearUseStrategyPO = gearUseStrategyPO,
                                 useSDM = useSDM)

    # Run this scenario
    rLst <- runScenario(simObj, cellsObj, BG, domainObj, plotObj) 

    # Save this scenario's results 
    resLstAll[[i]] <- rLst$resLst
    if ( dim(rLst$resLst$errors)[1] > 0 ) {
      message(paste0("Errors occurred in scenario ",i))
    } else {
      message(paste0("Completed scenario ",i))
    }

    # End the diversion of output to file, if necessary, and change back any options.
    message("Do clean up")
    doSimCleanUp(simObj)
  }

  # Finish time for last scenario.
  timeEnd <- Sys.time()
  timeTotal <- format(Sys.time() - timeTotal)
  retLst <- list(plotObj=plotObj, simObj=rLst$simObj, resLstAll=resLstAll, scenariosObj=scenariosObj,
                 timeStart=timeStart, timeEnd=timeEnd, timeTotal=timeTotal)
  
  # Return result.
  return(retLst)
  
}

#-----------------------------------------------------------------------------------------

runScenario <- function(simObj, cellObj, BG, domainObj, plotObj) {
  
  # Run a single scenario from the information in the arguments.
  # Repeats a simulation simObj$numRuns times and collects the results from each run.
  
  # Scenario specific settings ...

  # How many cores are there?  How many runs?
  numCoresUse <- simObj$numCoresUse
  if ( is.null(numCoresUse) || numCoresUse < 1 || numCoresUse >= detectCores() ) {
    numCoresUse <- max(1, detectCores() - 1)
  } 
  
  # Make the scenario specific parts of the cells object (that don't change per sim run).
  cellObj <- makeProbObs(cellObj, simObj$bias.formula, simObj$gamma, simObj$delta)
  if ( plotObj$probObs )
    plotCellValues(plotObj, domainObj$mask, cellObj$cells, cellObj$truePrObs,
                   fileNameStart = "CheckCellProbObs",  titleTxts = "True probability of observation for ")
  if ( plotObj$biasIntensity)
    plotCellValues(plotObj, domainObj$mask, cellObj$cells, cellObj$trueLambda * cellObj$truePrObs,
                   fileNameStart = "CheckCellBiasedIntensity",  titleTxts = "True biased intensity for ")
  cellObj <- setProbDet(cellObj, simObj$zeta)
  
  # Dump data necessary for an SDM run (for analysis if there is an error).
  fileName <- paste0("DataDump-AllRuns")
  fileName <- makeFileName(fileName, simObj$dataDumpDir, "RData")
  save(simObj, cellObj, BG, domainObj, file = fileName)
  
  
  # Run the simulation ...

  # Get the time.
  timeStartSim <- Sys.time()
  
  # # Turn warnings into errors for the duration of the simulation runs.
  # oldOptionsWarn <- getOption("warn")
  # options(warn=2)
  
  # Run repeat simulations using parallelisation.  (NB: use par1Res for debugging only!)
  #par1Res <- runSim(1, simObj, domainObj, cellObj, plotObj, BG)
  parRes <- mclapply(1:(simObj$numRuns), runSim, simObj, domainObj, cellObj, 
                     plotObj, BG, mc.cores=simObj$numCoresUse)
  
  # # Turn warnings back to their old setting.
  #options(warn=oldOptionsWarn)
  
  
  # Collect the results ...

  message(paste("Begin organise results from runs", Sys.time()))
  
  # Initialise the result list object.
  resLst <- initResults(simObj$numRuns, simObj$numCoeffs, simObj$numBiases, 
                        cellObj$namesSpecies, cellObj$namesGears)
  resLst <- setResTryMe(simObj, resLst)
  
  # Decode the results.
  for ( i in 1:(simObj$numRuns) ) {
    # Collect results.
    if ( is.null(parRes[[i]]) ) {
      message(paste0("No results for run ",i))
      resLst <- setResError(resLst, i, msg = "No results returned from core for this run.")
      
    } else if ( inherits(parRes[[i]], "try-error")) {
      message(paste0("No results for run ",i))
      resLst <- setResError(resLst, i, msg = parRes[[i]]$errors[ ,4])
      
    } else {
      message(paste0("Results for run ", i))
      resLst <- setResOneToMany(resLst, parRes[[i]], i)
    }
  }
  message(paste("End organise results from runs", Sys.time()))
  
  # Get the time.
  timeEndSim <- Sys.time()
  
  # Save times in results.
  resLst <- setResTimes(resLst, timeStartSim, timeEndSim, 2)
  
  # Save the results.
  SavedResults <- makeFileName("Results", simObj$resultsDir, "RData")
  save(resLst, simObj, cellObj, file=SavedResults)
  SavedSettings <- makeFileName("SimSettings", simObj$resultsDir, "txt")
  write.list(simObj, file = SavedSettings)
  SavedResults <- makeFileName("Results", simObj$resultsDir, "txt")
  strTextResults <- c("numRuns", "numSpecies", "namesSpecies", "errors", "timeSetup",
                      "timeSim", "namesValidSDM")
  write.list(resLst[strTextResults], file = SavedResults)
  
  
  # Finish ...

  # Return value.
  return(list(resLst=resLst, simObj=simObj, cellObj=cellObj, plotObj=plotObj))
  
}

#-----------------------------------------------------------------------------------------

loadScenariosResults <- function(scenariosDir, nameResultsFile="Results") {
  
  # Use the given scenarios directory to load the results for these scenarios.
  # Assumes same directory structure as created by scenarios (see initialiseScenarios).

  # Directory where the data was saved ...  
  savedDir <- paste0(scenariosDir,"/Data")
  if ( ! dir.exists(savedDir) ) 
    stop(paste0("Unable to load data as scenarios directory doesn't exist.\n", 
                "  Directory tried: '", scenariosDir,"'"))
  
  # File where the data was saved ...
  savedFile <- makeFileName(nameResultsFile, savedDir,  "RData")
  
  # Load back into function's environment.
  currentEnvir <- environment()
  loadedObjectNames <- load(savedFile, envir = currentEnvir)
  
  # Reformat into a list for return from function.
  retLst <- mget(loadedObjectNames, envir = currentEnvir)
  # retLst <- list(plotObj=plotObj, simObj=simObj, resLstAll=resLstAll, scenariosObj=scenariosObj)
  
  # Return value.
  return(retLst)
  
}

#-----------------------------------------------------------------------------------------

saveScenariosResults <- function(scenariosDir, retLst, nameResultsFile="Results") {
  
  # Use the given scenarios directory to save the results for these scenarios.
  # Assumes same directory structure as created by scenarios (see initialiseScenarios).
  
  # Directory where the data will be saved ...  
  saveDir <- paste0(scenariosDir,"/Data")
  if ( ! dir.exists(saveDir) ) {
    dir.create(saveDir, recursive = TRUE)
  }
  
  # File where the data will be saved ...
  saveFile <- makeFileName(nameResultsFile, saveDir,  "RData")
  if ( file.exists(saveFile) ) {
    userAns <- readline("Saved results file exists.  Do you wish to overwrite?  Answer Y or N: ")
    if ( tolower(userAns) != "y" ) stop("Results not saved!")
  }
  
  # Assign return list objects as separate objects within the function's environment.
  saveObjects <- names(retLst)
  currentEnvir <- environment()
  for ( obj in saveObjects ) {
    assign(obj, retLst[[obj]], envir=currentEnvir)
  }
  
  # Save these objects to the file.
  save(list=saveObjects, file=saveFile, envir = currentEnvir)
  
}

#-----------------------------------------------------------------------------------------

numSuccessfulRunsScenarios <- function(scenariosObj, resLstAll, asPercentage=FALSE) {
  
  # Returns the number of successful runs (i.e. that produced coefficient estimates) per 
  # scenario x SDM for each species.
  #
  # Arguments ...
  # scenariosObj: A scenarios object (see initialiseScenarios)
  # resLstAll:    A vector of lists, one per scenario, that contains the results from the runs.
  # asPercentage: Whether or not to return the results as the number of successful runs or
  #               The percentage of successful runs (of the total number of runs)
  
  # Numbers and names of things.
  numScenarios <- scenariosObj$numScenarios
  namesScenarios <- scenariosObj$namesScenarios
  namesSpecies <- resLstAll[[1]]$namesSpecies
  numSpecies <- length(namesSpecies)
  numSDMs <- resLstAll[[1]]$numSDMs
  namesSDMs <- resLstAll[[1]]$namesValidSDM
  validSDMs <- resLstAll[[1]]$validSDM
  
  # Initialise return value.
  namesScenariosSDMs <- paste0(rep(namesSDMs, each=numScenarios), "-" ,
                               rep(namesScenarios, times=numSDMs))
  numSuccesses <- matrix(nrow=numScenarios*numSDMs, ncol=numSpecies, 
                         dimnames = list(namesScenariosSDMs, namesSpecies))
  
  # Calculate number of successful runs per scenario
  for ( i in 1:numScenarios ) {
    for ( s in 1:numSDMs ) {
      # Row in return value.
      thisRow <- i + (s-1)*numScenarios
      
      # Number of successful runs for this scenario x SDM combo for each species.
      resLstAll[[i]][[validSDMs[s]]] <- makeNumSuccessfulRuns(resLstAll[[i]][[validSDMs[s]]],
                                                              resLstAll[[i]][[validSDMs[s]]]$coeffs$beta)
      numSuccesses[thisRow,namesSpecies] <- resLstAll[[i]][[validSDMs[s]]]$numSuccessfulRuns[1,namesSpecies]
    }
  }
  
  # Return a percentage?
  if ( asPercentage ) {
    numRuns <- resLstAll[[1]]$numRuns
    numSuccesses <- numSuccesses / numRuns
  }
  
  # Return value.
  return(numSuccesses)
  
}

#-----------------------------------------------------------------------------------------

plotScenariosSummaries2 <- function(scenariosObj, trueCoeffs, resLstAll, 
                                    samePage=TRUE, ylimVals=NULL, plotSDMs=NULL, 
                                    absStats=TRUE, vioPlots=FALSE, diffPORanges=c(TRUE,TRUE),
                                    ...) {

  # Make the summary level plots.  Plots to assess unbiasedness (mean is near true solution) 
  # and efficiency (variance is not too great) per SDM, for all species x coefficients 
  # (scaled so that they make sense on the same plot).
  
  # Names and numbers of things.  FYI: assume all objects are required dimensions, no checking done!
  namesCoeffs <- rownames(trueCoeffs)
  numCoeffs <- length(namesCoeffs)
  numRuns <- scenariosObj$numRuns
  numScenarios <- scenariosObj$numScenarios
  namesScenarios <- scenariosObj$prettyNamesScenarios
  namesSpecies <- colnames(trueCoeffs)
  numSpecies <- length(namesSpecies)
  namesValidSDMs <- resLstAll[[1]]$namesValidSDM
  validSDMs <- resLstAll[[1]]$validSDM        # Assumes all scenarios use the same SDMs!
  numValidSDMs <- length(validSDMs)           # Assumes namesValidSDMs and validSDMs are same length!
  
  # Plot directory for scenario level plots.
  plotDir <- paste0(scenariosObj$scenariosDir, "/Plots")
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Initialise statistics to summarise the scenarios results (per SDM!).
  unbiasedness <- array(dim=c(numScenarios, numCoeffs, numSpecies, numValidSDMs), 
                        dimnames=list(namesScenarios, namesCoeffs, namesSpecies, namesValidSDMs))
  efficiencySD <- unbiasedness
  efficiencySE <- unbiasedness
  
  # Calculate the statistics ...
  for ( sdm in 1:numValidSDMs ) {
    nameSDM <- validSDMs[sdm]
    
    for ( species in namesSpecies ) {
      
      for ( j in 1:numCoeffs ) {
        nameCoeff <- namesCoeffs[j]
        
        # True coefficient.
        trueCoeff <- trueCoeffs[nameCoeff,species]
        
        # Will need the estimates of this coefficient for the given species.
        # coeffEsts <- matrix(nrow=numRuns, ncol=numScenarios, 
        #                     dimnames = list(1:numRuns, namesScenarios))
        
        # Get data from each scenario.
        for (sc in 1:numScenarios) {
          resLst <- resLstAll[[sc]]
          
          # All estimates of this species x this coeff for this scenario (i.e. estimates from each run)
          coeffEsts <- resLst[[nameSDM]]$coeffs[nameCoeff,species, ]
          
          # Which runs were successful for this species and this coefficient?
          indSuccessfulRuns <- which(! is.na(coeffEsts))

          # Statistics data.
          if ( length(indSuccessfulRuns) == 0 ) {
            msg <- paste0("No statistics available as all runs have failed: \n",
                          "    SDM      = ", namesValidSDMs[sdm], "\n",
                          "    species  = ", species, "\n",
                          "    coeff    = ", nameCoeff, "\n",
                          "    scenario = ", namesScenarios[sc])
            warning(msg)
          } else {
            # unbiasedness and efficiency using standard errors ...
            avgCoeffEst <- mean(coeffEsts[indSuccessfulRuns])
            unbiasedness[sc,nameCoeff,species,sdm] <- (avgCoeffEst - trueCoeff) / trueCoeff
            avgCoeffSE <-  mean(resLst[[nameSDM]]$SE[nameCoeff,species,indSuccessfulRuns])
            efficiencySE[sc,nameCoeff,species,sdm] <- avgCoeffSE / trueCoeff
            
            # Can form an efficiency from the standard deviation of estimates?
            if ( length(indSuccessfulRuns) > 1 ) {
              efficiencySD[sc,nameCoeff,species,sdm] <- sd(coeffEsts[indSuccessfulRuns]) / trueCoeff
            } else {
              msg <- paste0("Unable to form a standard deviation from only one successful run: \n",
                            "    SDM      = ", namesValidSDMs[sdm], "\n",
                            "    species  = ", species, "\n",
                            "    coeff    = ", nameCoeff, "\n",
                            "    scenario = ", namesScenarios[sc])
              warning(msg)
            }
          }
        }
      }
    }
  }
  
  ### Plotting ...
  
  # Which SDMs to plot.
  if ( is.null(plotSDMs) ) plotSDMs <- namesValidSDMs
  
  # Plot for unbiasedness
  internalPlotScenariosStat2(unbiasedness, "unbiasedness", samePage = samePage, plotSDMs = plotSDMs, 
                             ylimVals = ylimVals, diffPORanges = diffPORanges, absStats = absStats,
                             vioPlots = vioPlots, plotDir = plotDir, ...)
  
  # Plot for efficiency using SD of estimates.
  internalPlotScenariosStat2(efficiencySD, "efficiencySD", samePage = samePage, plotSDMs = plotSDMs, 
                             ylimVals = ylimVals, diffPORanges = c(FALSE,FALSE), absStats = absStats,
                             vioPlots = vioPlots, plotDir = plotDir, ...)
  
  # Plot for efficiency using SEs from results.
  internalPlotScenariosStat2(efficiencySE, "efficiencySE", samePage = samePage, plotSDMs = plotSDMs, 
                             ylimVals = ylimVals, diffPORanges = c(FALSE,FALSE), absStats = absStats,
                             vioPlots = vioPlots, plotDir = plotDir, ...)
  
# 
#   # Return statistics ...
#   return(list(accuracy=unbiasedness, efficiencySD=efficiencySD, efficiencySE=efficiencySE))
}

#-----------------------------------------------------------------------------------------

internalPlotScenariosStat <- function(stats, statName, plotObj, plotDir, samePage=FALSE, 
                                      ylimVals=NULL, redAlphas=FALSE, labels=FALSE,
                                      boxPlots=FALSE){
  
  # Internal plotting function to plot the summary statistics for scenarios.
  #
  # Arguments ...
  # stats:     A 4D array that contains the values to be plotted (numScenarios x numCoeffs 
  #            x numSpecies x numValidSDMs)
  # statName:  Name of statistic (used in file names and titles).
  # samePage:  Whether the plots are all on the same page (one file) or separate pages 
  #            (multiple files)  
  # ylimVals:  Indicates the range of the y-axis to be plotted.  Generally, when not NULL, 
  #            the y-axis is to be truncated to less than the actual values present.  When
  #            NULL, ylimVals is set to range(stats) so that all plots have same axis.
  # redAlphas: Colour only the points for the alpha coefficients red (otherwise use black)
  # labels:    Use labels (coefficient x species) instead of points for points outside
  #            y-axis range of [-1,1]
  # boxPlots:  
  
  # Numbers and names of things.
  dimStats <- dim(stats)
  numScenarios <- dimStats[1]
  numCoeffs <- dimStats[2]
  numSpecies <- dimStats[3]
  numSDMs <- dimStats[4]
  dimnamesStats <- dimnames(stats)
  namesScenarios <- dimnamesStats[[1]]
  namesCoeffs <- dimnamesStats[[2]]
  namesSDMs <- dimnamesStats[[4]]
  namesSpecies <- dimnamesStats[[3]]
  acceptableStatVal <- 1
  
  # Is the y-axis truncated (assumes ylimVals are less than actual range)
  if ( is.null(ylimVals) ) {
    fileTrunc <- ""
    titleTrunc <- ""
    ylimVals <- range(stats)   # Same y-axis for all plots.
  } else {
    fileTrunc <- "-trunc"
    titleTrunc <- " - truncated"
    ylimVals <- ylimVals
  }

  # There will be a plot per valid SDM.  Are they on the same page (in the same file)?
  if ( samePage ) {
    numPages <- 1
    numPlotsPerPage <- numSDMs
    fileNames <- makeFileName(paste0(statName, fileTrunc), plotDir, plotObj$device)
    titleTxt <- paste0(statName," of all estimates", titleTrunc)
  } else {
    numPages <- numSDMs
    numPlotsPerPage <- 1
    fileNames <- makeFileName(paste0(namesSDMs, statName, fileTrunc), 
                                     plotDir, plotObj$device)
    titleTxt <- paste0(statName, " of all estimates for ", namesSDMs, " SDM", titleTrunc)
  }
  
  # Labels for points, should they be needed.
  if ( labels ) {
    # Create combination of coefficient and species as labels.
    shortCoeffNames <- c("a", paste0("b", 1:(numCoeffs-1)))
    labelMat <- matrix(nrow=numCoeffs, ncol=numSpecies, 
                       dimnames=list(shortCoeffNames, namesSpecies))
    labelMat[] <- paste0(rep(shortCoeffNames, times=numSpecies), rep(namesSpecies, each=numCoeffs))
  }
  
  # Mean of the statistics for each scenario x SDM combination.
  meanStats <- apply(stats, c(1,4), mean)

  # Start plots.
  indSDM <- 0
  opar <- par(mfcol=c(1,numPlotsPerPage))
  for ( page in 1:numPages ) {
    # Start plot layers for this page ...
    layers <- setPlotLayers()
    xVals <- 1:(numScenarios * numPlotsPerPage)
    xlimVals <- c(0.5, max(xVals)+0.5)
    layers <- addPlotData(layers, xlimVals, ylimVals, plotfunc="plot", 
                          type="n", xaxs="i", xaxt="n", xlab="", ylab="", ylim=ylimVals)
    xLabs <- rep(namesScenarios, numPlotsPerPage)
    layers <- addPlotData(layers, 1, xVals, "axis", labels=xLabs, las=2, cex=0.5)
    layers <- addPlotData(layers, list(h=0), plotfunc="abline", lty=3)
    
    # Do layers that contain stats.
    for ( plot in 1:numPlotsPerPage ) {
      # Index for SDM to be used.
      indSDM <- indSDM + 1
      
      # X-axis values to be used for this plot.
      theseXVals <- (1:numScenarios) + ((plot - 1) * numScenarios)
      
      # Box plots or points and means?
      if ( boxPlots ) {
        # Add stats as box plots.
        for ( i in 1:numScenarios ) {
          # Add box plot layer for this scenario and this SDM.
          statsForBox <- as.vector(stats[i, , ,indSDM])
          layers <- addPlotData(layers, statsForBox, plotfunc="boxplot", at=theseXVals[i])
          
          # Highlight alpha coefficients by adding points and making them red?
          if ( redAlphas ) {
            layers <- addPlotData(layers, rep(theseXVals[i], times=numSpecies), 
                                  stats[i,1, ,indSDM], plotfunc="points", col="red")
          }
          
          # Do we need to plot labels?
          if ( labels ) {
            # Is this species x coefficients absolute value large enough 
            indLargeStat <- which(abs(stats[i, , ,indSDM]) > acceptableStatVal, arr.ind=TRUE)
            
            # Plot labels as well, if necessary.
            numLargeStats <- dim(indLargeStat)[1]
            if ( numLargeStats > 0) {
              # Get the approriate labels
              labelnames <- labelMat[indLargeStat]
              labelPos <- rep(3, numLargeStats)     # label above point (see pos in text help)
              statsMat <- stats[i, , ,indSDM]
              indNeg <- statsMat[indLargeStat] < (-acceptableStatVal)
              labelPos[indNeg] <- 1                 # label below point 
              layers <- addPlotData(layers, rep(theseXVals[i], times=numLargeStats), 
                                    statsMat[indLargeStat], plotfunc="text", labels=labelnames,
                                    pos=labelPos, cex=0.65)
            }
          } # if (labels )
        } # for ( i in 1:numScenarios )
        
      } else {
        # Add stats as points.
        for ( species in namesSpecies ) {
          for ( coeff in 1:numCoeffs ) {
            if ( redAlphas && coeff == 1 ) {
              # Separate out the alpha coefficients by making them red.
              layers <- addPlotData(layers, theseXVals, stats[ ,coeff,species,indSDM], 
                                    plotfunc="points", col="red")
            } else {
              layers <- addPlotData(layers, theseXVals, stats[ ,coeff,species,indSDM], 
                                    plotfunc="points", col="black")
            }
            
            # Do we need to plot labels?
            if ( labels ) {
              # What is this species x coefficient's label?
              indLargeStat <- which(abs(stats[ ,coeff,species,indSDM]) > acceptableStatVal)
              
              # Plot labels as well, if necessary.
              numLargeStats <- length(indLargeStat)
              if ( numLargeStats > 0) {
                label <- labelMat[coeff,species]
                labelPos <- rep(3, numLargeStats)       # label above point (see pos in text help)
                indNeg <- stats[indLargeStat,coeff,species,indSDM] < (-acceptableStatVal)
                labelPos[indNeg] <- 1                   # label below point 
                layers <- addPlotData(layers, theseXVals[indLargeStat], 
                                      stats[indLargeStat,coeff,species,indSDM], 
                                      plotfunc="text", labels=rep(label,times=numLargeStats),
                                      pos=labelPos, cex=0.65)
              }
            } # if (labels )
          } # for (coeff ...)
        } # for ( species ...)
        
        # Add mean line for each scenario.
        layers <- addPlotData(layers, theseXVals, meanStats[ ,indSDM], plotfunc="points", 
                              pch="_", col="green")
      }
      
      # Add division between SDM's if on same plot.
      if ( plot > 1 ) {
        xLine <- theseXVals[1] - 0.5
        layers <- addPlotData(layers, list(v=xLine), plotfunc="abline", lty="solid")
      }
    }
    
    # Do I need to add SDM titles to plots on same page?
    if ( numPlotsPerPage > 1 ) {
      x <- ((1:numSDMs) * numScenarios) + 0.5
      y <- rep(ylimVals[2], numSDMs)
      layers <- addPlotData(layers, x, y, plotfunc = "text", labels = namesSDMs, pos = 2)
    }
      
    # Plot layers ...
    plot.plotLayers(layers, plotObj$device, fileNames[page], plotObj$fileHeight, 
                    plotObj$fileWidth, titleTxt[page], xlab="")
    
  }

  # Return plotting to original settings.
  par(opar)
  
}

#-----------------------------------------------------------------------------------------

plotScenariosSummaries <- function(scenariosObj, plotObj, trueCoeffs, resLstAll, 
                                   truncPlots=FALSE, useBoxPlots=FALSE, redAlphas=FALSE, 
                                   SDMsOnSamePlot=FALSE) {

  # Make the summary level plots.  Plots to assess unbiasedness (mean is near true solution) 
  # and efficiency (variance is not too great) per SDM, for all species x coefficients 
  # (scaled so that they make sense on the same plot).
  
  # Names and numbers of things.  FYI: assume all objects are required dimensions, no checking done!
  namesCoeffs <- rownames(trueCoeffs)
  numCoeffs <- length(namesCoeffs)
  numRuns <- scenariosObj$numRuns
  numScenarios <- scenariosObj$numScenarios
  namesScenarios <- scenariosObj$namesScenarios
  namesSpecies <- colnames(trueCoeffs)
  numSpecies <- length(namesSpecies)
  namesValidSDMs <- resLstAll[[1]]$namesValidSDM
  validSDMs <- resLstAll[[1]]$validSDM        # Assumes all scenarios use the same SDMs!
  numValidSDMs <- length(validSDMs)           # Assumes namesValidSDMs and validSDMs are same length!
  
  # Plot directory for scenario level plots.
  plotDir <- paste0(scenariosObj$scenariosDir, "/Plots")
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Want a plot for each sdm x each statistic.
  for ( s in 1:numValidSDMs ) {
    nameSDM <- validSDMs[s]
    
    # Initialise statistics to summarise the scenarios results (per SDM!).
    unbiasedness <- array(dim=c(numScenarios, numCoeffs, numSpecies), 
                          dimnames=list(namesScenarios, namesCoeffs, namesSpecies))
    efficiencySD <- unbiasedness
    efficiencySE <- unbiasedness
    
    for ( species in namesSpecies ) {
      
      for ( j in 1:numCoeffs ) {
        nameCoeff <- namesCoeffs[j]
        
        # True coefficient.
        trueCoeff <- trueCoeffs[nameCoeff,species]
        
        # Will need the estimates of this coefficient for the given species.
        # coeffEsts <- matrix(nrow=numRuns, ncol=numScenarios, 
        #                     dimnames = list(1:numRuns, namesScenarios))

        # Get data from each scenario.
        for (i in 1:numScenarios) {
          resLst <- resLstAll[[i]]
          
          # All estimates of this species x this coeff for this scenario (i.e. estimates from each run)
          coeffEsts <- resLst[[nameSDM]]$coeffs[nameCoeff,species, ]

          # Which runs were successful for this species and this coefficient?
          indSuccessfulRuns <- which(! is.na(coeffEsts))
          
          # Statistics data.
          avgCoeffEst <- mean(coeffEsts[indSuccessfulRuns])
          unbiasedness[i,nameCoeff,species] <- (avgCoeffEst - trueCoeff) / trueCoeff
          efficiencySD[i,nameCoeff,species] <- sd(coeffEsts[indSuccessfulRuns]) / trueCoeff
          avgCoeffSE <-  mean(resLst[[nameSDM]]$SE[nameCoeff,species,indSuccessfulRuns])
          efficiencySE[i,nameCoeff,species] <- avgCoeffSE / trueCoeff
        }
      }
    }
  
    ### Plot for unbiasedness
    
    # Truncate plot?
    if ( truncPlots ) {
      thisFileName <- makeFileName(paste0(nameSDM,"unbiasedness-trunc"), plotDir, plotObj$device)
      titleTxt <- paste0("Unbiasedness of all estimates for ", nameSDM, " SDM - truncated")
      ylimVals <- c(-1.5,1.5)
    } else { 
      thisFileName <- makeFileName(paste0(nameSDM,"unbiasedness"), plotDir, plotObj$device)
      titleTxt <- paste0("Unbiasedness of all estimates for ", nameSDM, " SDM")
      ylimVals <- range(unbiasedness)
    }
    
    # Make layers for the plot.
    layers <- setPlotLayers()
    layers <- addPlotData(layers, c(1,numScenarios), range(unbiasedness), plotfunc="plot", 
                          type="n", xaxt="n", xlab="", ylab="", ylim=ylimVals)
    layers <- addPlotData(layers, 1, 1:numScenarios, "axis", labels=namesScenarios, las=2, cex=0.5)
    layers <- addPlotData(layers, list(h=0), plotfunc="abline", lty=3)
    layers <- makePlotLayersResStats(unbiasedness, layers, useBoxPlots, redAlphas)
    plot.plotLayers(layers, plotObj$device, thisFileName, plotObj$fileHeight, 
                    plotObj$fileWidth, titleTxt, xlab="")

    
    ### Plot for efficiency using SD of estimates.
    
    # Truncate plot?
    if ( truncPlots ) {
      thisFileName <- makeFileName(paste0(nameSDM,"efficiencySD-trunc"), plotDir, plotObj$device)
      titleTxt <- paste0("Efficiency of all estimates using SD for ", nameSDM, " SDM - truncated")
      ylimVals <- c(-10,10)
    } else {
      thisFileName <- makeFileName(paste0(nameSDM,"efficiencySD"), plotDir, plotObj$device)
      titleTxt <- paste0("Efficiency of all estimates using SD for ", nameSDM, " SDM")
      ylimVals <- range(efficiencySD)
    }
    
    # Make layers for the plot.
    layers <- setPlotLayers()
    layers <- addPlotData(layers, c(1,numScenarios), range(efficiencySD), plotfunc="plot", 
                          type="n", xaxt="n", xlab="", ylab="", ylim=ylimVals)
    layers <- addPlotData(layers, 1, 1:numScenarios, "axis", labels=namesScenarios, las=2, cex=0.5)
    layers <- addPlotData(layers, list(h=0), plotfunc="abline", lty=3)
    layers <- makePlotLayersResStats(efficiencySD, layers, useBoxPlots, redAlphas)
    plot.plotLayers(layers, plotObj$device, thisFileName, plotObj$fileHeight, 
                    plotObj$fileWidth, titleTxt, xlab="")
    
    ### Plot for efficiency using SEs from results.
    
    # Truncate plot?
    if ( truncPlots ) {
      thisFileName <- makeFileName(paste0(nameSDM,"efficiencySE-trunc"), plotDir, plotObj$device)
      titleTxt <- paste0("Efficiency of all estimates using SE for ", nameSDM, " SDM - truncated")
      ylimVals <- c(-10,10)
    } else {
      thisFileName <- makeFileName(paste0(nameSDM,"efficiencySE"), plotDir, plotObj$device)
      titleTxt <- paste0("Efficiency of all estimates using SE for ", nameSDM, " SDM")
      ylimVals <- range(efficiencySE)
    }
    
    # Make layers for the plot.
    layers <- setPlotLayers()
    layers <- addPlotData(layers, c(1,numScenarios), range(efficiencySE), plotfunc="plot", 
                          type="n", xaxt="n", xlab="", ylab="", ylim=ylimVals)
    layers <- addPlotData(layers, 1, 1:numScenarios, "axis", labels=namesScenarios, las=2, cex=0.5)
    layers <- addPlotData(layers, list(h=0), plotfunc="abline", lty=3)
    layers <- makePlotLayersResStats(efficiencySE, layers, useBoxPlots, redAlphas)
    plot.plotLayers(layers, plotObj$device, thisFileName, plotObj$fileHeight, 
                    plotObj$fileWidth, titleTxt, xlab="")
  }
  
}

#-----------------------------------------------------------------------------------------

plotScenarioCoeffEsts <- function(scenariosObj, plotObj, trueCoeffs, resLstAll, 
                                  plotSDMs = resLstAll[[1]]$namesValidSDM, 
                                  fileNameStart="EstCoeffsSpecies") {
  
  # Plot for each coefficient x species.  Each plot contains a boxplot for each sdm x scenario.
  
  # Names and numbers of things.  FYI: assume all objects are required dimensions, no checking done!
  namesCoeffs <- rownames(trueCoeffs)
  numCoeffs <- length(namesCoeffs)
  numRuns <- scenariosObj$numRuns
  numScenarios <- scenariosObj$numScenarios
  namesScenarios <- scenariosObj$namesScenarios
  namesSpecies <- colnames(trueCoeffs)
  numSpecies <- length(namesSpecies)
  namesValidSDMs <- resLstAll[[1]]$namesValidSDM
  validSDMs <- resLstAll[[1]]$validSDM
  numValidSDMs <- length(validSDMs)           # Assumes namesValidSDMs and validSDMs are same length!

  # Which SDMs to plot.
  # Figure out the number of SDMs to plot.
  numPlotSDMs <- length(plotSDMs)
  if ( ! all(plotSDMs %in% namesValidSDMs) ) {
    stop("Unrecognised SDM requested in plots.")
  }
  
  # Figure out the validSDMs needed from the plotSDMs requested. 
  # NB: validSDMs are the names used by the result object to store each SDM results.  This
  #     is different to the names used by the user to request various SDMs (clunky work 
  #     around as I changed the public names of the SDMs towards the end of the coding and
  #     it was too much trouble to change the internal names).
  indPlotSDMs <- match(plotSDMs, namesValidSDMs)
  plotValidSDMs <- validSDMs[indPlotSDMs] 
  
  # x-axis labels for the plots.
  namesSDMScenarios <- paste0(rep(plotSDMs, each=numScenarios),
                              rep(namesScenarios, times=numPlotSDMs))
  
  # Plot directory for scenario level plots.
  plotDir <- paste0(scenariosObj$scenariosDir, "/Plots")
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Want a plot for each species and each coefficient.
  for ( species in namesSpecies ) {
    
    for ( j in 1:numCoeffs ) {
      nameCoeff <- namesCoeffs[j]
      
      # True coefficient.
      trueCoeff <- trueCoeffs[nameCoeff,species]
      
      # Create storage for the data required for the boxplot into a single matrix.
      coeffEsts <- matrix(nrow=numRuns, ncol=numPlotSDMs*numScenarios, 
                          dimnames = list(1:numRuns, namesSDMScenarios))
      
      # Get data from each SDM x scenario combination.
      for (i in 1:numScenarios) {
        # This scenario's results.
        scenarioRes <- resLstAll[[i]]
        
        for ( s in 1:numPlotSDMs ) {
          nameSDM <- plotValidSDMs[s]
          
          # Index in storage matrix.
          indCol <- i + ((s-1)*numScenarios)
          
          # All estimates of this species x this coeff for this SDM x scenario (i.e. estimates from each run)
          coeffEsts[ ,indCol] <- scenarioRes[[nameSDM]]$coeffs[nameCoeff,species, ]
        }
      }

      # File name for the plot.
      thisFileName <- paste0(fileNameStart, "-", species, "-", nameCoeff)  
      thisFileName <- makeFileName(thisFileName, plotDir, plotObj$device)
      
      # Make layers for the plot.
      layers <- setPlotLayers()
      layers <- addPlotData(layers, coeffEsts, plotfunc="boxplot", las=2, cex=0.25, 
                            names=rep(namesScenarios,numPlotSDMs), 
                            xlim=c(1,(numPlotSDMs*numScenarios)) )
      layers <- addPlotData(layers, list(h=trueCoeff), plotfunc="abline", lty=3)
      
      # Do we need to add vertical line/s to separate SDMs into sub-plots?
      if ( numPlotSDMs > 1 ) {
        for ( i in 2:numPlotSDMs) {
          xVal <- ((i-1) * numScenarios ) + 0.5
          layers <- addPlotData(layers, list(v=xVal), plotfunc="abline", lty=1)
        }
      }
      
      # Add name of SDM to subplots.
      xVal <- (1:numPlotSDMs) * numScenarios
      yVal <- rep(max(coeffEsts, na.rm=TRUE), times=numPlotSDMs)
      layers <- addPlotData(layers, xVal, yVal, plotfunc="text", labels=plotSDMs)
      
      
      # Plot this coefficient for this species ...
      titleTxt <- paste0("True and estimated ", nameCoeff, " for ", species, " data.")
      
      # Plot the layers.
      plot.plotLayers(layers, plotObj$device, thisFileName, plotObj$fileHeight, 
                      plotObj$fileWidth, titleTxt, xlab="")
      
    }
  }

}

#-----------------------------------------------------------------------------------------

showAllScenariosErrors <- function(resLstAll) {
  
  # Prints to the console all the errors that are recorded in the results, for each scenario.
  
  # Number of scenarios
  numScenarios <- length(resLstAll)
  
  # Cycle through the results for each sceanrio.
  for ( i in 1:numScenarios ) {
    scenarioErrors <- resLstAll[[i]]$errors
    
    # Are there errors for this scenario?
    if ( dim(scenarioErrors)[1] > 0 ) {
      message("Errors that occurred in scenario ", i)
      print(scenarioErrors)
    } else {
      message("No errors recorded for scenario ", i)
    }
  }
  
}

#-----------------------------------------------------------------------------------------

showAllScenariosWarnings <- function(resLstAll) {
  
  # Prints to the console all the warnings that are recorded in the results, for each scenario.
  
  # Number of scenarios
  numScenarios <- length(resLstAll)
  
  # Cycle through the results for each sceanrio.
  for ( i in 1:numScenarios ) {
    scenarioWarnings <- resLstAll[[i]]$warnings
    
    # Are there warnings for this scenario?
    if ( dim(scenarioWarnings)[1] > 0 ) {
      message("Warnings that occurred in scenario ", i)
      print(scenarioWarnings)
    } else {
      message("No warnings recorded for scenario ", i)
    }
  }
  
}

#-----------------------------------------------------------------------------------------

makePlotLayersResStats <- function(resStats, layers, useBoxPlots=FALSE, redAlphas=FALSE){

  # Add the data layers for the given result statistic to the given layers that setup the plot.
  #
  # Arguments ...
  # resStats:    an array containing the result statistic (numScenarios x numCoeffs x numSpecies)
  # layers:      plot layers already setup to start plot (this function will add layers to this)
  # useBoxplots: by default, points of data are plotted but if this is true, boxplots will be used
  # redAlphas:   by default, no distinction is made between any of the coefficients.  
  #              However, if this is TRUE, then the alpha coefficient data points will be
  #              coloured red OR, if useBoxplots = TRUE, plotted as separate boxes.

  # Numbers and names of things.
  numCoeffs <- dim(resStats)[2]
  namesCoeffs <- dimnames(resStats)[[2]]
  numSpecies <- dim(resStats)[3]
  namesSpecies <- dimnames(resStats)[[3]]
  numScenarios <- dim(resStats)[1]
  namesScenarios <- dimnames(resStats)[[1]]

  if ( useBoxPlots ) {
    # Boxplots ... organise data into correct form.
    plotData <- as.data.frame(matrix(nrow = numSpecies*numCoeffs, ncol = numScenarios),
                              stringsAsFactors=FALSE)
    names(plotData) <- namesScenarios
    for ( j in 1:numCoeffs ) {
      for ( k in 1:numSpecies) {
        species <- namesSpecies[k]
        ind <- (j-1)*numSpecies + k
        plotData[ind, ] <- resStats[ ,j,species]
      }
    }
    
    # Do we need to separate out the alpha coefficients (i.e. make them red)?
    if ( redAlphas ) {
      # NB: alpha coefficients are the first numSpecies rows of data.
      layers <- addPlotData(layers, plotData[1:numSpecies, ], plotfunc="boxplot", 
                            names=rep("",numScenarios), border="red")
      layers <- addPlotData(layers, plotData[-(1:numSpecies), ], plotfunc="boxplot",
                            names=rep("",numScenarios))
    } else {
      layers <- addPlotData(layers, plotData, plotfunc="boxplot", names=rep("",numScenarios))
    }
    
  } else {
    # Points ... use data as is.
    for ( species in namesSpecies ) {
      for ( j in 1:numCoeffs ) {
        if ( redAlphas && j == 1 ) {
          # Separate out the alpha coefficients by making them red.
          layers <- addPlotData(layers, 1:numScenarios, resStats[ ,j,species], 
                                plotfunc="points", col="red")
        } else {
          layers <- addPlotData(layers, 1:numScenarios, resStats[ ,j,species], plotfunc="points")
        }
      }
    }
  }
  
  # Return plot layers.
  return(layers)
  
}

#-----------------------------------------------------------------------------------------

internalPlotScenariosStat2 <- function(stats, statName="statistics", samePage=TRUE,
                                       ylimVals=NULL, plotSDMs=dimnames(stats)[[4]], 
                                       diffPORanges=c(FALSE,FALSE), absStats=TRUE, vioPlots=FALSE,
                                       plotDevice="RStudioGD", plotDir=getwd(), 
                                       plotUnits="cm", plotHeight=9, plotWidth=12.25, 
                                       plotRes=600, doTitle=FALSE, outlierLabels=FALSE,
                                       xAxisTitle=NULL){
  
  # Internal plotting function to plot the summary statistics for scenarios.
  #
  # Arguments ...
  # stats:      A 4D array that contains the values to be plotted (numScenarios x numCoeffs 
  #             x numSpecies x numValidSDMs)
  # statName:   Name of statistic (used in file names and titles).
  # samePage:   Whether the plots are all on the same page (one file) or separate pages 
  #             (multiple files) NOT CODED YET!  TO DO !!!
  # ylimVals:   Indicates the range of the y-axis to be plotted.  Generally, when not NULL, 
  #             the y-axis is to be truncated to less than the actual values present.  When
  #             NULL, ylimVals is set to range(stats) so that all plots have same axis.
  # plotSDMs:   Which SDMs (and what order when on same page) to plot results from.
  # diffPORanges: if true, use different ranges for the y-axis in the PO SDM plot of 
  #             the alphas for diffPORanges[1] and/or betas for diffPORanges[2].
  # absStats:   if true, plot the absolute value of the statistics.
  # vioPlots:   if true, include violin plots behind the points.
  # plotDevice: What type of graphics device is to be used to plot.  "RStudioGD" is the 
  #             default and will display the plot in the RStudio plots tab.  See the 
  #             "grDevices" package for more help with devices but "png" is a good one!
  # plotDir:    Only when plotDevice not equal to "RStudioGD".  Directory where to save
  #             plot file.  File name is created as <statName>.<plotDevice>
  # plotUnits:  Only when plotDevice not equal to "RStudioGD".  Units in which height and
  #             width are given.
  # plotHeight: Only when plotDevice not equal to "RStudioGD".  Height of the plot in the 
  #             given units.
  # plotWidth:  Only when plotDevice not equal to "RStudioGD".  Width of the plot in the 
  #             given units.
  # plotRes:    Only when plotDevice not equal to "RStudioGD".  Resolution of the plot in
  #             pixels per inch (ppi).
  # doTitle:    Whether or not to add title to the plot.
  # outlierLabels: Whether or not to add outlier labels to plot.
  # xAxisTitle: Title to add to x-axis (at bottom of plot).  If NULL, nothing added.
  
  # Numbers and names of things.
  dimStats <- dim(stats)
  numScenarios <- dimStats[1]
  numCoeffs <- dimStats[2]
  numSpecies <- dimStats[3]
  numSDMs <- dimStats[4]
  dimnamesStats <- dimnames(stats)
  namesScenarios <- dimnamesStats[[1]]
  namesCoeffs <- dimnamesStats[[2]]
  namesSDMs <- dimnamesStats[[4]]
  namesSpecies <- dimnamesStats[[3]]

  # Colours for plots.
  pointsCol <- "black"
  violinCol <- "grey"
  pointsColPO <- colourBlindRGB("blue")
  violinColPO <- colourBlindRGB("skyBlue")
  
  # Figure out the number of SDMs to plot.
  numPlotSDMs <- length(plotSDMs)
  if ( ! all(plotSDMs %in% namesSDMs) ) {
    stop("Unrecognised SDM requested in plots.")
  }
  
  # Is the y-axis truncated (assumes ylimVals are less than actual range)
  if ( is.null(ylimVals) ) {
    fileTrunc <- ""
    titleTrunc <- ""
  } else {
    fileTrunc <- "-trunc"
    titleTrunc <- " - truncated"
  }
  
  # Which device are we printing to?
  if ( plotDevice == "RStudioGD" ) {
    # plot to the R studio plot window.
    plotToFile <- FALSE
  } else {
    fileName <- makeFileName(paste0(statName, fileTrunc), plotDir, plotDevice)
    argList <- list(filename = fileName, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
    plotToFile <- TRUE
  }
  
  # Divide plotting page into 2 rows (one for alphas and one for betas) and numSDMs columns.
  omaPar <- par("mar")         # use default single plot per page margins
  marPar <- c(0.0,0.0,1.3,0.0) # small margin on top to separate rows and allow for column titles.
  if ( is.null(xAxisTitle) ) {
    # No x-axis title, change margins to gain plotting area.
    omaPar[1] <- omaPar[1] - 1.0  # leave space for x-axis ticks and labels.  
  }
  if ( ! doTitle ) {
    # No main title, change margins to gain plotting area.
    omaPar[3] <- 0.0  # NB: there is already a bit at top due to mar settings.
  }
  opar <- par(mfrow=c(2,numPlotSDMs), oma = omaPar, mar = marPar, cex=0.66)
  
  # Convert stats to absolute values?
  if ( absStats ) stats <- abs(stats)
  
  # Plot the alpha coefficients in the first row.
  internalPlotRowScenarios(stats, plotSDMs, whichCoeffs=c(1), ylimVals = ylimVals, 
                           diffPORanges = diffPORanges[1], vioPlots = TRUE, 
                           rowTitle = "alphas", columnTitles = TRUE)
  
  # Plot the beta coefficients in the second row.
  internalPlotRowScenarios(stats, plotSDMs, whichCoeffs=2:numCoeffs, ylimVals = ylimVals, 
                           diffPORanges = diffPORanges[2], vioPlots = TRUE, rowTitle = "betas", 
                           xAxisLabels = namesScenarios)
  
  # Add a title to the plot.
  if ( doTitle ) 
    title(paste0("Coefficient ", statName, " for the given scenarios with the given SDMs", 
               titleTrunc), outer=TRUE)

  # Add an x-axis title to the plot.
  if ( ! is.null(xAxisTitle) ) title(xlab = xAxisTitle, outer=TRUE, line=min(4,omaPar[1]))
    
  # Turn off the plotting device, if it is a file.
  if ( plotToFile ) dev.off()
  
  # Return plotting to original settings.
  par(opar)
  
}

#-----------------------------------------------------------------------------------------

internalPlotRowScenarios <- function(stats, plotSDMs = dimnames(stats)[[4]], 
                                     whichCoeffs = 1:dim(stats)[2], ylimVals = NULL,
                                     diffPORanges = FALSE, vioPlots = FALSE, 
                                     rowTitle="stats", columnTitles=FALSE, 
                                     xAxisLabels=NULL, outlierLabels=FALSE) {
  
  # Plot a row (i.e. either alphas or betas) of the matrix of plots to compare the 
  # scenario results.

  # Numbers and names of things.
  dimStats <- dim(stats)
  numScenarios <- dimStats[1]
  numCoeffs <- dimStats[2]
  numSpecies <- dimStats[3]
  numSDMs <- dimStats[4]
  dimnamesStats <- dimnames(stats)
  namesScenarios <- dimnamesStats[[1]]
  namesCoeffs <- dimnamesStats[[2]]
  namesCoeffsShort <- c("a", paste0("b",1:(numCoeffs-1)))
  namesSDMs <- dimnamesStats[[4]]
  namesSpecies <- dimnamesStats[[3]]
  
  # Colours for plots.
  pointsCol <- "black"
  violinCol <- "grey"
  pointsColPO <- colourBlindRGB("blue")
  violinColPO <- colourBlindRGB("skyBlue")
  
  # Figure out the number of SDMs to plot.
  numPlotSDMs <- length(plotSDMs)
  if ( ! all(plotSDMs %in% namesSDMs) ) {
    stop("Unrecognised SDM requested in plots.")
  }
  
  # Y-axis range ... 
  yaxisRange <- matrix(nrow=2, ncol=numPlotSDMs, dimnames = list(1:2,plotSDMs))
  indPOSDM <- length(0)
  indOtherSDMs <- 1:numPlotSDMs
  if ( is.null(ylimVals) ) {
    if ( diffPORanges ) {
      # Different range for PO SDM.
      if ( "PO" %in% plotSDMs ) 
        yaxisRange[ ,"PO"] <- range(stats[ ,whichCoeffs, ,"PO"], na.rm = TRUE)
      
      # Same range for other SDMs.
      indPOSDM <- which(plotSDMs %in% "PO" )
      if ( length(indPOSDM) > 0 ) indOtherSDMs <- indOtherSDMs[-indPOSDM]
      if ( length(indOtherSDMs) > 0 ) {
        tmp <- range(stats[ ,whichCoeffs, ,plotSDMs[indOtherSDMs]], 0, na.rm = TRUE)  # include zero when stats are positive.
        yaxisRange[ ,plotSDMs[indOtherSDMs]] <- rep(tmp, length(indOtherSDMs))
      }
    } else {
      # Same range for all.
      tmp <- range(stats[ ,whichCoeffs, ,plotSDMs], 0, na.rm = TRUE)  # include zero when stats are positive.
      yaxisRange[ ,plotSDMs] <- rep(tmp, numPlotSDMs)
    }
  } else {
    # Same for all, whatever is specified.
    yaxisRange[1, ] <- rep(ylimVals[1], numPlotSDMs)
    yaxisRange[2, ] <- rep(ylimVals[2], numPlotSDMs)
  }
  
  # Calculate widths for violin plots in row.
  if ( vioPlots ) {
    # Get the maximum density for each SDM x scenario combination.
    maxDensity <- matrix(nrow=numScenarios, ncol=numPlotSDMs, dimnames=list(1:numScenarios, plotSDMs))
    for ( sdm in plotSDMs ) {
      
      for ( sc in 1:numScenarios ) {
        # What is the maximum of the density for each scenario.
        thisStats <- stats[sc,whichCoeffs, ,sdm]
        maxDensity[sc,sdm] <- max(density(thisStats)$y)
      }
    }
    
    # Scale maximum densities so that widths are between [0,1]
    widthsViolins <- maxDensity
    if ( ! is.null(ylimVals) || ! diffPORanges ) {
      # All sdm plots are on the same scale.
      maxMaxDensity <- max(maxDensity[])
      widthsViolins <- maxDensity/maxMaxDensity
      
    } else { # if ( is.null(ylimVals) && diffPOAlphaRanges )
      # PO sdm is on a different scale (NB: and ylimVals are not set by user).
      if ( length(indPOSDM) > 0 ) {
        maxMaxDensity <- max(maxDensity[ ,"PO"])
        widthsViolins[ ,"PO"] <- maxDensity[ ,"PO"] / maxMaxDensity
      }
      
      # Other SDMs are on the same scale (i.e. y-axis range).
      if ( length(indOtherSDMs) > 0 ) {
        maxMaxDensity <- max(maxDensity[ ,plotSDMs[indOtherSDMs]])
        widthsViolins[ ,plotSDMs[indOtherSDMs]] <- maxDensity[ ,plotSDMs[indOtherSDMs]] /
                                                     maxMaxDensity
      }
    }
  }
  
  # Work out what the labels would be for each species x coefficient combination.
  if ( outlierLabels ) {
    allLabels <- matrix(nrow=numCoeffs, ncol=numSpecies, dimnames=list(1:numCoeffs, namesSpecies))
    for ( i in 1:numCoeffs ) {
      allLabels[i, ] <- paste0(namesCoeffsShort[i], namesSpecies)
    }
  }
  
  # Plot row ...
  firstPlotInRow <- TRUE
  for ( sdm in plotSDMs ) {
    # Plot statistics for each SDM (i.e. column) for this row.
    plot(c(0.5,numScenarios+0.5),yaxisRange[ ,sdm], type="n", xaxs="i", xaxt="n", xlab="", 
         yaxt="n", ylab="")
    for ( sc in 1:numScenarios ) {
      # All species' coefficient stats at this scenario and SDM!
      thisStats <- as.vector(stats[sc,whichCoeffs, ,sdm])
      if ( outlierLabels ) thisLabels <- as.vector(allLabels[whichCoeffs, ])
      indNA <- which(is.na(thisStats))
      if ( length(indNA) > 0 ) {
        warning(paste0("One or more NA values in ", statName, " for SDM = ", sdm, 
                       " and scenario = ", namesScenarios[sc]))
        thisStats <- thisStats[-indNA]
        if ( outlierLabels ) thisLabels <- thisLabels[-indNA]
      }
      
      # If PO is to be a different range, do different y-axis for PO.
      if ( is.null(ylimVals) && diffPORanges && sdm == "PO" ) {
        if ( vioPlots ) vioplot(thisStats, col=violinColPO, border=violinColPO, add=TRUE,
                                drawRect = FALSE, at=sc, wex=widthsViolins[sc,sdm])
        axis(4, labels=TRUE, col=pointsColPO)
        stripchart(thisStats, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsColPO)
      } else {
        if ( vioPlots ) vioplot(thisStats, col=violinCol, border=violinCol, add=TRUE,
                                drawRect = FALSE, at=sc, wex=widthsViolins[sc,sdm])
        stripchart(thisStats, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsCol)
      }
      
      # Labels on outliers?
      if ( outlierLabels ) {
        # Use the interquartile range method to identify outliers.
        quartiles <- quantile(thisStats, probs=c(0.25,0.5,0.75), type = 8)
        intQuartRange <- quartiles[3] - quartiles[1]
        notOutlierRange <- c(quartiles[1] - 1.5*intQuartRange,
                             quartiles[3] + 1.5*intQuartRange)
        indOutliers <- which(thisStats < notOutlierRange[1] | thisStats > notOutlierRange[2])
        numOutliers <- length(indOutliers)
        if ( numOutliers > 0 ) {
          # FYI: cex=0.3 will only be readable if the plot resolution (i.e. pixels per inch) 
          #      is high enough (at least 600ppi).
          text(sc, thisStats[indOutliers], labels=thisLabels[indOutliers], pos=4, cex=0.3,
               offset=0.2)
        }
        
        # Median and boundaries between outliers and "normal" values.
        points(sc, quartiles[2], pch="_", col="green")
        points(sc, notOutlierRange[1], pch="_", col="blue")
        points(sc, notOutlierRange[2], pch="_", col="blue")
      }
    }
    abline(h=0, lty="dotted")
    
    
    # If this is the first column, do y-axis.
    if ( firstPlotInRow ) {
      axis(2, labels=TRUE)
      midyRange <- ((yaxisRange[2,sdm] - yaxisRange[1,sdm])/2.0) + yaxisRange[1,sdm]
      axis(2, at=midyRange, labels=rowTitle, tick=FALSE, padj=-2.0)
      firstPlotInRow <- FALSE
    }
    
    # If this is the PO column (this will cause all sorts of problems if PO is the first column!)
    if ( is.null(ylimVals) && diffPORanges && sdm == "PO" && !firstPlotInRow ) {
      axis(4, labels=TRUE, col=pointsColPO)
    } else if ( is.null(ylimVals) && diffPORanges && sdm == "PO" && firstPlotInRow) {
      warning("No code for correct axes when PO is the first SDM and a different range is required!")
    }
    
    # Do SDM labels at top x-axes?
    if ( columnTitles ) 
      axis(3, (numScenarios+1)/2.0, lwd.ticks = 0, labels=paste(sdm,"SDM"), padj=1.5)

    # Do x-axis for this row's column?
    if ( ! is.null(xAxisLabels) ) axis(1, 1:numScenarios, xAxisLabels, las=2)
  }

}  

#-----------------------------------------------------------------------------------------

plotWarningErrorSummary <- function(scenariosObj, retLst) {
  
  # Simple count of warnings and errors per species x scenario combination.  Plots to console!
  # Visual version of numSuccessfulRunsScenarios function output (I think?).
  
  # numbers of things
  numScenarios <- scenariosObj$numScenarios
  namesSpecies <- retLst$simObj$namesSpecies
  numSpecies <- length(namesSpecies)

  # Get number of warnings per species for each scenario.
  y <- matrix(0,nrow=numScenarios, ncol=numSpecies, dimnames=list(1:numScenarios,namesSpecies))
  for ( i in 1:numScenarios) {
    tmp <- table(retLst$resLstAll[[i]]$warnings$species)
    y[i,names(tmp)] <- tmp
    tmp <- table(retLst$resLstAll[[i]]$errors$species)
    if ( ! is.null(names(tmp)) ) y[i,names(tmp)] <- y[i,names(tmp)] + tmp
  }
  
  # Plot ...
  plotCol <- rainbow(numScenarios)
  plot(c(1,numSpecies), range(y, na.rm=TRUE), type="n", xaxt="n", xlab="", ylab="num warnings + errors")
  for ( i in 1:numScenarios ) {
    points(1:numSpecies, y[i,], col=plotCol[i], pch=i)
  }
  legend("topleft", scenariosObj$namesScenarios, col=plotCol, pch=1:numScenarios)
  axis(1, 1:numSpecies, namesSpecies, las=2)
  title("Number of warnings and errors per scenario for each species", 
        sub=paste0("numRuns=", scenariosObj$numRuns))
  abline(h=0, lty="dashed")

}

#-----------------------------------------------------------------------------------------

makeScenarioStats <- function(scenariosObj, resLstAll, whichStats = c(1,2,3), 
                              numCoresUse = max(1,detectCores()-1), 
                              whichSDMs = resLstAll[[1]]$namesValidSDM, 
                              whichScenarios = 1:scenariosObj$numScenarios, 
                              useContrasts = FALSE) {
  
  # Make the scenario statistics (per species x scenario x sdm x coeff) for the expected number of 
  # values per cell (i.e. look at the accuracy and precision of the estimated lambda in each cell).
  # FYI: loads each scenario's cellObj and simObj from the data already saved in each scenario's
  #      directory.
  # 
  # Arguments ...
  # whichStats:     which statistics/metrics to be calculated
  #                   1 = difference between the true and estimated coefficients 
  #                   2 = correlation between true and estimated lambda[cell_i]
  #                   3 = difference between true and estimated for total expected number 
  #                       of species (i.e. overall abundance for whole domain).  
  #                       This is sum_i(lambda[cell_i] * area[cell_i])!
  # whichSDMs:      for which SDMs are these to be calculated?
  # whichScenarios: for which scenarios are these to be calculated?
  # useContrasts:   whether or not to use true zeta_kg as zeta_kg (FALSE) or the contrast 
  #                 treatment version of zeta_kg - zeta_k1 for g > 1 (TRUE) which will be 
  #                 closer to estimated value.  Also has contrast treatment values for 
  #                 alpha, and gamma but these are dependent on which SDM stats are being
  #                 calculated for.  Only effects statistic 1 (i.e. intercepts of model)!  
  #                 Other statistics use true lambda in cells object.
  
  # Names and numbers of things.  
  numRuns <- scenariosObj$numRuns
  namesScenarios <- scenariosObj$namesScenarios
  prettyNamesScenarios <- scenariosObj$prettyNamesScenarios
  namesSpecies <- resLstAll[[1]]$namesSpecies
  numSpecies <- length(namesSpecies)
  namesValidSDMs <- resLstAll[[1]]$namesValidSDM
  validSDMs <- resLstAll[[1]]$validSDM        # Assumes all scenarios use the same SDMs!
  validStats <- c(1,2,3)
  numGears <- resLstAll[[1]]$numGears         # Assumes all scenarios use the same numGears!
    
  # Which SDMs results to use.
  numRequestedSDMs <- length(whichSDMs)
  if ( ! all(whichSDMs %in% namesValidSDMs) ) {
    stop("Unrecognised SDM requested in argument 'whichSDMs'.")
  } else if ( numRequestedSDMs == 0 ) {
    stop("No SDM has been specified in the 'whichSDMs' argument.")
  }
  
  # Check statistics requested are valid.
  numStats <- length(whichStats)
  if ( numStats == 0 ) {
    stop("No statistics have been requested in the 'whichStats' argument.")
  } else if ( ! all(whichStats %in% validStats) ) {
    stop("Unrecognised statistic requested in argument 'whichStats'.")
  }
  
  # Which scenarios are we interested in?
  numRequestedScenarios <- length(whichScenarios)
  if ( ! all(whichScenarios %in% (1:(scenariosObj$numScenarios))) ) {
    stop("Unrecognised scenario index requested in 'whichScenarios' argument.")
  } else if ( numRequestedScenarios == 0 ) {
    stop("No scenarios have been requested in 'whichScenarios' argument.")
  }
  
  # All true coefficients ...  
  sdm <- validSDMs[1] 
  namesBetas <- rownames(resLstAll[[1]][[sdm]]$coeffs$beta)
  numBetas <- length(namesBetas)    # includes intercept alpha
  namesZetas <- paste0("zeta",1:numGears)
  numZetas <- length(namesZetas)    
  namesGamma <- "gamma"
  numGamma <- 1                    
  namesDeltas <- rownames(resLstAll[[1]][[sdm]]$coeffs$delta)
  numDeltas <- length(namesDeltas)  # not species specific
  namesAllCoeffs <- c(namesBetas, namesZetas, namesGamma, namesDeltas)
  numAllCoeffs <- length(namesAllCoeffs)
  allTrueCoeffs <- array(data = NA, 
                         dim = c(numAllCoeffs, numSpecies, numRequestedSDMs),
                         dimnames = list(namesAllCoeffs, namesSpecies, whichSDMs))
#  allTrueCoeffs <- matrix(nrow = numAllCoeffs, ncol = numSpecies, 
#                          dimnames = list(namesAllCoeffs, namesSpecies))
  
  # Initialise the return value.
  stats <- vector("list", numStats)
  for ( istat in 1:numStats ) {
    thisStat <- whichStats[istat]
    if ( thisStat == 1 ) {
      # Might need extra dimensions for requested coefficients.
      tmp <- array(dim=c(numRequestedScenarios, numSpecies, numRequestedSDMs, numAllCoeffs), 
                   dimnames=list(prettyNamesScenarios[whichScenarios], namesSpecies, 
                                 whichSDMs, namesAllCoeffs))
      stats[[istat]] <- list(avg=tmp, sd=tmp)
    } else {
      tmp <- array(dim=c(numRequestedScenarios, numSpecies, numRequestedSDMs, 1), 
                   dimnames=list(prettyNamesScenarios[whichScenarios], namesSpecies, whichSDMs, 1))
      stats[[istat]] <- list(avg=tmp, sd=tmp)
    }
    stats[[istat]]$usedContrasts <- useContrasts
  }
  
  # Calculate each statistic's mean and sd (across numRuns) for each sdm x species x scenario.
  for ( reqSc in 1:numRequestedScenarios ) {
    # What is the actual scenario number for the reqSc^th requested scenario?
    sc <- whichScenarios[reqSc]
    nameScenario <- namesScenarios[sc]
    resLstScenario <- resLstAll[[sc]]
    
    # Start timing.
    message("Begin statistics calculation for scenario ", nameScenario)
    timeBegin <- Sys.time()
    message(paste0("    Begin: ", timeBegin))
    
    # Load required data (i.e. cellObj and simObj which are different for each scenario)
    load(paste0(scenariosObj$scenariosDir, "/", nameScenario,"/DataDumps/DataDump-AllRuns.RData"))
    

    # True coefficient values for this scenario.
    beta <- as.matrix(simObj$beta[ ,namesSpecies])
    alpha <- beta[1, ]
    zeta <- as.matrix(simObj$zeta[ ,namesSpecies])
    gamma <- matrix(simObj$gamma[namesSpecies],nrow=1,ncol=numSpecies)
    rowsZeta <- (numBetas + 1):(numBetas + numZetas)
    rowsGamma <- numBetas + numZetas + 1
    rowsDelta <- (numBetas + numZetas + numGamma + 1):numAllCoeffs
    
    for ( reqSDM in 1:numRequestedSDMs ) {    
      # Store true coefficient values for scenario.
      allTrueCoeffs[1:(numBetas+numZetas+numGamma),namesSpecies,reqSDM] <- rbind(beta,zeta,gamma)
      allTrueCoeffs[rowsDelta,1,reqSDM] <- simObj$delta
      
      # Change "true" coefficient value to account for the contrast treatment, if requested.
      if ( useContrasts && numZetas > 1 ) {
        nameSDM <- whichSDMs[reqSDM]
        
        # Calculate the sum, across g, of the zeta_kg ...
        logMeanExpZeta <- log(apply(exp(zeta[ ,namesSpecies]),2,mean))
        
        # alpha_k ...
        if ( whichSDMs[reqSDM] == "MsPP") {
          allTrueCoeffs[1,namesSpecies,reqSDM] <- as.matrix(alpha + logMeanExpZeta)
        } else {
          allTrueCoeffs[1,namesSpecies,reqSDM] <- as.matrix(alpha + zeta[1, ])
        }
        
        # gamma_k ...
        if ( nameSDM == "Gear" ) {
          allTrueCoeffs[rowsGamma,namesSpecies,reqSDM] <- as.matrix(gamma + logMeanExpZeta
                                                                    - zeta[1, ])
        } else if ( nameSDM == "MsPP"){
          allTrueCoeffs[rowsGamma,namesSpecies,reqSDM] <- as.matrix(gamma)
        }
        
        # zeta_kg for g > 1
        if ( nameSDM %in% c("PA","Gear") ) {
          # make these alpha_k + zeta_kg for g > 1 to get around problem of values being close to zero.
          for ( g in 2:numZetas ) {
            # True zeta values ...
            allTrueCoeffs[rowsZeta[g],namesSpecies,reqSDM] <- alpha + zeta[g, ]
            
            # Estimates of zeta values ...
            indNameSDM <- which(nameSDM == namesValidSDMs)
            validSDM <- validSDMs[indNameSDM]
            resLstScenario[[validSDM]]$coeffs$zeta[g-1,namesSpecies, ] <- 
                                resLstScenario[[validSDM]]$coeffs$beta[1,namesSpecies, ]
               + resLstScenario[[validSDM]]$coeffs$zeta[g-1,namesSpecies, ]
          } 
        }
      } 
    } # reqSDM in 1:numRequestedSDMs
    
    # Run parallelised version of statistics calculations (sdm x run x stat loops inside here)
    # parRes <- internalParallelScenarioStats(namesSpecies, resLstScenario, whichSDMs, 
    #                                         whichStats, allTrueCoeffs, cellObj, 
    #                                         simObj$lambda.formula)
    
    parRes <- mclapply(namesSpecies, internalParallelScenarioStats, resLstScenario,
                       whichSDMs, whichStats, allTrueCoeffs, cellObj, simObj$lambda.formula,
                       mc.cores = numCoresUse)
    
    # Collect results into return value.
    for (k in 1:numSpecies ) {
      # Was there an error?
      if ( is.null(parRes[[k]]) ) {
        message(paste0("Error for scenario ", nameScenario, ":\n", 
                       "Results are not available as core returned NULL."))

      } else if ( parRes[[k]]$numErrors > 0 ) {
        message(paste0("Error for scenario ", nameScenario, ":"))
        for ( err in parRes[[k]]$numErrors ) {
          message(parRes[[k]]$errorMsgs[err])
        }
      } else {
        species <- parRes[[k]]$species
        for ( istat in 1:numStats) {
          thisStat <- whichStats[istat]
          if ( thisStat == 1 ) {
            for ( j in 1:numAllCoeffs ) {
              stats[[istat]]$avg[reqSc,species, ,j] <- parRes[[k]]$avg[istat, ,j]
              stats[[istat]]$sd[reqSc,species, ,j] <- parRes[[k]]$sd[istat, ,j]
            }
          } else {
            stats[[istat]]$avg[reqSc,species, ,1] <- parRes[[k]]$avg[istat, ,1]
            stats[[istat]]$sd[reqSc,species, ,1] <- parRes[[k]]$sd[istat, ,1]
          }
        }
      }
      
      # Were there warnings?
      if ( parRes[[k]]$numWarnings > 0 ) {
        message(paste0("Warning for scenario ", nameScenario,":"))
        for ( wrn in parRes[[k]]$numWarnings ) {
          message(parRes[[k]]$warningMsgs[wrn])
        }
      }
    }    
    timeEnd <- Sys.time()
    message(paste0("    End:   ", timeEnd))
  }
  
  # Return value.
  return(stats)
  
}

#-----------------------------------------------------------------------------------------

plotScenariosSummaryLambda <- function(stat, nameStat="", whichCoeffs=1:dim(stat$avg)[[4]],
                                       ylimVals=NULL, plotSDMs=dimnames(stat$avg)[[3]], 
                                       vioPlots=TRUE, absMeans=TRUE, plotDevice="RStudioGD", 
                                       plotDir=getwd(), plotUnits="cm", plotHeight=9, 
                                       plotWidth=12.25, plotRes=600, doTitle=FALSE,
                                       xAxisLabels=NULL, xAxisTitle=NULL, xAxisLabelStyle=2,
                                       bigNumThreshold = 10, horizontalLines=c(0,0),
                                       diffPORange=c(FALSE,FALSE), outlierLabels=FALSE,
                                       whichScenarios=1:dim(stat$avg)[[1]]) {

  # Plot the scenarios' summary for a single overall statistic (i.e. look at the accuracy/avg
  # and precision/sd of the given statistic). Plotting function for "makeScenarioStats" but
  # call this function once for *each* statistic made in "makeScenarioStats".

  # Names and numbers of things.  
  dimnamesStat <- dimnames(stat$avg)
  namesScenarios <- dimnamesStat[[1]][whichScenarios]
  numScenarios <- length(namesScenarios)
  namesSpecies <- dimnamesStat[[2]]
  numSpecies <- length(namesSpecies)
  namesValidSDMs <- dimnamesStat[[3]]
  numValidSDMs <- length(namesValidSDMs)
  namesForthDim <- dimnamesStat[[4]]            # May be coefficients, may not be!
  numForthDim <- length(namesForthDim)
  
  # Check whichCoeffs are within valid range.
  numReqCoeffs <- length(whichCoeffs)
  if ( numReqCoeffs == 0 && numForthDim > 1 ) {
    stop("No coefficients have been requested in the 'whichCoeffs' argument.")
  } else if ( ! all(whichCoeffs %in% 1:numForthDim) ) {
    stop("Unrecognised coefficient requested in argument 'whichCoeffs'.")
  }

  # Were x-axis labels provided (i.e. pretty names for scenarios)?
  if ( is.null(xAxisLabels) ) {
    xAxisLabels <- namesScenarios
  } else if ( length(xAxisLabels) == numScenarios && 
              length(unique(xAxisLabels)) == numScenarios ) {
    # use provided labels.
  } else {
    stop("Invalid x-axis labels given in argument 'xAxisLabels'.")
  }

  # Figure out the number of SDMs to plot.
  numPlotSDMs <- length(plotSDMs)
  if ( ! all(plotSDMs %in% namesValidSDMs) ) {
    stop("Unrecognised SDM requested in plots.")
  }
  
  # Plot directory for scenario level plots.
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Are the absolute value of the means to be used?
  if ( absMeans ) stat$avg <- abs(stat$avg)
  
  # Threshold setting for plotting extremely big numbers (or big negative numbers)!
  for ( sdm in plotSDMs ) {
    for ( sc in whichScenarios ) {
      for ( d4 in whichCoeffs) {
        # For this sdm and this scenario, the species accuracy statistics are:
        thisStat <- stat$avg[sc, ,sdm,d4]
        indBigNum <- which(abs(thisStat) > bigNumThreshold)
        if ( length(indBigNum) > 0 ) {
          signBigNum <- sign(thisStat[indBigNum])
          stat$avg[sc,indBigNum,sdm,d4] <- signBigNum * bigNumThreshold
        }
        
        # For this sdm and this scenario, the species precision statistics are:
        thisStat <- as.vector(stat$sd[sc, ,sdm,d4])
        indBigNum <- which((abs(thisStat) > bigNumThreshold))
        if ( length(indBigNum) > 0) {
          signBigNum <- sign(thisStat[indBigNum])
          stat$sd[sc,indBigNum,sdm,d4] <- signBigNum * bigNumThreshold
        }
        
        # For NaN which comes out of the sd function (eg. sd(c(1:5,Inf)) = NaN)
        indNaN <- which(is.nan(thisStat))
        if ( length(indNaN) > 0) stat$sd[sc,indNaN,sdm,d4] <- NA
      }
    }
  } 
  
  
  # Is the y-axis truncated (assumes ylimVals are less than actual range)
  if ( is.null(ylimVals) ) {
    fileTrunc <- ""
    titleTrunc <- ""
  } else {
    fileTrunc <- "-trunc"
    titleTrunc <- " - truncated"
  }
  
  # Which device are we printing to?
  if ( plotDevice == "RStudioGD" ) {
    # plot to the R studio plot window.
    plotToFile <- FALSE
  } else {
    fileName <- makeFileName(paste0(nameStat, fileTrunc), plotDir, plotDevice)
    argList <- list(filename = fileName, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
    plotToFile <- TRUE
  }
  
  # Divide plotting into the number of columns required.
  omaPar <- par("mar")         # use default single plot per page margins
  marPar <- c(0.0,0.0,1.3,0.0) # small margin on top to separate rows and allow for column titles.
  if ( is.null(xAxisTitle) ) {
    # No x-axis title, change margins to gain plotting area.
    omaPar[1] <- omaPar[1] - 1.0  # leave space for x-axis ticks and labels.  
  }
  if ( ! doTitle ) {
    # No main title, change margins to gain plotting area.
    omaPar[3] <- 0.0  # NB: there is already a bit at top due to mar settings.
  }
  opar <- par(mfrow=c(2,numPlotSDMs), oma = omaPar, mar = marPar, cex=0.66)
  
  # Colours for plotting.
  pointsCol <- "black"
  violinCol <- "grey"
  pointsColPO <- colourBlindRGB("blue")
  violinColPO <- colourBlindRGB("skyBlue")
  
  # Y-axis range ...
  yaxisRange <- array(dim=c(2,2,numPlotSDMs), 
                      dimnames = list(c("mean","sd"), c("min","max"), plotSDMs))  
  if ( is.null(ylimVals) ) {
    # What is the position of the "PO" SDM in the SDMs to be plotted?
    indPOSDM <- which(plotSDMs %in% "PO")
    indOtherSDMs <- 1:numPlotSDMs
    if ( length(indPOSDM) > 0) indOtherSDMs <- indOtherSDMs[-indPOSDM]
    
    # Plot the PO SDM with a different y-axis range for the first row?
    if ( diffPORange[1] && "PO" %in% plotSDMs ) {
      # Yes (and it is one of the plots!!!)
      yaxisRange[1,1,indPOSDM] <- min(stat$avg[whichScenarios, ,plotSDMs[indPOSDM],whichCoeffs])
      yaxisRange[1,2,indPOSDM] <- max(stat$avg[whichScenarios, ,plotSDMs[indPOSDM],whichCoeffs])
      if ( length(indOtherSDMs) > 0 ) {
        # There are other SDMs to be plotted.
        yaxisRange[1,1,indOtherSDMs] <- rep(min(stat$avg[whichScenarios, ,plotSDMs[indOtherSDMs],whichCoeffs]), 
                                            times=length(indOtherSDMs))
        yaxisRange[1,2,indOtherSDMs] <- rep(max(stat$avg[whichScenarios, ,plotSDMs[indOtherSDMs],whichCoeffs]), 
                                            times=length(indOtherSDMs))
      }
    } else {
      # No, use whatever the values suggest for the y axis range (same for all SDMs)
      yaxisRange[1,1, ] <- rep(min(stat$avg[whichScenarios, ,plotSDMs,whichCoeffs]), times=numPlotSDMs)
      yaxisRange[1,2, ] <- rep(max(stat$avg[whichScenarios, ,plotSDMs,whichCoeffs]), times=numPlotSDMs)
    }
    
    # Plot the PO SDM with a different y-axis range for the second row?
    yaxisRange[2,1, ] <- rep(0, times=numPlotSDMs)
    if ( diffPORange[2] && "PO" %in% plotSDMs ) {
      # Yes (and it is one of the plots!!!)
      yaxisRange[2,2,indPOSDM] <- max(stat$sd[whichScenarios, ,indPOSDM,whichCoeffs])
      if ( length(indOtherSDMs) > 0 ) {
        # There are other SDMs to be plotted.
        yaxisRange[2,2,indOtherSDMs] <- rep(max(stat$sd[whichScenarios, ,indOtherSDMs,whichCoeffs]), 
                                            times=length(indOtherSDMs))
      }
    } else {
      # No, use whatever the values suggest for the y axis range (same for all SDMs)
      yaxisRange[2,2, ] <- rep(max(stat$sd[whichScenarios, ,plotSDMs,whichCoeffs], na.rm=TRUE), times=numPlotSDMs)
    }
    
  } else if ( length(ylimVals) == 4 ) {
    # A different range for mean and sd, eg. c(min mean, max mean, min sd, max sd).
    yaxisRange[1,1, ] <- ylimVals[1]
    yaxisRange[1,2, ] <- ylimVals[2]
    yaxisRange[2,1, ] <- ylimVals[3]
    yaxisRange[2,2, ] <- ylimVals[4]
  } else if ( length(ylimVals) == 2 ){
    # Same for mean and sd, eg. c(min, max)
    yaxisRange[1,1, ] <- ylimVals[1]
    yaxisRange[1,2, ] <- ylimVals[2]
    yaxisRange[2, , ] <- yaxisRange[1, , ]    
  } else {
    stop("Unrecognised version of 'ylimVals' argument.")
  }
  
  # Start plots.
  rowTitles <- c("mean", "sd")
  rowDataNames <- c("avg","sd")
  for ( row in 1:2 ) {
    # What is the row title.
    rowTitle <- rowTitles[row]
    
    # When are these things added?
    if ( row == 1 ) {
      columnTitles <- TRUE
      xAxisTicks <- FALSE
    } else {
      columnTitles <- FALSE
      xAxisTicks <- TRUE
    }
    firstPlotInRow <- TRUE
    
    # Calculate widths for violin plots in row.
    if ( vioPlots ) {
      # Get the maximum density for each SDM x scenario combination (i.e. each violin).
      maxDensity <- matrix(nrow=numScenarios, ncol=numPlotSDMs, dimnames=list(namesScenarios, plotSDMs))
      for ( sdm in plotSDMs ) {
        
        for ( sc in 1:numScenarios ) {
          # For this sdm and this scenario, the species statistics are:
          thisStat <- as.vector(stat[[rowDataNames[row]]][whichScenarios[sc], ,sdm, whichCoeffs])

          # Check for NA values.
          indNA <- which(is.na(thisStat))
          if ( length(indNA) > 0 ) {
            warning(paste0(length(indNA), " NA values in ", rowDataNames[row], 
                           " for SDM = ", sdm, 
                           " and scenario = ", namesScenarios[sc]))
            thisStat <- thisStat[-indNA]
            #if ( outlierLabels ) thisLabels <- thisLabels[-indNA]
          }
          
          # What is the maximum of the density for each scenario x sdm.
          maxDensity[sc,sdm] <- max(density(thisStat)$y)
        }
      }
      
      # Scale maximum densities so that widths are between [0,1]
      maxMaxDensity <- max(maxDensity[])
      widthsViolins <- maxDensity/maxMaxDensity
    }
    
    # Do a plot for each SDM.
    for ( sdm in plotSDMs ) {
      # Plot statistics for each SDM (i.e. column) for this row.
      plot(c(0.5,numScenarios+0.5), yaxisRange[row, ,sdm], type="n", xaxs="i", xaxt="n", xlab="", 
           yaxt="n", ylab="")
      for ( sc in 1:numScenarios ) {
        # All species' stats at this scenario and SDM!
        thisStats <- as.vector(stat[[rowDataNames[row]]][whichScenarios[sc], ,sdm,whichCoeffs])
        if ( outlierLabels ) thisLabels <- rep(namesSpecies, numReqCoeffs)
        indNA <- which(is.na(thisStats))
        if ( length(indNA) > 0 ) {
          # warning(paste0(length(indNA), " NA values in ", rowDataNames[row], 
          #                " for SDM = ", sdm, 
          #                " and scenario = ", namesScenarios[sc]))
          thisStats <- thisStats[-indNA]
          #if ( outlierLabels ) thisLabels <- thisLabels[-indNA]
        }
        
        # Plot for this sdm x scenario.
        # If PO is to be a different range, do different y-axis for PO.
        if ( is.null(ylimVals) && diffPORange[row] && sdm == "PO" ) {
          if ( vioPlots ) vioplot(thisStats, col=violinColPO, border=violinColPO, add=TRUE,
                                  drawRect = FALSE, at=sc, wex=widthsViolins[sc,sdm])
          axis(4, labels=TRUE, col=pointsColPO)
          stripchart(thisStats, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsColPO)
        } else {
          if ( vioPlots ) vioplot(thisStats, col=violinCol, border=violinCol, add=TRUE,
                                  drawRect = FALSE, at=sc, wex=widthsViolins[sc,sdm])
          stripchart(thisStats, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsCol)
        }

        
        # Labels on outliers?
        if ( outlierLabels ) {
          # Use the interquartile range method to identify outliers.
          quartiles <- quantile(thisStats, probs=c(0.25,0.5,0.75), type = 8)
          intQuartRange <- quartiles[3] - quartiles[1]
          notOutlierRange <- c(quartiles[1] - 1.5*intQuartRange,
                               quartiles[3] + 1.5*intQuartRange)
          indOutliers <- which(thisStats < notOutlierRange[1] | thisStats > notOutlierRange[2])
          numOutliers <- length(indOutliers)
          if ( numOutliers > 0 ) {
            # FYI: cex=0.3 will only be readable if the plot resolution (i.e. pixels per inch)
            #      is high enough (at least 600ppi).
            text(sc, thisStats[indOutliers], labels=thisLabels[indOutliers], pos=4, cex=0.3,
                 offset=0.2)
          }

          # Median and boundaries between outliers and "normal" values.
          points(sc, quartiles[2], pch="_", col="green")
          points(sc, notOutlierRange[1], pch="_", col="blue")
          points(sc, notOutlierRange[2], pch="_", col="blue")
        }
      }
      abline(h=horizontalLines[row], lty="dotted")
      #abline(h=0, lty="dotted")
      #if ( row == 1 ) abline(h=1, lty="dotted")
  
      # If this is the first column, do y-axis.
      if ( firstPlotInRow ) {
        axis(2, labels=TRUE)
        midyRange <- ((yaxisRange[row,2,sdm] - yaxisRange[row,1,sdm])/2.0) + yaxisRange[row,1,sdm]
        axis(2, at=midyRange, labels=rowTitle, tick=FALSE, padj=-2.0)
        firstPlotInRow <- FALSE
      }
      
      # If this is the PO column (this will cause all sorts of problems if PO is the first column!)
      # if ( is.null(ylimVals) && diffPORanges && sdm == "PO" && !firstPlotInRow ) {
      #   axis(4, labels=TRUE, col=pointsColPO)
      # } else if ( is.null(ylimVals) && diffPORanges && sdm == "PO" && firstPlotInRow) {
      #   warning("No code for correct axes when PO is the first SDM and a different range is required!")
      # }
      
      # Do SDM labels at top x-axes?
      if ( columnTitles ) 
        axis(3, (numScenarios+1)/2.0, lwd.ticks = 0, labels=paste(sdm,"SDM"), padj=1.0)
      
      # Do x-axis for this row's column?
      if ( xAxisTicks ) {
        axis(1, 1:numScenarios, xAxisLabels, las=xAxisLabelStyle)
      }
    }
  }

  # Add a title to the plot.
  if ( doTitle ) 
    title(paste0("Similarity of intensities (", nameStat, 
                 ") for the given scenarios with the given SDMs", titleTrunc), outer=TRUE)
  
  # Add an x-axis title to the plot.
  if ( ! is.null(xAxisTitle) ) title(xlab = xAxisTitle, outer=TRUE, line=min(4,omaPar[1]))
  
  # Turn off the plotting device, if it is a file.
  if ( plotToFile ) dev.off()
  
  # Return plotting to original settings.
  par(opar)

}

#-----------------------------------------------------------------------------------------

internalParallelScenarioStats <- function(species, resLstScenario, 
                                          whichSDMs = resLstScenario$namesValidSDM, 
                                          whichStats = c(1,2,3), 
                                          trueCoeffs = NULL, 
                                          cellObj = NULL, 
                                          lambda.formula = NULL) {

  # Runs part of the scenario stats function in parallel to make it faster, hopefully.
  # No error checking as assume function that calls this has done error checking.
  #
  # Arguments ...
  # species:        the name of the species to run here.
  # resLstScenario: the results from the required scenario (i.e. all the coefficient estimates)
  # whichSDMs:      only work out statistics for these SDMs.
  # whichStats:     a vector giving the stats that are to be performed (by an index number)
  #                     1 = difference between true and estimated coefficients.
  #                     2 = correlation of true to estimated lambda value at each cell
  #                     3 = difference between true and estimated total abundance.
  # trueCoeffs:     an array of the true coefficients for all species (row order: alpha, 
  #                 betas, zetas, gamma and deltas; numAllCoeffs x numSpecies x numSDMs).  
  #                 Note that zeta1 will be included here (unlike estiamtes).  Only needed  
  #                 if whichStats includes 1!
  # cellObj:        the cell object that corresponds to the required scenario.  Only needed
  #                 if whichStats includes 2 or 3.  Should contain trueLambda values,
  #                 cell area and covariate values, for the cells in the domain.
  # lambda.formula: formula from which lambda values are created (as in glm formula).
  #                 Only needed if whichStats includes 2 or 3.

  # Names and numbers of things. 
  numRuns <- resLstScenario$numRuns
  namesValidSDMs <- resLstScenario$namesValidSDM
  validSDMs <- resLstScenario$validSDM        # Assumes all scenarios use the same SDMs!
#  numValidSDMs <- length(validSDMs)           # Assumes namesValidSDMs and validSDMs are same length!
  numStats <- length(whichStats)
  areaCell <- cellObj$areaCell
  validStats <- c(1,2,3)
  numSDMs <- length(whichSDMs)
  
  # Figure out the validSDMs needed from the whichSDMs requested (storage name can be different!)
  indWhichSDMs <- match(whichSDMs, namesValidSDMs)
  whichValidSDMs <- validSDMs[indWhichSDMs] 
  numBetas <- dim(resLstScenario[[whichValidSDMs[1]]]$coeffs$beta)[[1]]
  
  # Check correct arguments are present.
  if ( 2 %in% whichStats || 3 %in% whichStats ) {
    if ( is.null(lambda.formula) || is.null(cellObj) ) 
      stop("Necessary arguments for statistics 2 or 3 are missing.")
  }
  
  # Initialise return value (just uses 1st index in 3rd dimension if stat != 1)
  if ( 1 %in% whichStats ) {
    if ( is.null(trueCoeffs) ) {
      stop("Argument 'trueCoeffs' needs to be provided when whichStats contains 1.")
    } else {
      namesCoeffs <- dimnames(trueCoeffs)[[1]]
      numCoeffs <- length(namesCoeffs)
    }
  } else {
    namesCoeffs <- dimnames(resLstScenario[[whichValidSDMs[1]]]$coeffs$beta)[[1]]
    numCoeffs <- length(namesCoeffs)
  }
  tmp <- array(dim=c(numStats,numSDMs,numCoeffs), dimnames=list(whichStats, whichSDMs, namesCoeffs))
  stats <- list(avg=tmp, sd=tmp, species=species, numErrors=0, errorMsgs=NULL, 
                numWarnings=0, warningMsgs=NULL)

  # True values intensity and total abundance statistics.
  if ( 2 %in% whichStats || 3 %in% whichStats ) speciesTrueLambda <- cellObj$trueLambda[ ,species]

  for ( isdm in 1:numSDMs ) {
    # Name of the SDM in the results list.
    nameSDM <- whichValidSDMs[isdm]
    
    # True values for coefficients.
    if ( 1 %in% whichStats ) speciesTrueCoeffs <- trueCoeffs[ ,species,whichSDMs[isdm]] 
    
    # Calculate the statistics for each run.  
    statsRuns <- array(dim=c(numRuns,numStats,numCoeffs), 
                       dimnames=list(1:numRuns,whichStats,namesCoeffs))

    # Calculate stats
    for ( i in 1:numRuns ) {
      # Get the estimated coefficients for this run.  
      speciesEstCoeffs <- c(resLstScenario[[nameSDM]]$coeffs$beta[ ,species,i],
                            NA,                                                   # Dummy for zeta1.
                            resLstScenario[[nameSDM]]$coeffs$zeta[ ,species,i],
                            resLstScenario[[nameSDM]]$coeffs$gamma[species,i],
                            resLstScenario[[nameSDM]]$coeffs$delta[ ,i])

      # Are these valid coefficients (i.e. has this run produced an estimate)? 
      if ( all(! is.na(speciesEstCoeffs[1:numBetas])) ) {
        # Get the estimated intensity for this run, if necessary.
        if ( 2 %in% whichStats || 3 %in% whichStats )
          speciesEstLambda <- lambda.cell(lambda.formula, speciesEstCoeffs[1:numBetas], cellObj$covars)
          

        # Calculate the statistics for this run.
        for ( istat in 1:numStats ) {
          thisStat <- whichStats[istat]
          if ( thisStat == 1 ) {
            # Difference between true and estimated coefficients.
            statsRuns[i,istat, ] <- (speciesTrueCoeffs - speciesEstCoeffs)/speciesTrueCoeffs
            
          } else if ( thisStat == 2 ) {
            # Correlation between true and this run's estimated lambda.
            statsRuns[i,istat,1] <- cor(speciesTrueLambda, speciesEstLambda)

          } else if ( thisStat == 3 ) {
            # Difference in overall expected abundance between the true and estimated solution.
            # (Scale by true abundance to get percentage of difference per species)
            # NB: areaCell cancels out as it is in both the numerator and denominator!
            trueTotalAbundance <- sum(speciesTrueLambda)
            estTotalAbundance <- sum(speciesEstLambda)
            if ( is.infinite(estTotalAbundance) ) {
              msg <- paste0("  Infinite estimate for total abundance for run ", i, 
                               ", species ", species, ", and SDM ", nameSDM, ".")
              stats$numWarnings <- stats$numWarnings + 1
              stats$warningMsgs <- c(stats$warningMsgs, msg)
            } else {
              statsRuns[i,istat,1] <- (trueTotalAbundance - estTotalAbundance)/trueTotalAbundance
            }
          } else {
            stats$numErrors <- stats$numErrors + 1
            stats$errorMsgs <- c(stats$errorMsgs, 
                                 "Unrecognised statistic requested.  Check argument 'whichStats'.")
            return(stats)
          } 
        }
      } else {
        statsRuns[i, , ] <- rep(NA, numStats*numCoeffs)
      }
    }

    # Summarise across runs (i.e. what is the average performance, what is the variance of this).
    for ( istat in 1:numStats ) {
      thisStat <- whichStats[istat]
      if ( thisStat == 1 ) {
        num3Dim <- numCoeffs
      } else {
        num3Dim <- 1
      }
      for ( id3 in 1:num3Dim ) {  #statsRuns[i,stat,1]
        stats$avg[istat,isdm,id3] <- mean(statsRuns[ ,istat,id3], na.rm=TRUE)
        stats$sd[istat,isdm,id3] <- sd(statsRuns[ ,istat,id3], na.rm = TRUE)
        if ( is.infinite(stats$sd[istat,isdm,id3]) ) {
          msg <- paste0("  Infinite value for SD of statistic ", whichStats[istat], 
                        ", species ", species, ", and SDM ", nameSDM, ".")
          stats$numWarnings <- stats$numWarnings + 1
          stats$warningMsgs <- c(stats$warningMsgs, msg)
        }
        #print(id3)
      }
      #print(istat)
    }
    #print(isdm)
  }

  # Return value.
  #print("hello")
  return(stats)
  
}

#-----------------------------------------------------------------------------------------

saveScenariosSummaryStats <- function(stats, retLst, 
                                      namesStats=c("Coefficients","Response", "Total Abundance"),
                                      scenariosDir = retLst$scenariosObj$scenariosDir, 
                                      nameResultsFile="Results") {
  
  # Add the summary stats to the saved results.
  
  # Get the number of statistics
  numStats <- length(stats)
  if ( numStats != length(namesStats) ) stop("The 'namesStats' argument is incorrect length.")
  
  # Directory where the data was saved ...  
  savedDir <- paste0(scenariosDir,"/Data")
  if ( ! dir.exists(savedDir) ) stop("Unable to save data as scenarios directory doesn't exist.")
  
  # File where the data was saved ...
  savedFile <- makeFileName(nameResultsFile, savedDir,  "RData")
  
  # Save objects.
  statsObj <- list(numStats=numStats, namesStats=namesStats, stats=stats)
  plotObj <- retLst$plotObj
  simObj <- retLst$simObj
  resLstAll <- retLst$resLstAll
  scenariosObj <- retLst$scenariosObj
  save(plotObj, simObj, resLstAll, scenariosObj, statsObj, file=savedFile)
  
  # Return value.
  retLst$statsObj <- statsObj
  return(retLst)
  
}

#-----------------------------------------------------------------------------------------

plotExperimentComparisons <- function(whichExperiments, whichStat=1, whichCoeffs=1, whichSpecies=NULL,
                                      plotSDMs=c("PA","MsPP","Gear"), ylimVals=c(0,1,0,1), 
                                      vioPlots=TRUE, absMeans=TRUE, plotDevice="RStudioGD", 
                                      plotDir=getwd(), fileName="comparison", plotUnits="cm", 
                                      plotHeight=20, plotWidth=15.75, plotRes=600, plotTitle=NULL,
                                      xAxisLabels=NULL, xAxisTitle=NULL, xAxisLabelStyle=2,
                                      xAxisReverse=FALSE,  
                                      columnHeadings=paste("experiment",1:length(whichExperiments)),
                                      horizontalLines=c(0,0), diffPORange=NULL, medianLines=TRUE,
                                      accuracyPrecision = "both") {
  
  # Plot a single statistic for a group of experiments (e.g. plot correlation for all bias 
  # change experiments).   
  #
  # Arguments ...
  # whichExperiments: a string vector containing directory names, one for each experiment.
  # whichStat:        an integer scalar that indicates the statistic within each experiment 
  #                   that we are interested in.  See 'makeScenarioStats' for valid values.
  # whichCoeffs:      an integer vector which indicates which of the coefficients should 
  #                   be included if whichStat = 1.
  # whichSpecies:     which species to plot.  A value of NULL will use all species as 
  #                   specified in the first experiment.
  # plotSDMs:         which SDMs from the experiment to plot (and the order).  
  # ylimVals:         a vector of length 4 that gives a y-axis range for mean and sd, eg. 
  #                   c(min mean, max mean, min sd, max sd).
  # xAxisTitle:       can be a single string (in which case it will be placed in the centre
  #                   of the outer plot margins) or a vector of strings, one for each experiment.
  # xAxisReverse:     reverse the positions of the scenarios along the x-axis for all experiments
  #                   ( = TRUE) or for specified experiments (e.g. = c(TRUE, FALSE, FALSE))
  # columnHeadings:   a vector of strings, one for each experiment, to use as column headings.
  #                   Set to "NULL" if no column headings are required.
  # horizontalLines:  a vector of length 2 that given the y values at which horizontal lines
  #                   will be drawn for each accuracy (mean) and each precision (sd) plot.
  # diffPORange:      a vector of length 2 that gives a different y-axis range for PO SDM 
  #                   mean plots.  (Set to NULL to use same y-axis range as other SDMs)
  # medianLines:      True to include a line to indicate the median value per scenario of
  #                   each plot, false to not include.
  # accurayPrecision: whether to plot accuracy or precision or both (these are rows of plots).
  #                   Acceptable values are "accuracy", "precision", "both".
  
  # Get numbers of things to set up plotting layout.
  numExperiments <- length(whichExperiments)
  numPlotSDMs <- length(plotSDMs)
  substituteInfValue <- 1.0e+25
  
  # Check accuracy or precision or both.
  if ( accuracyPrecision == "accuracy" ) {
    doAccuracy <- TRUE
    doPrecision <- FALSE
  } else if ( accuracyPrecision == "precision" ) {
    doAccuracy <- FALSE
    doPrecision <- TRUE
  } else if ( accuracyPrecision == "both" ) {
    doAccuracy <- TRUE
    doPrecision <- TRUE
  } else {
    doAccuracy <- FALSE
    doPrecision <- FALSE
    stop("No valid comparisons have been requested, check argument 'accuracyPrecision'.")
  }
    
  # Colours for plotting.
  pointsCol <- "black"
  violinCol <- "grey"
  pointsColPO <- colourBlindRGB("blue")
  violinColPO <- colourBlindRGB("skyBlue")
  medianCol <- colourBlindRGB("orange")
  
  # Plot directory for experiment comparison level plots.
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }

  # Which device are we printing to?
  if ( plotDevice == "RStudioGD" ) {
    # plot to the R studio plot window.
    plotToFile <- FALSE
  } else {
    fileName <- makeFileName(fileName, plotDir, plotDevice)
    argList <- list(filename = fileName, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
    plotToFile <- TRUE
  }
  
  # Divide plotting into the number of rows and columns required.
  omaPar <- par("mar")         # use default single plot per page margins
  marPar <- c(0.0,0.5,1.3,0.0) # small margin on top to separate rows and allow for column headings.
  if ( is.null(xAxisTitle) ) {
    # No x-axis title, change margins to gain plotting area.
    omaPar[1] <- omaPar[1] - 1.0  # leave space for x-axis ticks and labels.  
  }
  if ( is.null(plotTitle) ) {
    # No main title, change margins to gain plotting area.
    omaPar[3] <- 0.0  # NB: there is already a bit at top due to mar settings.
  }
  if ( doAccuracy && doPrecision ) {
    # Both accuracy and precision is required.
    opar <- par(mfcol=c(numPlotSDMs*2,numExperiments), oma = omaPar, mar = marPar, cex=0.66)
  } else {
    # Either accuracy or precision is required but not both.
    opar <- par(mfcol=c(numPlotSDMs,numExperiments), oma = omaPar, mar = marPar, cex=0.66)
  }
  
  if ( length(xAxisReverse) == 1 ) {
    xAxisReverseUse <- rep(xAxisReverse, numExperiments)
  } else if ( length(xAxisReverse) == numExperiments ) {
    xAxisReverseUse <- xAxisReverse
  } else {
    stop("Argument 'xAxisReverse' is an unexpected length.")
  }
  
  # Load experiments (FYI: necessary to load and store so that violin width can be calculated)
  allExperimentsInfo <- vector(mode="list", length = numExperiments)
  maxNumScenarios <- 0
  for ( ex in 1:numExperiments ) {
    # Load in this experiment's results
    exDir <- whichExperiments[ex]
    retLst <- loadScenariosResults(exDir)
    
    # Check requested SDMs are contained within this experiment.
    if ( ! all(plotSDMs %in% retLst$resLstAll[[1]]$namesValidSDM) ) {
      stop("Unrecognised SDM requested from experiment '", exDir, "'.")
    }
    
    # Check requested species are contained within this experiment.
    exSpecies <- retLst$simObj$namesSpecies
    if ( is.null(whichSpecies) ) {
      # Use the species in the first experiment, other experiments will be checked against these.
      whichSpecies <- exSpecies
    } else if ( ! all(whichSpecies %in% exSpecies) ) {
      stop("Unrecognised species requested from experiment '", exDir, "'.")
    }

    # Store this experiment's information.
    allExperimentsInfo[[ex]] <- list(numScenarios = retLst$scenariosObj$numScenarios,
                                     namesScenarios = retLst$scenariosObj$prettyNamesScenarios,
                                     exDir = exDir, whichStat = whichStat, whichSpecies = whichSpecies,
                                     thisStatAvg=retLst$statsObj$stats[[whichStat]]$avg[ ,whichSpecies, , ,drop=FALSE], 
                                     thisStatSD=retLst$statsObj$stats[[whichStat]]$sd[ ,whichSpecies, , ,drop=FALSE])
    if ( absMeans ) allExperimentsInfo[[ex]]$thisStatAvg <- abs(allExperimentsInfo[[ex]]$thisStatAvg)
    if ( xAxisReverseUse[ex] ) {
      # Reverse the scenarios on the x-axis for this experiment.
      allExperimentsInfo[[ex]]$namesScenarios <- rev(allExperimentsInfo[[ex]]$namesScenarios)
      numScenarios <- allExperimentsInfo[[ex]]$numScenarios
      allExperimentsInfo[[ex]]$thisStatAvg <- allExperimentsInfo[[ex]]$thisStatAvg[numScenarios:1, , , ,drop=FALSE]
      allExperimentsInfo[[ex]]$thisStatSD <- allExperimentsInfo[[ex]]$thisStatSD[numScenarios:1, , , ,drop=FALSE]
    }
    
    # What is the maximum number of scenarios, so far?
    maxNumScenarios <- max(maxNumScenarios, allExperimentsInfo[[ex]]$numScenarios)
    
  }  # experiments (or plot columns).

  # Calculate widths for violin plots for means and SDs, if necessary.
  if ( vioPlots ) {
    # Get the maximum density for each SDM x scenario combination (i.e. each violin).
    typeNames <- c("avg","sd")
    widthsViolins <- array(dim=c(numExperiments, numPlotSDMs, maxNumScenarios, length(typeNames)), 
                           dimnames=list(1:numExperiments, plotSDMs, 1:maxNumScenarios, typeNames))
    for ( type in typeNames ) { 
      for ( ex in 1:numExperiments ) {
        # Number of scenarios in this experiment (could be different?)
        namesScenarios <- allExperimentsInfo[[ex]]$namesScenarios
        numScenarios <- length(namesScenarios)
        
        # Which type of values "avg" or "sd"        
        if ( type =="avg" ) {
          theseValues <- allExperimentsInfo[[ex]]$thisStatAvg
        } else if ( type == "sd" ) {
          theseValues <- allExperimentsInfo[[ex]]$thisStatSD
        } else {
          stop("Unrecognised type of statistic.  Look at code!")
        }
        
        for ( sdm in plotSDMs ) {
          for ( sc in 1:numScenarios ) {
            # For this experiment, sdm and scenario, the species statistics are:
            thisStat <- as.vector(theseValues[sc, ,sdm,whichCoeffs])
            
            # Check for Inf values
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) {
              thisStat[indInf] <- substituteInfValue
            }
            
            # Check for NA, NaN, or Inf values
            indNA <- which(is.na(thisStat) || is.nan(thisStat))
            if ( length(indNA) > 0 ) {
              warning(paste0(length(indNA), " NA or NaN values in ", type,
                             " for experiment = ", ex,
                             " for SDM = ", sdm, 
                             " and scenario = ", namesScenarios[sc]))
              thisStat <- thisStat[-indNA]
            }
            
            # What is the maximum of the density for each scenario x sdm.
            if ( length(thisStat) < 2 ) {
              message(paste0("Not enough values to form a density in stat '", type,
                             "', experiment = ", ex,
                             ", SDM = ", sdm, 
                             " and scenario = ", namesScenarios[sc]))
            }
            widthsViolins[ex,sdm,sc,type] <- max(density(thisStat, na.rm=TRUE)$y)
          }
        }
      }
      
      # Scale maximum densities so that violin widths are between [0,1]
      maxMaxDensity <- max(widthsViolins[ , , ,type], na.rm=TRUE)
      widthsViolins[ , , ,type] <- widthsViolins[ , , ,type]/maxMaxDensity
    }    
  }
  
  # Cycle through information and add a mean and/or SD plot per SDM x experiment combination.
  for ( ex in 1:numExperiments ) {
    #if ( ex == 2) message("experiment: ", ex)
    # Number and names of scenarios for this experiment.
    numScenarios <- allExperimentsInfo[[ex]]$numScenarios
    namesScenarios <- allExperimentsInfo[[ex]]$namesScenarios
    
    # Were x-axis labels provided (i.e. pretty names for scenarios)?
    if ( is.null(xAxisLabels) ) {
      xAxisLabelsUse <- namesScenarios
    } else if ( length(xAxisLabels) == numScenarios && 
                length(unique(xAxisLabels)) == numScenarios ) {
      # use provided labels.
      xAxisLabelsUse <- xAxisLabels
    } else {
      stop("Invalid x-axis labels given in argument 'xAxisLabels'.")
    }
    
    # Loop to plot accuracy and/or precision of each SDM for this experiment.    
    isFirstRow <- TRUE
    for ( sdm in plotSDMs) {
      #if ( ex == 2 ) message("    SDM: ", sdm)
      plotMedians <- vector(mode="double", length=numScenarios)
      if ( doAccuracy ) {
        # Plot accuracy statistics (a.k.a mean) for this SDM and this experiment.
        # Range of plot
        if ( sdm == "PO" && ! is.null(diffPORange) ) {
          yaxisRange <- diffPORange
        } else {
          yaxisRange <- ylimVals[1:2]
        }
        plot(c(0.5,numScenarios+0.5), yaxisRange, type="n", xaxs="i", xaxt="n", xlab="", 
             yaxt="n", ylab="")
        axis(1, tick=TRUE, labels=FALSE)
        abline(h=horizontalLines[1], lty="dotted")

        # Column headings ...
        if ( isFirstRow  ) {
          isFirstRow <- FALSE
          if  ( ! is.null(columnHeadings) ) {
            if ( length(columnHeadings) == numExperiments )  {
              midxRange <- (numScenarios/2.0) + 0.5
              axis(3, at=midxRange, labels=columnHeadings[ex], tick=FALSE, padj=1.0)
            } else {
              stop("The number of columnHeadings needs to be the same as the number of experiments.")
            }
          }
        } else {
          # Do nothing, no column headings needed when not first row of plot matrix.
        }
        
        # First plot of the row stuff ...
        if ( ex == 1) {
          # Y-axis tick marks, numbers and label.
          if ( all(yaxisRange == c(0.0,1.0))) {
            # Specifically for publication, really!
            yaxisTicks <- axTicks(2, axp=c(yaxisRange[1],yaxisRange[2],4))
            yaxisLabels <- vector(mode = "character", length = 5)
            yaxisLabels[c(1,3,5)] <- format(yaxisTicks[c(1,3,5)])
            yaxisLabels[c(2,4)] <- ""
            yaxisLabels[5] <- paste0("\u2265",yaxisLabels[5])
            axis(2, labels=yaxisLabels, at=yaxisTicks, las=1)
          } else {
            axis(2, labels=TRUE, las=1)
          }
          midyRange <- ((yaxisRange[2] - yaxisRange[1])/2.0) + yaxisRange[1]
          axis(2, at=midyRange, labels="accuracy", tick=FALSE, padj=-3.0)
          
          # SDM title. (NB: doAccuracy == TRUE already to get here!)
          if ( doPrecision ) {
            # Leave until plotting precision.
          } else {
            axis(2, at=midyRange, labels=paste0(sdm, " SDM"), tick=FALSE, padj=-4.5)
          }
        }
        
        # Plot accuracy of species for each scenario.
        for (sc in 1:numScenarios ) {
          # Statistic values for this scenario (will be plotted as a vertical bar).
          thisStat <- as.vector(allExperimentsInfo[[ex]]$thisStatAvg[sc, ,sdm,whichCoeffs])
          
          # Check for NA values.
          indNA <- c(which(is.na(thisStat)), which(is.nan(thisStat)))
          if ( length(indNA) > 0 ) {
            warning(paste0(length(indNA), " NA or NaN values in ", type,
                           " for experiment = ", ex,
                           " for SDM = ", sdm, 
                           " and scenario = ", namesScenarios[sc]))
            thisStat <- thisStat[-indNA]
          }
          
          # Plot ...
          if ( sdm == "PO" && ! is.null(diffPORange)) {
            violinColour <- violinColPO
            pointsColour <- pointsColPO
          } else {
            violinColour <- violinCol
            pointsColour <- pointsCol
          }
          if ( vioPlots ) {
            # Reset infinite values.
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) thisStat[indInf] <- substituteInfValue
            vioplot(thisStat, col=violinColour, border=violinColour, add=TRUE,
                    drawRect = FALSE, at=sc, wex=widthsViolins[ex,sdm,sc,"avg"])
          }
          stripchart(thisStat, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsColour)
          
          # Median value for this scenario.
          plotMedians[sc] <- median(thisStat, na.rm=TRUE)
        }    
        
        # Add median value for all scenarios in this plot.
        if ( medianLines ) {
          points(1:numScenarios, plotMedians, pch="_", col=medianCol)
          lines(1:numScenarios, plotMedians, lty="dotted", col=medianCol)
        }
      }

      if ( doPrecision ) {      
        # Plot precision statistics (a.k.a sd) for this SDM and this experiment.
        # Range of plot
        yaxisRange <- ylimVals[3:4]
        plot(c(0.5,numScenarios+0.5), yaxisRange, type="n", xaxs="i", xaxt="n", xlab="", 
             yaxt="n", ylab="")
        axis(1, tick=TRUE, labels=FALSE)
        abline(h=horizontalLines[2], lty="dotted")
        
        # Column headings ...
        if ( isFirstRow  ) {
          isFirstRow <- FALSE
          if  ( ! is.null(columnHeadings) ) {
            if ( length(columnHeadings) == numExperiments )  {
              midxRange <- (numScenarios/2.0) + 0.5
              axis(3, at=midxRange, labels=columnHeadings[ex], tick=FALSE, padj=1.0)
            } else {
              stop("The number of columnHeadings needs to be the same as the number of experiments.")
            }
          }
        } else {
          # Do nothing, no column headings needed when not first row of plot matrix.
        }
        
        # First plot of the row stuff ...
        if ( ex == 1) {
          # Y-axis tick marks, numbers and label.
          if ( all(yaxisRange == c(0.0,1.0))) {
            # Specifically for publication, really!
            yaxisTicks <- axTicks(2, axp=c(yaxisRange[1],yaxisRange[2],4))
            yaxisLabels <- vector(mode = "character", length = 5)
            yaxisLabels[c(1,3,5)] <- format(yaxisTicks[c(1,3,5)])
            yaxisLabels[c(2,4)] <- ""
            yaxisLabels[5] <- paste0("\u2265",yaxisLabels[5])
            axis(2, labels=yaxisLabels, at=yaxisTicks, las=1)
          } else {
            axis(2, labels=TRUE, las=1)
          }
          midyRange <- ((yaxisRange[2] - yaxisRange[1])/2.0) + yaxisRange[1]
          axis(2, at=midyRange, labels="precision", tick=FALSE, padj=-3.0)

          # SDM title. (NB: doPrecision == TRUE already to get here!)
          if ( doAccuracy ) {
            axis(2, at=yaxisRange[2], labels=paste0("       ", sdm, " SDM"), tick=FALSE, padj=-4.5)
          } else {
            axis(2, at=midyRange, labels=paste0(sdm, " SDM"), tick=FALSE, padj=-4.5)
          }
        }
        
        # Plot precision of species for each scenario.  
        for (sc in 1:numScenarios ) {
          # Statistic values for this scenario (will be plotted as a vertical bar).
          thisStat <- as.vector(allExperimentsInfo[[ex]]$thisStatSD[sc, ,sdm,whichCoeffs])
  
          # Check for NA values.
          indNA <- c(which(is.na(thisStat)), which(is.nan(thisStat)))
          if ( length(indNA) > 0 ) {
            warning(paste0(length(indNA), " NA or NaN values in standard deviation of ",
                           " stat = ", whichStat,
                           " for experiment = ", ex,
                           " for SDM = ", sdm, 
                           " and scenario = ", namesScenarios[sc]))
            thisStat <- thisStat[-indNA]
          }
          
          # Plot ...
          if ( vioPlots ) {
            # Reset infinite values.
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) thisStat[indInf] <- substituteInfValue
            vioplot(thisStat, col=violinCol, border=violinCol, add=TRUE,
                                drawRect = FALSE, at=sc, wex=widthsViolins[ex,sdm,sc,"sd"])
          }
          stripchart(thisStat, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsCol)
          
          # Median value for this scenario.
          plotMedians[sc] <- median(thisStat, na.rm=TRUE)
        }
        
        # Add median value for all scenarios in this plot.
        if ( medianLines ) {
          points(1:numScenarios, plotMedians, pch="_", col=medianCol)
          lines(1:numScenarios, plotMedians, lty="dotted", col=medianCol)
        }
      } # if doPrecision
 
    }  # sdm's (or plot rows).
    
    # Add x-axis for this column, only on bottom row.
    axis(1, 1:numScenarios, xAxisLabelsUse, las=xAxisLabelStyle)
    if ( ! is.null(xAxisTitle) ) {
      if ( length(xAxisTitle) == numExperiments ) {
        midxRange <- (numScenarios/2.0) + 0.5
        axis(1, at=midxRange, labels=xAxisTitle[ex], tick=FALSE, padj=4.0)
      } else if ( length(xAxisTitle) == 1 ) {
        # Will add below.
      } else {
        # Will warn below.
      }
    }
    
  } # experiments (or plot columns).
  
  # Add a title to the plot.
  if ( ! is.null(plotTitle) ) 
    title(plotTitle, outer=TRUE)
  
  # Add an x-axis title to the plot.
  if ( ! is.null(xAxisTitle) ) {
    if ( length(xAxisTitle) == numExperiments ) {
      # Should have already added them above!      
    } else if ( length(xAxisTitle) == 1 ) {
      title(xlab = xAxisTitle, outer=TRUE, line=min(4,omaPar[1]))
    } else {
      warning("Argument 'xAxisTitle' is an unexpected length.  No titles added.")
    }
  }
  
  # Turn off the plotting device, if it is a file.
  if ( plotToFile ) dev.off()
  
  # Return plotting to original settings.
  par(opar)
}  

#-----------------------------------------------------------------------------------------

plotStatisticsComparisons <- function(whichExperiment, whichStats=c(1,1,2,3), 
                                      whichCoeffs=list(1,-1,NULL,NULL), whichSpecies=NULL,
                                      plotSDMs=c("PA","MsPP","Gear"), ylimVals=c(0,1,0,1), 
                                      vioPlots=TRUE, absMeans=TRUE, plotDevice="RStudioGD", 
                                      plotDir=paste0(whichExperiment,"/Plots"), 
                                      fileName="comparison", plotTitle=NULL, plotUnits="cm", 
                                      plotHeight=20, plotWidth=15.75, plotRes=600, 
                                      xAxisLabels=NULL, xAxisTitle=NULL, xAxisLabelStyle=2,
                                      xAxisReverse=FALSE, 
                                      columnHeadings=paste("statistic",1:length(whichStats)),
                                      diffPOIntRange=NULL, medianLines=TRUE, 
                                      accuracyPrecision = "both") {
  
  # Plot multiple statistics for a single experiment (i.e. plot coeff difference, correlation, 
  # etc., for an experiment).  
  #
  # Arguments ...
  # whichExperiment:  a string vector containing a directory name for an experiment.
  # whichStats:       an integer vector that indicates the statistics to plot.  See 
  #                   'makeScenarioStats' for valid values.
  # whichCoeffs:      a list that has the same number of items as whichStats.  Each item
  #                   contains the coefficients to be plotted for its matching statistic. 
  #                   Use NULL per item to get all coeffs for that item.  Use negative 
  #                   values to remove the specified coefficient.
  # whichSpecies:     which species to plot.
  # plotSDMs:         which SDMs from the experiment to plot (and the order).  
  # ylimVals:         a vector of length 4 that gives a y-axis range for mean and sd, eg. 
  #                   c(min mean, max mean, min sd, max sd).
  # xAxisTitle:       can be a single string (in which case it will be placed in the centre
  #                   of the outer plot margins) or a vector of strings, one for each statistic.
  # xAxisReverse:     reverse the positions of the scenarios along the x-axis for all statistics.
  # columnHeadings:   a vector of strings, one for each statistic, to use as column headings.
  #                   Set to "NULL" if no column headings are required.
  # horizontalLines:  REPLACED BY HARDWIRED numbers based on statistics requested. = 0 unless
  #                   whichstats[i] == 2, then equal to 1 for accuracy only!
  #                   a numStats x 2 matrix that give the y values at which horizontal lines
  #                   will be drawn for each accuracy (mean) and each precision (sd) plot.
  # diffPOIntRange:   a vector of length 2 that gives a different y-axis range for PO SDM 
  #                   accuracy plots for intercept statistics (e.g. c(min mean, max mean)). 
  #                   Set to NULL to use same y-axis range as other SDMs
  # medianLines:      True to include a line to indicate the median value per scenario of
  #                   each plot, false to not include.
  # accurayPrecision: whether to plot accuracy or precision or both (these are rows of plots).
  #                   Acceptable values are "accuracy", "precision", "both".
  
  
  # Get numbers of things to set up plotting layout.
  numPlotStats <- length(whichStats)
  numPlotSDMs <- length(plotSDMs)
  substituteInfValue <- 1.0e+25
  
  # Check accuracy or precision or both.
  if ( accuracyPrecision == "accuracy" ) {
    doAccuracy <- TRUE
    doPrecision <- FALSE
  } else if ( accuracyPrecision == "precision" ) {
    doAccuracy <- FALSE
    doPrecision <- TRUE
  } else if ( accuracyPrecision == "both" ) {
    doAccuracy <- TRUE
    doPrecision <- TRUE
  } else {
    doAccuracy <- FALSE
    doPrecision <- FALSE
    stop("No valid comparisons have been requested, check argument 'accuracyPrecision'.")
  }
  
  # Colours for plotting.
  pointsCol <- "black"
  violinCol <- "grey"
  pointsColPO <- colourBlindRGB("blue")
  violinColPO <- colourBlindRGB("skyBlue")
  medianCol <- colourBlindRGB("orange")

  # Check dimensions are valid.  
  if ( length(xAxisReverse) == 1 ) {
    xAxisReverseUse <- xAxisReverse
  } else {
    stop("Argument 'xAxisReverse' should be a scalar.  Please check using correct function!")
  }
  
  # Check dimensions are valid.  
  if ( length(whichExperiment) > 1 ) {
    stop("Argument 'whichExperiment' should be a scalar. Please check using correct function!")
  }
  
  # Check dimensions are valid.  
  if ( length(whichCoeffs) != numPlotStats || class(whichCoeffs) != "list") {
    stop("Argument 'whichCoeffs' should be a list of the same length as 'whichStats'.")
  }
  
  # Load in this experiment's results
  exDir <- whichExperiment
  retLst <- loadScenariosResults(exDir)
  
  # Check requested statistics are contained within this experiment.
  if ( ! all(whichStats %in% c(1:retLst$statsObj$numStats)) ) {
    stop("Unrecognised statistic requested from experiment '", exDir, "'.")
  }
  
  # Check requested coefficients, per requested statistic, are contained within this experiment.
  namesCoeffs <- dimnames(retLst$statsObj$stats[[1]]$avg)[[4]]
  numCoeffs <- length(namesCoeffs)
  doDiffPOYAxis <- FALSE
  for ( stat in 1:numPlotStats ) {
    # What are the requested coefficients for this statistic?
    theseCoeffs <- whichCoeffs[[stat]]

    # What are all the possible coefficients for this statistics.
    whichStat <- whichStats[stat]
    allCoeffs <- 1:dim(retLst$statsObj$stats[[whichStat]]$avg)[[4]]
    
    # Are they valid?
    if ( is.null(theseCoeffs) ) {
      whichCoeffs[[stat]] <- allCoeffs
    } else if ( ! all( abs(theseCoeffs) %in% allCoeffs) ) {
      # NB: absoulte value as negatives can be used to remove coefficients too.
      stop("Unrecognised coefficient requested for statistic at item ",stat," in 'whichCoeffs'.")
    }
    
    # Do a different axis for PO intercept accuracy?
    if ( ! is.null(diffPOIntRange) && whichStat == 1 && 1 %in% whichCoeffs[[stat]] ) {
      doDiffPOYAxis <- TRUE
    }
  }
  
  # Check requested SDMs are contained within this experiment.
  if ( ! all(plotSDMs %in% retLst$resLstAll[[1]]$namesValidSDM) ) {
    stop("Unrecognised SDM requested from experiment '", exDir, "'.")
  }
  
  # Check requested species are contained within this experiment.
  if ( is.null(whichSpecies) ) {
    # all good, going to take whatever species are within the experiment.
    whichSpecies <- dimnames(retLst$statsObj$stats[[1]]$avg)[[2]]
  } else if ( ! all(whichSpecies %in% dimnames(retLst$statsObj$stats[[1]]$avg)[[2]]) ) {
    stop("Unrecognised species requested from experiment '", exDir, "'.")
  }

  # Do we need absolute values for the accuracy?  
  if ( absMeans ) {
    for ( stat in 1:retLst$statsObj$numStats ) {
      retLst$statsObj$stats[[stat]]$avg <- abs(retLst$statsObj$stats[[stat]]$avg)
    }
  }
  
  # Do we need to reverse the scenario direction?
  numScenarios = retLst$scenariosObj$numScenarios
  namesScenarios = retLst$scenariosObj$prettyNamesScenarios
  if ( xAxisReverseUse ) {
    # Reverse the scenarios on the x-axis for this experiment.
    namesScenarios <- rev(namesScenarios)
    for ( stat in 1:retLst$statsObj$numStats ) {
      retLst$statsObj$stats[[stat]]$avg <- retLst$statsObj$stats[[stat]]$avg[numScenarios:1, , , ,drop=FALSE]
      retLst$statsObj$stats[[stat]]$sd <- retLst$statsObj$stats[[stat]]$sd[numScenarios:1, , , ,drop=FALSE]
    }
  }
  
  # Plot directory for statistics comparison level plots.
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Which device are we printing to?
  if ( plotDevice == "RStudioGD" ) {
    # plot to the R studio plot window.
    plotToFile <- FALSE
  } else {
    fileName <- makeFileName(fileName, plotDir, plotDevice)
    argList <- list(filename = fileName, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
    plotToFile <- TRUE
  }
  
  # Divide plotting into the number of rows and columns required.
  omaPar <- par("mar")         # use default single plot per page margins
  marPar <- c(0.0,0.5,1.3,0.0) # small margin on top to separate rows and allow for column headings.
  if ( is.null(xAxisTitle) ) {
    # No x-axis title, change margins to gain plotting area.
    omaPar[1] <- omaPar[1] - 1.0  # leave space for x-axis ticks and labels.  
  }
  if ( is.null(plotTitle) ) {
    # No main title, change margins to gain plotting area.
    omaPar[3] <- 0.0  # NB: there is already a bit at top due to mar settings.
  }
  if ( doDiffPOYAxis && "PO" %in% plotSDMs && ylimVals[1:2] == c(0.0,1.0)) {
    #  Probably only for publication version where using >= sign in axis tick mark values.
    omaPar[4] <- 3.1
  } else {
    omaPar[4] <- omaPar[4]+0.5
  }
  if ( doAccuracy && doPrecision ) {
    # Both accuracy and precision is required.
    opar <- par(mfcol=c(numPlotSDMs*2,numPlotStats), oma = omaPar, mar = marPar, cex=0.66)
  } else {
    # Either accuracy or precision is required but not both.
    opar <- par(mfcol=c(numPlotSDMs,numPlotStats), oma = omaPar, mar = marPar, cex=0.66)
  }

  # Calculate widths for violin plots for avg (accuracy) and SDs (precision), if necessary.
  if ( vioPlots ) {
    # Get the maximum density for each SDM x scenario combination (i.e. each violin).
    typeNames <- c("avg","sd")
    widthsViolins <- array(dim=c(numPlotStats, numPlotSDMs, numScenarios, length(typeNames)), 
                           dimnames=list(1:numPlotStats, plotSDMs, 1:numScenarios, typeNames))
    for ( type in typeNames ) { 
      for ( stat in 1:numPlotStats ) {
        # What is the requested statistic for this position in whichStats?
        whichStat <- whichStats[stat]
        whichCoeffsForStat <- whichCoeffs[[stat]]

        # Which type of values "avg" or "sd"        
        if ( type =="avg" ) {
          theseValues <- retLst$statsObj$stats[[whichStat]]$avg[ ,whichSpecies, ,whichCoeffsForStat,drop=FALSE]
          
        } else if ( type == "sd" ) {
          theseValues <- retLst$statsObj$stats[[whichStat]]$sd[ ,whichSpecies, ,whichCoeffsForStat,drop=FALSE]
        } else {
          stop("Unrecognised type of statistic.  Code may require update!")
        }
        
        for ( sdm in plotSDMs ) {
          for ( sc in 1:numScenarios ) {
            # For this experiment, sdm and scenario, the species statistics are:
            thisStat <- as.vector(theseValues[sc, ,sdm, ])
            
            # Check for Inf values
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) {
              thisStat[indInf] <- substituteInfValue
            }
            
            # Check for NA, NaN, or Inf values
            indNA <- which(is.na(thisStat) || is.nan(thisStat))
            if ( length(indNA) > 0 ) {
              warning(paste0(length(indNA), " NA or NaN values in ", type,
                             " for statistic = ", whichStat,
                             " for SDM = ", sdm, 
                             " and scenario = ", namesScenarios[sc]))
              thisStat <- thisStat[-indNA]
            }
            
            # What is the maximum of the density for each scenario x sdm.
            if ( length(thisStat) < 2) {
              message(paste0("Less than two values left in statistic to form density in ", type,
                             " for statistic = ", whichStat,
                             " for SDM = ", sdm, 
                             " and scenario = ", namesScenarios[sc]))
            }
            widthsViolins[stat,sdm,sc,type] <- max(density(thisStat, na.rm=TRUE)$y)
          }
        }
        
        # Scale maximum densities so that violin widths are between [0,1]
        maxMaxDensity <- max(widthsViolins[stat , , ,type], na.rm=TRUE)
        widthsViolins[stat , , ,type] <- widthsViolins[stat , , ,type]/maxMaxDensity
      }
    }    
  }
  
  # Were x-axis labels provided (i.e. pretty names for scenarios)?
  if ( is.null(xAxisLabels) ) {
    xAxisLabelsUse <- namesScenarios
  } else if ( length(xAxisLabels) == numScenarios && 
              length(unique(xAxisLabels)) == numScenarios ) {
    # use provided labels.
    xAxisLabelsUse <- xAxisLabels
  } else {
    stop("Invalid x-axis labels given in argument 'xAxisLabels'.")
  }
  
  # Cycle through information and add a mean and SD plot per SDM x statistic combination.
  doForthAxis <- FALSE
  for ( stat in 1:numPlotStats ) {
    whichStat <- whichStats[stat]
    whichCoeffsForStat <- whichCoeffs[[stat]]
    
    if ( whichStat == 2 ) {
      # Assumes that this is the correlation statistic!
      horizontalLine <- 1.0
    } else {
      horizontalLine <- 0.0
    }
    
    # Loop to plot accuracy and precision of each SDM for this statistic.    
    isFirstRow <- TRUE
    for ( sdm in plotSDMs) {
      plotMedians <- vector(mode="double", length=numScenarios)
      if ( doAccuracy ) {
        # Plot accuracy statistics (a.k.a 'avg') for this SDM and this statistic.
        # Range of plot
        if ( sdm == "PO" && ! is.null(diffPOIntRange) && whichStat == 1 && 1 %in% whichCoeffsForStat) {
          yaxisRange <- diffPOIntRange
        } else {
          yaxisRange <- ylimVals[1:2]
        }
        plot(c(0.5,numScenarios+0.5), yaxisRange, type="n", xaxs="i", xaxt="n", xlab="", 
             yaxt="n", ylab="")
        axis(1, tick=TRUE, labels=FALSE)
        abline(h=horizontalLine, lty="dotted")
      
        # First plot of the row stuff ...
        if ( stat == 1) {
          # Y-axis tick marks, numbers and label.
          if ( sdm == "PO" && doDiffPOYAxis ) {
            axis(2, labels=TRUE, col=pointsColPO, col.axis=pointsColPO, las=1)
            doForthAxis <-TRUE
          } else {
            if ( all(yaxisRange == c(0.0,1.0))) {
              # Specifically only for publication, really!
              yaxisTicks <- axTicks(2, axp=c(yaxisRange[1],yaxisRange[2],4))
              yaxisLabels <- vector(mode = "character", length = 5)
              yaxisLabels[c(1,3,5)] <- format(yaxisTicks[c(1,3,5)])
              yaxisLabels[c(2,4)] <- ""
              yaxisLabels[5] <- paste0("\u2265",yaxisLabels[5])   
              axis(2, labels=yaxisLabels, at=yaxisTicks, las=1)
            } else {
              axis(2, labels=TRUE, las=1)
            }
          }
          midyRange <- ((yaxisRange[2] - yaxisRange[1])/2.0) + yaxisRange[1]
          axis(2, at=midyRange, labels="accuracy", tick=FALSE, padj=-3.0)
          # SDM title. (NB: doAccuracy == TRUE already to get here!)
          if ( doPrecision ) {
            # Leave until plotting precision.
          } else {
            axis(2, at=midyRange, labels=paste0(sdm, " SDM"), tick=FALSE, padj=-4.5)
          }
        }
      
        # Last plot of row, if necessary.
        if ( stat == numPlotStats && doForthAxis && sdm == "PO") {
          if ( all(yaxisRange == c(0.0,1.0))) {
            # Specifically only for publication, really!
            yaxisTicks <- axTicks(2, axp=c(yaxisRange[1],yaxisRange[2],4))
            yaxisLabels <- vector(mode = "character", length = 5)
            yaxisLabels[c(1,3,5)] <- format(yaxisTicks[c(1,3,5)])
            yaxisLabels[c(2,4)] <- ""
            yaxisLabels[5] <- paste0("\u2265",yaxisLabels[5])   
            axis(4, labels=yaxisLabels, at=yaxisTicks, las=1)
          } else {
            axis(4, labels=TRUE, las=1)
          }
          doForthAxis <- FALSE
        }
        
        # Column headings ...
        if ( isFirstRow  ) {
          isFirstRow <- FALSE
          if  ( ! is.null(columnHeadings) ) {
            if ( length(columnHeadings) == numPlotStats )  {
              midxRange <- (numScenarios/2.0) + 0.5
              axis(3, at=midxRange, labels=columnHeadings[stat], tick=FALSE, padj=1.0)
            } else {
              stop("The length of 'columnHeadings' needs to be the same as the number of requested statistics.")
            }
          }
        } else {
          # Do nothing, no column headings needed when not first row of plot matrix.
        }
      
        # Plot accuracy of species for each scenario.
        for (sc in 1:numScenarios ) {
          # Statistic values for this scenario (will be plotted as a vertical bar).
          thisStat <- as.vector(retLst$statsObj$stats[[whichStat]]$avg[sc,whichSpecies,sdm,whichCoeffsForStat])
  
          # Check for NA values.
          indNA <- c(which(is.na(thisStat)), which(is.nan(thisStat)))
          if ( length(indNA) > 0 ) {
            warning(paste0(length(indNA), " NA or NaN values in ", type,
                           " for statistic = ", whichStat,
                           " for SDM = ", sdm, 
                           " and scenario = ", namesScenarios[sc]))
            thisStat <- thisStat[-indNA]
          }
          
          # Plot ...
          if ( sdm == "PO" && ! is.null(diffPOIntRange) && whichStat == 1 && 1 %in% whichCoeffsForStat) {
            violinColour <- violinColPO
            pointsColour <- pointsColPO
          } else {
            violinColour <- violinCol
            pointsColour <- pointsCol
          }
          if ( vioPlots ) {
            # Reset infinite values.
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) thisStat[indInf] <- substituteInfValue
            vioplot(thisStat, col=violinColour, border=violinColour, add=TRUE,
                    drawRect = FALSE, at=sc, wex=widthsViolins[stat,sdm,sc,"avg"])
          }
          stripchart(thisStat, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsColour)
          
          # Median value for this scenario.
          plotMedians[sc] <- median(thisStat, na.rm=TRUE)
        }    
        
        # Add median value for all scenarios in this plot.
        if ( medianLines ) {
          points(1:numScenarios, plotMedians, pch="_", col=medianCol)
          lines(1:numScenarios, plotMedians, lty="dotted", col=medianCol)
        }
      } # if ( doAccuracy )       
      
      if ( doPrecision ) {
        # Plot precision statistics (a.k.a sd) for this SDM and this statistic.
        # Range of plot
        yaxisRange <- ylimVals[3:4]
        plot(c(0.5,numScenarios+0.5), yaxisRange, type="n", xaxs="i", xaxt="n", xlab="", 
             yaxt="n", ylab="")
        axis(1, tick=TRUE, labels=FALSE)
        abline(h=0.0, lty="dotted")         ## Hardwired.
        
        # Column headings ...
        if ( isFirstRow  ) {
          isFirstRow <- FALSE
          if  ( ! is.null(columnHeadings) ) {
            if ( length(columnHeadings) == numPlotStats )  {
              midxRange <- (numScenarios/2.0) + 0.5
              axis(3, at=midxRange, labels=columnHeadings[stat], tick=FALSE, padj=1.0)
            } else {
              stop("The length of 'columnHeadings' needs to be the same as the number of requested statistics.")
            }
          }
        } else {
          # Do nothing, no column headings needed when not first row of plot matrix.
        }
        
        # First of row stuff ...
        if ( stat == 1) {
          # Y-axis tick marks, numbers and label.
          if ( all(yaxisRange == c(0.0,1.0))) {
            # Specifically for publication, really!
            yaxisTicks <- axTicks(2, axp=c(yaxisRange[1],yaxisRange[2],4))
            yaxisLabels <- vector(mode = "character", length = 5)
            yaxisLabels[c(1,3,5)] <- format(yaxisTicks[c(1,3,5)])
            yaxisLabels[c(2,4)] <- ""
            yaxisLabels[5] <- paste0("\u2265",yaxisLabels[5])
            axis(2, labels=yaxisLabels, at=yaxisTicks, las=1)
          } else {
            axis(2, labels=TRUE, las=1)
          }
          midyRange <- ((yaxisRange[2] - yaxisRange[1])/2.0) + yaxisRange[1]
          axis(2, at=midyRange, labels="precision", tick=FALSE, padj=-3.0)
          
          # SDM title.
          if ( doAccuracy ) {
            axis(2, at=yaxisRange[2], labels=paste0("       ", sdm, " SDM"), tick=FALSE, padj=-4.5)
          } else {
            axis(2, at=midyRange, labels=paste0(sdm, " SDM"), tick=FALSE, padj=-4.5)
          }
        }
  
            
        # Plot precision of species for each scenario.  
        for (sc in 1:numScenarios ) {
          # Statistic values for this scenario (will be plotted as a vertical bar).
          thisStat <- as.vector(retLst$statsObj$stats[[whichStat]]$sd[sc,whichSpecies,sdm,whichCoeffsForStat])
          
          # Check for NA values.
          indNA <- c(which(is.na(thisStat)), which(is.nan(thisStat)))
          if ( length(indNA) > 0 ) {
            warning(paste0(length(indNA), " NA or NaN values in standard deviation of ",
                           " for statistic = ", whichStat,
                           " for SDM = ", sdm, 
                           " and scenario = ", namesScenarios[sc]))
            thisStat <- thisStat[-indNA]
          }
          
          # Plot ...
          if ( vioPlots ) {
            # Reset infinite values.
            indInf <- which(is.infinite(thisStat))
            if ( length(indInf) > 0 ) thisStat[indInf] <- substituteInfValue
            vioplot(thisStat, col=violinCol, border=violinCol, add=TRUE,
                    drawRect = FALSE, at=sc, wex=widthsViolins[stat,sdm,sc,"sd"])
          }
          stripchart(thisStat, at=sc, vertical=TRUE, pch="_", add=TRUE, col=pointsCol)
        
          # Median value for this scenario.
          plotMedians[sc] <- median(thisStat, na.rm=TRUE)
        }    
        
        # Add median value for all scenarios in this plot.
        if ( medianLines ) {
          points(1:numScenarios, plotMedians, pch="_", col=medianCol)
          lines(1:numScenarios, plotMedians, lty="dotted", col=medianCol)
        }
      } # if ( doPrecision )
    }  # sdm's (or plot rows).
    
    # Add x-axis for this column only on bottom row.
    axis(1, 1:numScenarios, xAxisLabelsUse, las=xAxisLabelStyle)
    if ( ! is.null(xAxisTitle) ) {
      if ( length(xAxisTitle) == numPlotStats ) {
        midxRange <- (numScenarios/2.0) + 0.5
        axis(1, at=midxRange, labels=xAxisTitle[stat], tick=FALSE, padj=4.0)
      } else if ( length(xAxisTitle) == 1 ) {
        # Will add below.
      } else {
        # Will warn below.
      }
    }
    
  } # statistics (or plot columns).
  
  # Add a title to the plot.
  if ( ! is.null(plotTitle) ) 
    title(plotTitle, outer=TRUE)
  
  # Add an x-axis title to the plot.
  if ( ! is.null(xAxisTitle) ) {
    if ( length(xAxisTitle) == numPlotStats ) {
      # Should have already added them above!      
    } else if ( length(xAxisTitle) == 1 ) {
      title(xlab = xAxisTitle, outer=TRUE, line=min(4,omaPar[1]))
    } else {
      warning("Argument 'xAxisTitle' is an unexpected length.  No x-axis titles added.")
    }
  }
  
  # Turn off the plotting device, if it is a file.
  if ( plotToFile ) dev.off()
  
  # Return plotting to original settings.
  par(opar)
}  

#-----------------------------------------------------------------------------------------

makeDensityPlots <- function(envirVals, scenariosObj, 
                             whichScenarios=1:scenariosObj$numScenarios, whichRuns=c(1,5000), 
                             legendTitle=NULL, xLabel="covariate value") {
  
  # Plot the density function for all the environment values (within the domain) and then
  # each successive scenario's density function for the environment values at the survey
  # sample sites (specified as cells in the PA data for each scenario)
  
  # envirVals:      a vector of the selected environment variable's values for cells within
  #                 the domain.  It is assumed that envirVals and the results from the 
  #                 scenariosObj runs are based on the same domain (i.e. same cells).
  # scenariosObj:   a scenarios object that contains information about the scenarios run in 
  #                 an experiment.
  # whichScenarios: which scenarios to use in plotting
  # whichRuns:      which simulation run to use the survey sample locations from.  Usually
  #                 only have the option of 1 or 5000.
  
  # Check scenarios requested are valid.
  if ( all(whichScenarios %in% 1:scenariosObj$numScenarios) ) {
    numReqScenarios <- length(whichScenarios)
  } else {
    stop("Invalid scenarios requested.")
  }
  
  # Going to need the range of the y-axis in the plot to be based on all scenarios densities.
  # So, will need to load and store the sample sites for each scenario.
  experimentDir <- scenariosObj$scenariosDir
  numRuns <- length(whichRuns)
  densities <- vector(mode="list",length = numReqScenarios*numRuns)
  yRange <- NULL
  for ( sc in 1:numReqScenarios ) {
    # Get the directory for each scenario in this experiment.
    thisScenarioNum <- whichScenarios[sc]
    scenarioDir <- paste0(experimentDir,"/",scenariosObj$namesScenarios[thisScenarioNum],"/")
    
    for ( run in 1:numRuns ) {
      # Make the filename for requested run's survey sample sites.
      thisRunNum <- whichRuns[run]
      dataFile <- paste0(scenarioDir, "DataDumps/DataExample-run", thisRunNum,".RData")
      
      # Load the requested data (i.e. the survey object)
      load(dataFile)
    
      # Save the density for each of the scenarios.
      densityRes <- density(envirVals[surveyObj$rowsInCells])
      indList <- (sc-1)*numRuns + run
      densities[[indList]] <- list(x = densityRes$x, y=densityRes$y)
      
      # While here, get the y-axis range required by the scenario densities.
      yRange <- range(c(densityRes$y, yRange))
    }
  }
  
  # Get the density function for all the environment values.
  densityEnv <- density(envirVals)
  
  # Get the required range for all densities.
  xRange <- range(densityEnv$x)
  yRange <- range(c(densityEnv$y, yRange))
  
  # Plot densities.
  plot(densityEnv, ylim=yRange, main="", xlab=xLabel, ylab="density")
  abline(h = 0.0, lty="dotted")
  ltyLines <- c("dotdash","dashed")
  scenarioPal <- brewer.pal(numReqScenarios,"Dark2")
  for ( sc in 1:numReqScenarios ) {
    for ( run in 1:numRuns ) {
      indList <- (sc-1)*numRuns + run
      lines(densities[[indList]]$x, densities[[indList]]$y, col=scenarioPal[sc], lty=ltyLines[run])
    }
  }
  
  # Legend for scenarios
  legend("topleft", legend=c("envir", scenariosObj$prettyNamesScenarios[whichScenarios]), 
         col=c("#000000", scenarioPal), lty=c("solid",rep("solid",numReqScenarios)), 
         title=legendTitle)

  # Legend for runs
  if ( numRuns > 1 )
    legend("topright", legend=paste0("run ",whichRuns), lty=ltyLines, title="simulation")
}

#-----------------------------------------------------------------------------------------

plotScenarioMaps <- function(whichExperiment, whichSpecies=NULL, whichSDM=NULL, 
                             whichScenarios=NULL, plotTrue=TRUE, plotDevice="RStudioGD", 
                             plotDir=getwd(), fileName="mapCompare", plotUnits="cm", 
                             plotHeight=20, plotWidth=15.75, plotRes=600, plotNames=NULL,
                             numPlotCols=1, coastLine=NULL, colLand=NULL, nameLand=NULL, 
                             posNameLand=NULL, useAlpha=TRUE, 
                             plotDiff=FALSE, useValues = "numInd", colNA=NA, ...){
  
  
  # Plot scenario maps from a particular experiment for the given species.  Can also plot
  # the true map for this experiment and species, in which case, it occupies the first plot.
  #
  # Arguments ...
  # whichExperiment: a string vector containing a directory name of an experiment.
  # whichSpecies:    which species to plot.  
  # whichSDM:        which SDM's estimate from the scenario to plot.  If a scalar, all 
  #                  scenarios will use estiamtes from this SDM.  If a vector the same 
  #                  length as whichScenarios, the pair of SDM x scenario will be used.
  # whichScenarios:  a vector of scenario directory names from which to source results.
  # plotTrue:        if true, plots the true values for the given species in addition to 
  #                  the given scenarios.
  # plotNames:       a vector of character strings that will be plotted within each plot
  #                  to distinguish which plot is which. Include a name for the true plot,
  #                  in the first position, if plotTrue = TRUE.
  # numPlotCols:     number of columns of plots to put on the page (will figure out the 
  #                  number of rows automatically).
  # coastLine:       a spatial polygon giving the coast line for the map 
  # colLand:         colour of the land defined by coastLine, default is transparent.
  # nameLand:        name with which to label land mass, default is no label.
  # posNameLand:     position at which to centre the label of land mass.
  # useAlpha:        if true then will use the true and mean estimates for this coefficient,
  #                  otherwise, will set all alphas to zero to remove constant difference 
  #                  in scale.
  # plotDiff:        if true then plot the difference of the estimated cell value from the 
  #                  true cell value (in this case plot of true values will not be plotted).
  # useValues:       Use these values to form the map:
  #                    "intensity" - plot intensity of cell (lambda_ck)
  #                    "numInd"    - plot number of expected individuals per cell 
  #                                  (lambda_ck * areaCell)
  # colNA:           colour of the background (or masked cells), default is transparent

  # Error check.
  if ( is.null(whichExperiment) || length(whichExperiment) > 1 )
    stop("Problem with experiment directory name, please check.")
  if ( is.null(whichSpecies) || length(whichSpecies) > 1 ) 
    stop("Problem with specified species.  Need one and only one!")
  if ( is.null(whichSDM) ) stop("No SDM given.  At least one is necessary.")
  if ( is.null(whichScenarios) ) stop("No scenarios given.")

  # Numbers of things.
  numScenarios <- length(whichScenarios)
  if ( length(whichSDM) == 1 && numScenarios > 1) whichSDM <- rep(whichSDM, numScenarios)
  numPlots <- numScenarios 
  if ( plotTrue && ! plotDiff) numPlots <- numScenarios + 1
  numPlotRows <- ceiling(numPlots/numPlotCols)
  
  # More error checking.
  if ( ! is.null(plotNames) ) {
    if ( length(plotNames) != numPlots )
      stop("Problem with number of plot names provided, please check.")
  }
  if ( length(whichSDM) != numScenarios ) 
    stop("One SDM only to be specified for each scenario (assuming not same for all).")

  
  # Plot directory 
  if ( ! dir.exists(plotDir) ) {
    dir.create(plotDir, recursive = TRUE)
  }
  
  # Which device are we printing to?
  if ( plotDevice == "RStudioGD" ) {
    # plot to the R studio plot window.
    plotToFile <- FALSE
  } else {
    fileName <- makeFileName(fileName, plotDir, plotDevice)
    argList <- list(filename = fileName, width = plotWidth, height = plotHeight,
                    units = plotUnits, res = plotRes)
    do.call(plotDevice, argList)
    plotToFile <- TRUE
  }
  
  # Divide plotting into the number of rows and columns required.
  omaPar <- par("mar")         
  omaPar[1] <- 2.4             # Space at bottom to fit x-axis tick mark labels.
  omaPar[2] <- 1.0
  omaPar[3] <- 0.5             # Small space to fit top of y-axis tick mark labels.
  omaPar[4] <- 0.0
  
  marPar <- c(1.0,4.0,0.0,1.0) # small margin on bottom to separate rows and allow for tick marks.
                               # margin on right side to allow for raster legend.
  opar <- par(mfrow=c(numPlotRows,numPlotCols), oma = omaPar, mar = marPar, cex=0.66)
  message("oma = ",paste0(omaPar,collapse=", "))
  message("mar = ",paste0(marPar,collapse=", "))
  iPlotRow <- 0
  
  # Load cell object data and domain.
  savedDir <- paste0(whichExperiment, "/", whichScenarios[1], "/DataDumps")
  if ( ! dir.exists(savedDir) ) 
    stop(paste0("Unable to load data for map as directory doesn't exist.\n", 
                "  Directory tried: '", savedDir,"'"))
  tmp <- load(paste0(savedDir,"/DataDump-AllRuns.RData"))
  rlPlot <- domainObj$mask
  
  # If plot names are to be provided, get position.  HARDWIRED to be top right of plot!
  xyName <- as.vector(simObj$ext)[c(2,4)] * 0.90
  
  # Do we want to compare intensities or expected numbers of individuals?
  if ( useValues == "intensity" ) {
    areaCell <- 1
    legendTextTrue <- "True (per km^2)"
    if ( plotDiff ) {
      legendText <- "Difference (per km^2)"
    } else {
      legendText <- "Estimated (per km^2)"
    }
      
  } else if ( useValues == "numInd" ){
    areaCell <- cellObj$areaCell
    legendTextTrue <- "True (per cell)"
    if ( plotDiff ) {
      legendText <- "Difference (per cell)"
    } else {
        legendText <- "Estimated (per cell)"
    }

  } else {
    stop("Unrecognised value for argument 'useValues'.  Please check!")
  }
  
  # Plot raster of true species distribution values.  
  trueCoeffs <- simObj$beta[ ,whichSpecies]
  if ( ! useAlpha ) trueCoeffs[1] <- 0
  trueValues <- as.matrix(lambda.cell(simObj$lambda.formula, trueCoeffs, cellObj$covars)) * areaCell

  # Plot map of true values, if required.
  if ( plotTrue ) {
    if ( ! plotDiff ) {
      rlPlot[cellObj$cells] <- trueValues 
      plot(rlPlot, asp=1, xaxt="n", ylab="Northing", legend=FALSE, colNA=colNA, ...)
      axis(1, tick=TRUE, labels=FALSE)  # Will always be plot at top of column.
      if ( ! is.null(plotNames) ) text(xyName[1], xyName[2], labels=plotNames[1], adj=c(1,1)) # Hardwired adjustement
      if ( ! is.null(coastLine) ) plot(coastLine, add=TRUE, col=colLand)
      if ( ! is.null(posNameLand) && ! is.null(nameLand) ) 
        text(posNameLand[1], posNameLand[2], nameLand, cex=0.8)
      plot(rlPlot, asp=1, xaxt="n", legend.only=TRUE, 
           legend.args=list(text=legendTextTrue, side=2,cex=0.50), ...)
      iPlotRow <- iPlotRow + 1
    }
  }

  # Plot scenario maps.
  for ( sc in 1:numScenarios ) {
    # Load saved data relevant to given scenario for this experiment.
    savedDir <- paste0(whichExperiment, "/", whichScenarios[sc], "/Results")
    if ( ! dir.exists(savedDir) ) 
      stop(paste0("Unable to load data for scenario map as directory doesn't exist.\n", 
                  "  Directory tried: '", savedDir,"'"))
    tmp <- load(paste0(savedDir,"/Results.RData"))
    
    # Get the estimated cell values for this scenario.
    indSDM <- which(resLst$namesValidSDM == whichSDM[sc])
    internalNameSDM <- resLst$validSDM[indSDM]
    estCoeffs <- resLst[[internalNameSDM]]$coeffs$beta[ ,whichSpecies, ]
    meanEstCoeffs <- apply(estCoeffs,1,mean, na.rm=TRUE)
    if ( ! useAlpha ) {
      meanEstAlpha <- meanEstCoeffs[1]
      meanEstCoeffs[1] <- 0  
    }
    cellValues <- as.matrix(lambda.cell(simObj$lambda.formula, meanEstCoeffs, cellObj$covars)) * areaCell

    # Plot map of this scenario's first estimate for the given species.
    plotName <- NULL
    if ( plotDiff ) {    
      rlPlot[cellObj$cells] <- trueValues - cellValues
      if ( ! is.null(plotNames) ) plotName <- plotNames[sc]
    } else {
      rlPlot[cellObj$cells] <- cellValues
      if ( ! is.null(plotNames) ) plotName <- plotNames[as.integer(plotTrue) + sc]
    }
    plot(rlPlot, asp=1, xaxt="n", ylab="Northing", legend=FALSE, colNA=colNA, ...)
    iPlotRow <- iPlotRow + 1
    if ( iPlotRow > numPlotRows ) iPlotRow <- 1
    if ( iPlotRow < numPlotRows ) {
      # Not bottom plot of column 
      axis(1, tick=TRUE, labels=FALSE)
    } else {
      # Bottom plot, do x-axis
      axis(1, tick=TRUE, labels=TRUE)
      xRange <- range(cellObj$xy$x)
      midxRange <- ((xRange[2] - xRange[1])/2.0) + xRange[1]
      axis(1, at=midxRange, labels="Easting", tick=FALSE, padj=2.0) # x-axis labels
    }
    if ( ! is.null(plotName) ) text(xyName[1], xyName[2], labels=plotName, adj=c(1,1)) # Hardwired adjustement
    if ( ! is.null(coastLine) ) plot(coastLine, add=TRUE, col=colLand)
    if ( ! is.null(posNameLand) && ! is.null(nameLand) ) 
      text(posNameLand[1], posNameLand[2], nameLand, cex=0.8)
    legend.text <- "Difference (per km^2)"
    plot(rlPlot, asp=1, xaxt="n", ylab="Northing", legend.only=TRUE,
         legend.args=list(text=legendText, side=2, cex=0.50), ...)
  }
  
  # Turn off the plotting device, if it is a file.
  if ( plotToFile ) dev.off()
  
  # Return plotting to original settings.
  par(opar)
  
  # Return last cellValues.
  return(rlPlot)
}  

#-----------------------------------------------------------------------------------------

makeScenarioZeta <- function(zeta, zetaMultiplier, returnZeta=TRUE){
  
  # Scale the given zeta based on the zetaMultiplier value (centre using mean first).
  # Returns the scaled zeta with the same type, dimensions and names as zeta.
  # Can return probability of detection = exp(zeta) if returnZeta is FALSE.
  #
  # Arguments ...
  # zeta:           a matrix (numGears x numSpecies) where probDet = exp(zeta).
  # zetaMultiplier: a scalar that is used to increase or decrease the "distance" of the
  #                 gear coefficients from the mean, per species.
  
  # Get number of gears.
  dimZeta <- dim(zeta)
  if ( is.null(dimZeta) ) {
    numGears <- length(zeta)
  } else {
    numGears <- dimZeta[1]
  }
  
  # Centre the values at zero using the mean for each species.
  muSpecies <- apply(zeta, 2, mean)
  muSpecies.mat <- zeta                  # Initialising with same size and names.
  muSpecies.mat[] <- rep(muSpecies, each=numGears)
  zeta.cent <- zeta - muSpecies.mat
  
  # Apply the scaling 
  zeta.mult <- zeta.cent * zetaMultiplier
  
  # Un-centre the new values so that the centre is the same as the original.
  zeta.new <- zeta.mult + muSpecies.mat
  
  # Return value.
  if ( returnZeta ) {
    return(zeta.new)
  } else {
    return(exp(zeta.new))
  }
  
}  

#-----------------------------------------------------------------------------------------

runExperiment <- function(scenariosDir, numRuns = 1000,
                          numPA = 3500, numClusters = 0, 
                          gammaMultiplier = 1.0, deltaMultiplier = 0.0,
                          zetaMultiplier = 1.0, gearUseStrategyPA = "rand",
                          useSDMs = c("PA", "MsPP", "Gear"),
                          doExp = TRUE, doStats = TRUE, doPlots = TRUE,
                          scenariosPrettyNames=NULL, xAxisTitle="scenarios"){
  
  # Runs the code that I seem to be running for every experiment, i.e.
  #   - create scenarios
  #   - run scenarios
  #   - save results
  #   - create statistics from results
  #   - plot statistics
  # FYI: need to have run the setupSim.r script once before this.  This will create the 
  # simObj, cellsObj, BG, etc. that are used from the global environment in runScenarios.
  
  # Run scenarios of experiment ...
  if ( doExp ) {
    scenariosObj <- initialiseScenarios(scenariosDir, numPA, gammaMultiplier, deltaMultiplier, 
                                        zetaMultiplier, numClusters, numRuns, 
                                        prettyNames = scenariosPrettyNames)
    if ( is.null(simObj) || is.null(simObj$initBeta) ) stop("Please run 'setupSim.r' first.")
    retLst <- runScenarios(scenariosObj, simObj, cellsObj, useSDMs, BG, 
                           domainObj, plotObj, gearUseStrategyPA, NULL,
                           randomSeed = randomSeed)  
    
    # Save results.
    saveScenariosResults(scenariosObj$scenariosDir, retLst, scenariosObj$nameResultsFile)
    
    # Were there errors ...
    deSink()
    showAllScenariosErrors(retLst$resLstAll) 
    #showAllScenariosWarnings(retLst$resLstAll)
    numSuccesses <- numSuccessfulRunsScenarios(scenariosObj, retLst$resLstAll)
  }
  
  # Plot to check gear assignment distribution for PA data.
  if ( doPlots ) {
    # Need to load saved experiment results?
    if ( ! doExp ) {
      # Reload these results
      retLst <- loadScenariosResults(scenariosDir, nameResultsFile = "Results")
      scenariosObj <- retLst$scenariosObj
      scenariosObj$scenariosDir <- scenariosDir
    }
    
    # Plot directory for scenario level plots.
    # plotDir <- paste0(scenariosObj$scenariosDir, "/Plots")
    # if ( ! dir.exists(plotDir) ) {
    #   dir.create(plotDir, recursive = TRUE)
    # }
    # thisFileName <- makeFileName("gearDist", plotDir, plotObj$device)
    
    # Load data.
    nSc <- scenariosObj$numScenarios
    load(paste0(scenariosDir, "/", scenariosObj$namesScenarios[nSc], "/DataDumps/DataDump-AllRuns.RData"))
    load(paste0(scenariosDir, "/", scenariosObj$namesScenarios[nSc], "/DataDumps/DataExample-run1.RData"))
    
    # Covariate setup stuff.
    nameCovar <- names(cellObj$covars)[1]
    cellCovar <- cellObj$covars[ ,nameCovar]
    rangeCovar <- range(cellCovar)
    sampleCovar <- cellObj$covars[surveyObj$rowsInCells,nameCovar]
    
    # Start plot.
    opar <- par(mfrow=c(2,2))
    plot(density(cellCovar))
    points(density(sampleCovar), pch=".", col="red")
    for ( g in 1:surveyObj$numGears ) {
      indGear <- surveyObj$gears == g
      hist(sampleCovar[indGear], main=paste0("gear ", g), xlim=rangeCovar, ylim=c(0,700))
    }
    par(opar)
  }
  
  # Plot summary statistics ...
  if ( doStats ) {
    # Need to load saved experiment results?
    if ( ! doExp ) {
      # Reload these results
      retLst <- loadScenariosResults(scenariosDir, nameResultsFile = "Results")
      scenariosObj <- retLst$scenariosObj
      scenariosObj$scenariosDir <- scenariosDir
      retLst$scenariosObj <- scenariosObj             # Just in case I have renamed directory since saving results!
    }
    
    # Calculate statistics.
    stats <- makeScenarioStats(scenariosObj, retLst$resLstAll, whichStats = c(1,2,3),
                               whichSDMs=useSDMs, useContrasts = FALSE)
    retLst <- saveScenariosSummaryStats(stats, retLst)
  }
    
  if ( doPlots ) {
    # Need to load saved results?
    if ( ! doExp && ! doStats ) {
      # Reload these results
      retLst <- loadScenariosResults(scenariosDir, nameResultsFile = "Results")
      scenariosObj <- retLst$scenariosObj
      scenariosObj$scenariosDir <- scenariosDir
      stats <- retLst$statsObj$stats

      # Check stats for this experiment have been calculated.
      if ( is.null(retLst$statsObj) ) 
        stop("Unable to plot statistics for given experiment as no saved statistics.")
      
    } else if ( ! doExp && doStats ) { 
      # Do nothing, results loaded in statistics step and stats calculated.  
      
    } else if ( doExp && ! doStats ) {
      stop("Unable to plot statistics as statistics have not been calculated for this experiment.")
      
    } else {
      # Do nothing, experiment and statistics steps both done.
    }
    
    # Plot statistics.
    plotScenariosSummaryLambda(stats[[1]], "Alpha", whichCoeffs = 1, plotSDMs = useSDMs,
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle, outlierLabels = TRUE)
    plotScenariosSummaryLambda(stats[[1]], "Beta1", whichCoeffs = 2, plotSDMs = useSDMs,
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle, outlierLabels = TRUE)
    plotScenariosSummaryLambda(stats[[1]], "Beta2", whichCoeffs = c(3), plotSDMs = useSDMs,
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle, outlierLabels = TRUE)
    plotScenariosSummaryLambda(stats[[1]], "Beta3", whichCoeffs = c(4), plotSDMs = useSDMs,
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle, outlierLabels = TRUE)
    plotScenariosSummaryLambda(stats[[1]], "Zetas", whichCoeffs = 6:7, plotSDMs = c("PA","Gear"),
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle)
    plotScenariosSummaryLambda(stats[[1]], "Gamma", whichCoeffs = 8, plotSDMs = c("MsPP","Gear"),
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle)
    plotScenariosSummaryLambda(stats[[1]], "Delta", whichCoeffs = 9, plotSDMs = c("MsPP","Gear"),
                               diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle, vioPlots = FALSE, ylimVals = c(0,1))

    plotScenariosSummaryLambda(stats[[2]], "Intensity",
                               plotSDMs = useSDMs, vioPlots = TRUE, plotDevice="png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle)
    plotScenariosSummaryLambda(stats[[3]], "Total Abundance",
                               plotSDMs = useSDMs, vioPlots = TRUE, plotDevice="png",
                               plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                               xAxisTitle = xAxisTitle)

    whichExperiment <- scenariosObj$scenariosDir
    plotStatisticsComparisons(whichExperiment, whichStats=c(1,1,1,2,3),
                              whichCoeffs=list(1,2,c(3:4),NULL,NULL), whichSpecies=NULL,
                              plotSDMs=useSDMs, xAxisTitle = xAxisTitle,
                              plotDevice=plotObj$device, plotWidth = 16.8, #20
                              plotDir=paste0(whichExperiment,"/Plots"), 
                              columnHeadings = c(expression(alpha[k]),
                                                 expression(beta[1*k]), 
                                                 expression(list(beta[2*k],beta[3*k])),
                                                 expression(lambda[ck]),
                                                 expression(Sigma[c]*lambda[ck])))
    
    # Check estimated coefficients against what should have been estimated (difference from true that we expect!).
    # stats <- makeScenarioStats(scenariosObj, retLst$resLstAll, whichStats = 1,
    #                            whichSDMs=useSDMs, useContrasts = TRUE)
    # plotScenariosSummaryLambda(stats[[1]], "Alpha-cont", whichCoeffs = 1, plotSDMs = useSDMs,
    #                            diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
    #                            plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
    #                            xAxisTitle = xAxisTitle)
    # plotScenariosSummaryLambda(stats[[1]], "Alpha+Zeta-cont", whichCoeffs = c(1,6,7), plotSDMs = c("PA","Gear"),
    #                            diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
    #                            plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
    #                            xAxisTitle = xAxisTitle, outlierLabels = TRUE)
    # plotScenariosSummaryLambda(stats[[1]], "Gamma-cont", whichCoeffs = 8, plotSDMs = c("MsPP","Gear"),
    #                            diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
    #                            plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
    #                            xAxisTitle = xAxisTitle)
  }

  # Return value.
  if ( !doExp ) numSuccesses = "unknown"
  return(list(scenariosObj = scenariosObj, 
              retLst = retLst,
              numSuccesses = numSuccesses))
  
}

#-----------------------------------------------------------------------------------------
