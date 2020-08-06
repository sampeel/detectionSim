# ***************
#
# If have not already run setup for the simulation (i.e. simObj doesn't exist), 
# then use script in 'setupSim.r' to do this before using this script.
#
# ***************

#----------------------------------
# General Setup for experiments ...
#----------------------------------

highNumPA <- 3500L                       # Lots of survey locations
lowNumPA <- 500L                         # Not many survey locations
highGammaMultiplier <- 1.0               # More PO points per species
lowGammaMultiplier <- 1.3                # Less PO points per species
highDeltaMultiplier <- 2.0               # Lots of sample bias
lowDeltaMultiplier <- 0.0                # No sample bias
highNumClusters <- 500L                  # Slightly clustered survey locations
# (NB: no clustering or random sample locations is numClusters = 0)
lowNumClusters <- 20L                    # Highly clustered survey locations
tinyNumClusters <- 5L                    # Extremely clustered survey locations
testNumRuns <- 10L    
smallNumRuns <- 100L
lowNumRuns <- 1000L
highNumRuns <- 5000L
hugeNumRuns <- 10000L
truncRange <- c(0,1)                    # Truncate the y-axis to this range on summary plots
randomSeed <- 10


#---------------------------
# Experiment 0.0 - Benchmark - test MsPP and gearGLM perform same without gear!
#---------------------------

### Setup ...

gearUseStrategyPA <- "rand"
gearUseStrategyPO <- NULL
deltaMultiplier <- seq(from=lowDeltaMultiplier, to=highDeltaMultiplier, by=0.5)
numClusters <- as.integer(0)
numPA <- highNumPA      #as.integer(seq(from=500, to=5000, by=500)) 
gammaMultiplier <- 1.0
zetaMultiplier <- 1.0
numRuns <- 100
scenariosDir <- paste0(getwd(), "/", "Output-BenchmarkBias-runs", format(numRuns))
scenariosPrettyNames <- as.character(format(deltaMultiplier))
xAxisTitle <- "sample bias"

# Run scenarios ...
scenariosObj <- initialiseScenarios(scenariosDir, numPA, gammaMultiplier, deltaMultiplier, 
                                    zetaMultiplier, numClusters, numRuns, 
                                    prettyNames = scenariosPrettyNames)
simObjTmp <- simObj
simObjTmp$gear.formula <- NULL  # turn off gear stuff in gearGLM!
retLst <- runScenarios(scenariosObj, simObjTmp, cellsObj, c("PA","MsPP", "Gear"), BG, 
                       domainObj, plotObj, gearUseStrategyPA, gearUseStrategyPO,
                       randomSeed = randomSeed)  

# Save results.
saveScenariosResults(scenariosObj$scenariosDir, retLst, scenariosObj$nameResultsFile)

# Were there errors ...
deSink()
showAllScenariosErrors(retLst$resLstAll) 
showAllScenariosWarnings(retLst$resLstAll)
numSuccessfulRunsScenarios(scenariosObj, retLst$resLstAll)

# Plot summary statistics ...
stats <- makeScenarioStats(scenariosObj, retLst$resLstAll, whichStats = c(1,2,3),
                           whichSDMs=c("MsPP","Gear"), useContrasts = FALSE)
retLst <- saveScenariosSummaryStats(stats, retLst)
plotScenariosSummaryLambda(stats[[1]], "Alpha", whichCoeffs = 1, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle)
plotScenariosSummaryLambda(stats[[1]], "Betas", whichCoeffs = c(2:4), plotSDMs = c("MsPP","Gear"),
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
                           plotSDMs = c("MsPP","Gear"), vioPlots = TRUE, plotDevice="png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, horizontalLines = c(1,0))
plotScenariosSummaryLambda(stats[[3]], "Total Abundance",
                           plotSDMs = c("MsPP","Gear"), vioPlots = TRUE, plotDevice="png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle)

whichExperiment <- scenariosObj$scenariosDir
plotStatisticsComparisons(whichExperiment, whichStats=c(1,1,1,2,3), 
                          whichCoeffs=list(1,2,c(3:4),NULL,NULL), whichSpecies=NULL,
                          plotSDMs=c("MsPP","Gear"), xAxisTitle = xAxisTitle,
                          plotDevice=plotObj$device, plotWidth = 20,
                          plotDir=paste0(whichExperiment,"/Plots"), 
                          columnHeadings = c(expression(alpha[k]),
                                             expression(beta[1*k]), 
                                             expression(list(beta[2*k],beta[3*k])),
                                             expression(lambda[ik]),
                                             expression(Sigma[i]*lambda[ik])))


#---------------------------
# Experiment 0.1 - Benchmark - test MsPP and gearGLM perform same without gear!
#                              looking at possible delta bias in MsPP that is affecting 
#                              certain coeffs worse than others (i.e. Beta2-sp27, Beta3-sp22)
#---------------------------

### Setup ...

numRuns <- 5000
gearUseStrategyPA <- "rand"
gearUseStrategyPO <- NULL
deltaMultiplier <- 0.0               # no bias
numClusters <- as.integer(0)         # random samples
numPA <- highNumPA                   #as.integer(seq(from=500, to=5000, by=500)) 
gammaMultiplier <- seq(from=1.0, to=1.3, by=0.05)
zetaMultiplier <- 0.0                # no gear effect
scenariosDir <- paste0(getwd(), "/", "Output-BenchmarkNumPO-runs", format(numRuns))
scenariosPrettyNames <- c("13500","8000","4800","2900","1700","1000","630")
xAxisTitle <- "approx num PO"

# Run scenarios ...
scenariosObj <- initialiseScenarios(scenariosDir, numPA, gammaMultiplier, deltaMultiplier, 
                                    zetaMultiplier, numClusters, numRuns, 
                                    prettyNames = scenariosPrettyNames)
simObjTmp <- simObj
simObjTmp$gear.formula <- NULL  # turn off gear stuff in gearGLM!
retLst <- runScenarios(scenariosObj, simObjTmp, cellsObj, c("MsPP", "Gear"), BG, 
                       domainObj, plotObj, gearUseStrategyPA, gearUseStrategyPO,
                       randomSeed = randomSeed)  

# Save results.
saveScenariosResults(scenariosObj$scenariosDir, retLst, scenariosObj$nameResultsFile)

# Were there errors ...
deSink()
showAllScenariosErrors(retLst$resLstAll) 
showAllScenariosWarnings(retLst$resLstAll)
numSuccessfulRunsScenarios(scenariosObj, retLst$resLstAll)

# Plot summary statistics ...
stats <- makeScenarioStats(scenariosObj, retLst$resLstAll, whichStats = c(1,2,3),
                           whichSDMs=c("MsPP","Gear"), useContrasts = FALSE)
retLst <- saveScenariosSummaryStats(stats, retLst)
plotScenariosSummaryLambda(stats[[1]], "Alpha", whichCoeffs = 1, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle)
plotScenariosSummaryLambda(stats[[1]], "Beta1", whichCoeffs = 2, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, outlierLabels = TRUE)
plotScenariosSummaryLambda(stats[[1]], "Beta2", whichCoeffs = 3, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, outlierLabels = TRUE)
plotScenariosSummaryLambda(stats[[1]], "Beta3", whichCoeffs = 4, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, outlierLabels = TRUE)
plotScenariosSummaryLambda(stats[[1]], "Gamma", whichCoeffs = 8, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle)
plotScenariosSummaryLambda(stats[[1]], "Delta", whichCoeffs = 9, plotSDMs = c("MsPP","Gear"),
                           diffPORange = c(TRUE,FALSE), bigNumThreshold = 100, plotDevice = "png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, vioPlots = FALSE, ylimVals = c(0,1))

plotScenariosSummaryLambda(stats[[2]], "Intensity",
                           plotSDMs = c("MsPP","Gear"), vioPlots = TRUE, plotDevice="png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle, horizontalLines = c(1,0))
plotScenariosSummaryLambda(stats[[3]], "Total Abundance",
                           plotSDMs = c("MsPP","Gear"), vioPlots = TRUE, plotDevice="png",
                           plotDir = paste0(scenariosObj$scenariosDir, "/Plots"),
                           xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 1.1 - Random Gear - increasing gear difference
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomGear-runs", format(numRuns))
gearUseStrategyPA <- "rand"
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns,
                        gearUseStrategyPA = gearUseStrategyPA, 
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Reload these results, if necessary.
scenariosDir <- paste0(getwd(), "/", "Output-RandomGear-runs", format(numRuns))
retLst <- loadScenariosResults(scenariosDir, nameResultsFile = "Results")
scenariosObj <- retLst$scenariosObj
scenariosObj$scenariosDir <- scenariosDir


#---------------------------
# Experiment 1.2 - Random gear - decreasing sample number
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomNumPA-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- seq(from=highNumPA, to=lowNumPA, by=-500)
scenariosPrettyNames <- as.character(format(numPA))
xAxisTitle <- "num samples"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        numPA = numPA,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


#---------------------------
# Experiment 1.3 - Random gear - decreasing survey number
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomNumSurv-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numClusters <- c(0, 500, 100, 50, 20, 10, 5)
scenariosPrettyNames <- c(paste0("r",numPA), as.character(format(numClusters[-1])))
xAxisTitle <- "num surveys"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        numClusters = numClusters,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 1.4 - Random gear - decreasing numPO
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomNumPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
gammaMultiplier <- seq(from=1.0, to=1.3, by=0.05)
scenariosPrettyNames <- c("13500","8000","4800","2900","1700","1000","630")
xAxisTitle <- "approx num PO"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        gammaMultiplier = gammaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)



#---------------------------
# Experiment 1.5 - Random gear - increasing sample bias
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomBias-runs", format(numRuns))
gearUseStrategyPA <- "rand"
deltaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(deltaMultiplier))
xAxisTitle <- "sample bias multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        deltaMultiplier = deltaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 2.1 - Spatial gear - increasing gear difference
#---------------------------

# Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialGear2-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- highNumPA
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 2.2 - Spatial gear - decreasing sample number
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialNumPA-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- seq(from=highNumPA, to=lowNumPA, by=-500)
scenariosPrettyNames <- as.character(format(numPA))
xAxisTitle <- "num samples"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        numPA = numPA,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 2.3 - Spatial - decreasing survey number
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialNumSurv-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numClusters <- c(0, 500, 100, 50, 20, 10, 5)
scenariosPrettyNames <- c(paste0("r",highNumPA), as.character(format(numClusters[-1])))
xAxisTitle <- "num surveys"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        numClusters = numClusters,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

#---------------------------
# Experiment 2.4 - Spatial - decreasing numPO
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialNumPO-runs", format(numRuns))
gearUseStrategyPA <- "covar"
gammaMultiplier <- seq(from=1.0, to=1.3, by=0.05)
scenariosPrettyNames <- as.character(format(gammaMultiplier))
xAxisTitle <- "gamma multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        gammaMultiplier = gammaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

#---------------------------
# Experiment 2.5 - Spatial - increasing sample bias
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialBias-runs", format(numRuns))
gearUseStrategyPA <- "covar"
deltaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(deltaMultiplier))
xAxisTitle <- "sample bias multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, 
                        deltaMultiplier = deltaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

#---------------------------
# Experiment 3.1 - Random gear - increasing gear effect
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomAll-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- lowNumPA
numClusters <- 10
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 3.2 - Spatial gear - increasing gear effect
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialAll-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- lowNumPA
numClusters <- 10
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 4.1 - Random gear - increasing gear effect (without clustering)
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomNoClust-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- lowNumPA
numClusters <- 0
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0   
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 4.2 - Random gear - increasing gear effect (without clustering or lowPA)
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomLowPOBias-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- highNumPA
numClusters <- 0
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0   
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- "gear effect multiplier"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 4.3 - Random gear - what number of samples will work?
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomAllwrtNumPA-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- seq(from=1500, to=500, by=-200)
numClusters <- 10
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0
zetaMultiplier <- 1.0
scenariosPrettyNames <- as.character(format(numPA))
xAxisTitle <- "num samples"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

#---------------------------
# Experiment 4.4 - Random gear - what number of bigger samples will work?
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomAllwrtNumPABig-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- seq(from=3500, to=1500, by=-500)
numClusters <- 10
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0
zetaMultiplier <- 1.0
scenariosPrettyNames <- as.character(format(numPA))
xAxisTitle <- "num samples"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 4.5 - Spatial gear - what number of samples will work?
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-SpatialAllwrtNumPA-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- seq(from=3500, to=500, by=-500)
numClusters <- 10
gammaMultiplier <- lowGammaMultiplier
deltaMultiplier <- 1.0
zetaMultiplier <- 1.0
scenariosPrettyNames <- as.character(format(numPA))
xAxisTitle <- "num samples"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Experiment 4.6 - Random gear - will having more PO help?
#---------------------------

### Setup ...

numRuns <- 5000
scenariosDir <- paste0(getwd(), "/", "Output-RandomAllwrtNumPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 500
numClusters <- 10
gammaMultiplier <- seq(from=1.0, to=1.3, by=0.05)
deltaMultiplier <- 1.0
zetaMultiplier <- 1.0
scenariosPrettyNames <- c("15300","9100","5500","3300","2000","1200","720")
xAxisTitle <- "approx num PO"

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


#---------------------------
# Plot comparison beta_1 - numbers of data investigation
#---------------------------


experimentDirs <- c("Output-RandomAllwrtNumPABig-runs5000",
                    "Output-RandomAllwrtNumPA-runs5000",
                    "Output-RandomAllwrtNumPO-runs5000",
                    "Output-SpatialAllwrtNumPABig-runs5000")
experimentNames <- c("Experiment 4.4","Experiment 4.3","Experiment 4.6","Experiment 4.5")
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Output-Combined/Plot/NumDataComparison")
xAxisTitle <- c("num samples","num samples","approx num PO","num samples")
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          fileName = "beta1Coeff", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "both")


#------------------------------
# Test single alpha/beta metric
#------------------------------

species <- "sp253"
coeff <- "alpha"
scenario <- 1
namesSDMs <- expLst$retLst$resLstAll[[1]]$namesValidSDM
validSDMs <- expLst$retLst$resLstAll[[1]]$validSDM
theseResults <- NULL
for ( sdm in validSDMs ) {
  theseResults <- cbind(theseResults, expLst$retLst$resLstAll[[scenario]][[sdm]]$coeffs$beta[coeff,species, ])
}
colnames(theseResults) <- namesSDMs

boxplot(theseResults)
abline(h=expLst$retLst$simObj$initBeta[coeff,species],col="red",lty="dotted")
abline(h= expLst$retLst$simObj$initBeta[coeff,species] + expLst$retLst$simObj$zeta["zeta1",species],col="red",lty="dotted")
title(paste("Estimated values for species =", species, "and coeff =", coeff),sub=paste("nruns =", expLst$scenariosObj$numRuns))


#----------------------
# Test zeta accuaracy
#----------------------

### Check results with simple plot ...

calcZetaStat <- function(estZeta, trueZeta, whichGear=2) {
  
  namesSpecies <- dimnames(estZeta)[[2]]
  estZetaGear <- estZeta[whichGear-1,namesSpecies, ]
  trueZetaGear <- trueZeta[whichGear,namesSpecies] - trueZeta[1,namesSpecies]
  statZetaGear <- apply(abs((estZetaGear - trueZetaGear)), 1, mean)   # Repeats trueZetaGear for each column of estZetaGear
  return(statZetaGear)
  
}

# Create statistics ...
statZeta2PA <- calcZetaStat(retLst$resLstAll[[3]]$Cloglog$coeffs$zeta, simObj$initZeta, 2)
statZeta2Gear <- calcZetaStat(retLst$resLstAll[[3]]$Gear$coeffs$zeta, simObj$initZeta, 2)
statZeta3PA <- calcZetaStat(retLst$resLstAll[[3]]$Cloglog$coeffs$zeta, simObj$initZeta, 3)
statZeta3Gear <- calcZetaStat(retLst$resLstAll[[3]]$Gear$coeffs$zeta, simObj$initZeta, 3)

# Plot ...
numRuns <- dim(diffZeta2)[2]
namesSpecies <- dimnames(diffZeta2)[[1]]
numSpecies <- length(namesSpecies)
plot(c(0.5,4.0), range(c(statZeta2,statZeta3)), type="n")
points(rep(1,21), statZeta2[ ], pch="-")
text(rep(1.2,21), statZeta2[ ], labels = namesSpecies)
points(rep(3,21), statZeta3[ ], pch="-")
text(rep(3.2,21), statZeta3[ ], labels = namesSpecies)
# whichSpecies <- c("sp17","sp97","sp107","sp251")
# points(rep(1,length(whichSpecies)), statZeta2[whichSpecies], pch="-",col="red")
# text(rep(1.2,length(whichSpecies)), statZeta2[whichSpecies], labels = whichSpecies, col="red")
# points(rep(3,length(whichSpecies)), statZeta3[whichSpecies], pch="-",col="red")
# text(rep(3.2,length(whichSpecies)), statZeta3[whichSpecies], labels = whichSpecies, col="red")



f0 <- function() {
  r <- 0
  biases <- cellsObj$biases
  return(r)
}
f1 <- function(cellsObj) {
  r <- 1
  biases <- cellsObj$biases
  return(r)
}

n <- 10000000
startTime <- Sys.time()
for ( i in 1:n ) {
  res <- f0()
}
totalTime0 <- Sys.time() - startTime 

startTime <- Sys.time()
for ( i in 1:n ) {
  res <- f1(cellsObj)
}
totalTime1 <- Sys.time() - startTime 

totalTime0
totalTime1


