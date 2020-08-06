# Set working directory.
setwd("/perm_storage/home/sampeel/chap2/sim2")

# Load packages and source files.
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


# Re-load results (if needed)
retLst <- loadScenariosResults(scenariosDir, nameResultsFile = "Results")
scenariosObj <- retLst$scenariosObj
scenariosObj$scenariosDir <- scenariosDir



# #---------------------------
# # Figure 1 - Scale (intercept) issue
# #---------------------------
# 
# # Setup ...
# 
# numRuns <- 5000
# scenariosDir <- paste0(getwd(), "/", "Output-SpatialGear-runs", format(numRuns))
# gearUseStrategyPA <- "covar"
# numPA <- 3500
# numClusters <- 0
# deltaMultiplier <- 0.0
# gammaMultiplier <- 1.0
# zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
# scenariosPrettyNames <- as.character(format(zetaMultiplier))
# xAxisTitle <- expression("gear difference " * (zeta^x)) 
# 
# # Run experiment 
# expLst <- runExperiment(scenariosDir, numRuns, 
#                         gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
#                         numClusters = numClusters, deltaMultiplier = deltaMultiplier,
#                         gammaMultiplier = gammaMultiplier,
#                         zetaMultiplier = zetaMultiplier,
#                         scenariosPrettyNames = scenariosPrettyNames,
#                         xAxisTitle = xAxisTitle)
# 
# 
# # Plot figure.
# whichExperiment <- paste0(getwd(), "/", "Output-SpatialGear-runs", format(numRuns))
# xAxisTitle <- expression("gear difference " * (zeta^x)) 
# plotStatisticsComparisons(whichExperiment, whichStats=c(1,2,3),
#                           whichCoeffs = list(1,NULL,NULL), whichSpecies=NULL,
#                           plotSDMs = c("PA","MsPP","Gear"), 
#                           xAxisTitle = xAxisTitle,
#                           plotDevice = "png", plotWidth = 12, plotHeight = 12,
#                           plotDir=paste0(whichExperiment,"/Plots"), fileName = "figure1",
#                           columnHeadings = c(expression(alpha[k]),
#                                              expression(lambda[ck]),
#                                              expression(Sigma[c]*lambda[ck])),
#                           accuracyPrecision = "accuracy")


#---------------------------
# Appendix - number of samples
#---------------------------

# Experiment 1 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-RandomGear-nPA3500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 2 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-RandomGear-nPA2500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 2500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 3 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-RandomGear-nPA1500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 1500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 4 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-RandomGear-nPA500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Make figure 2.
experimentDirs <- c(paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA3500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA2500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA1500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA500-runs5000"))
experimentNames <- c(expression(n[samp] == 3500),expression(n[samp] == 2500),
                     expression(n[samp] == 1500),expression(n[samp] == 500))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "randNumPAComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")


# Make figure 1.
whichExperiment <- paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA3500-runs5000")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotDir <- paste0(getwd(),"/Paper/Main")
plotStatisticsComparisons(whichExperiment, whichStats=c(1,2,3),
                          whichCoeffs = list(1,NULL,NULL), whichSpecies=NULL,
                          plotSDMs = c("PA","MsPP","Gear"), 
                          xAxisTitle = xAxisTitle,
                          plotDevice = "png", plotWidth = 12, plotHeight = 12,
                          plotDir = plotDir, fileName = "scaleComparisonFig1",
                          columnHeadings = c(expression(alpha[k]),
                                             expression(lambda[ck]),
                                             expression(Sigma[c]*lambda[ck])),
                          accuracyPrecision = "accuracy")


# Experiment 5 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-SpatialGear-nPA3500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Experiment 6 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-SpatialGear-nPA2500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 2500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Experiment 7 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-SpatialGear-nPA1500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 1500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Experiment 8 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumSamps/", "Output-SpatialGear-nPA500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 500
numClusters <- 0
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Make figure 3
experimentDirs <- c(paste0(getwd(),"/Paper/AppNumSamps/Output-SpatialGear-nPA3500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-SpatialGear-nPA2500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-SpatialGear-nPA1500-runs5000"),
                    paste0(getwd(),"/Paper/AppNumSamps/Output-SpatialGear-nPA500-runs5000"))
experimentNames <- c(expression(n[samp] == 3500),expression(n[samp] == 2500),
                     expression(n[samp] == 1500),expression(n[samp] == 500))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "spatialNumPAComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")


#---------------------------
# App - Sampling bias
#---------------------------


# Experiment 1 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-RandomGear-Bias0-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 0.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


# Experiment 2 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-RandomGear-Bias1-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)



# Experiment 3 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-RandomGear-Bias2-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 2.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


# Make combined figure.
experimentDirs <- c(paste0(getwd(),"/Paper/AppSampBias/Output-RandomGear-Bias0-runs5000"),
                    paste0(getwd(),"/Paper/AppSampBias/Output-RandomGear-Bias1-runs5000"),
                    paste0(getwd(),"/Paper/AppSampBias/Output-RandomGear-Bias2-runs5000"))
experimentNames <- c(expression(delta^x == 0),expression(delta^x == 1),expression(delta^x == 2))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 12, plotHeight = 12, 
                          fileName = "randomSampBiasComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")



# Experiment 4 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-SpatialGear-Bias0-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 0.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Experiment 5 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-SpatialGear-Bias1-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)



# Experiment 6 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSampBias/Output-SpatialGear-Bias2-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
gammaMultiplier <- 1.1
deltaMultiplier <- 2.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, useSDMs = c("PA","MsPP","Gear"),
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA,
                        numClusters = numClusters, gammaMultiplier = gammaMultiplier,  
                        deltaMultiplier = deltaMultiplier, zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)


# Make combined figure.
experimentDirs <- c(paste0(getwd(),"/Paper/AppSampBias/Output-SpatialGear-Bias0-runs5000"),
                    paste0(getwd(),"/Paper/AppSampBias/Output-SpatialGear-Bias1-runs5000"),
                    paste0(getwd(),"/Paper/AppSampBias/Output-SpatialGear-Bias2-runs5000"))
experimentNames <- c(expression(delta^x == 0),expression(delta^x == 1),expression(delta^x == 2))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 12, plotHeight = 12, 
                          fileName = "spatialSampBiasComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")


#---------------------------
# App multiple surveys
#---------------------------

# Experiment 1 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-RandomGearClust-nPA500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 2 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-RandomGearClust-nPA1500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 1500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 3 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-RandomGearClust-nPA2500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 2500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 4 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-RandomGearClust-nPA3500-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppSurveys/Output-RandomGearClust-nPA3500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-RandomGearClust-nPA2500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-RandomGearClust-nPA1500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-RandomGearClust-nPA500-runs5000"))
experimentNames <- c(expression(n[samp] == 3500),expression(n[samp] == 2500),
                     expression(n[samp] == 1500),expression(n[samp] == 500))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "randNumPAwClustComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")


# Experiment 5 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-SpatialGearClust-nPA500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 6 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-SpatialGearClust-nPA1500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 1500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 7 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-SpatialGearClust-nPA2500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 2500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Experiment 8 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppSurveys/Output-SpatialGearClust-nPA3500-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 0.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp = FALSE, doStats = FALSE)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppSurveys/Output-SpatialGearClust-nPA3500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-SpatialGearClust-nPA2500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-SpatialGearClust-nPA1500-runs5000"),
                    paste0(getwd(),"/Paper/AppSurveys/Output-SpatialGearClust-nPA500-runs5000"))
experimentNames <- c(expression(n[samp] == 3500),expression(n[samp] == 2500),
                     expression(n[samp] == 1500),expression(n[samp] == 500))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "spatialNumPAwClustComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")

#---------------------------
# Number of PO random (with 2 x bias) 
#---------------------------

# Experiment 1 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-RandomGearBias-nPO10-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 2 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-RandomGearBias-nPO11-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle, doExp=FALSE, doStats=FALSE)

# Experiment 3 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-RandomGearBias-nPO12-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.2
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 4 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-RandomGearBias-nPO13-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppNumPOObs/Output-RandomGearBias-nPO10-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-RandomGearBias-nPO11-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-RandomGearBias-nPO12-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-RandomGearBias-nPO13-runs5000"))
experimentNames <- c(expression(delta^x == 0),expression(delta^x == 1),expression(delta^x == 2))
experimentNames <- c(expression(gamma^x == 1.0),expression(gamma^x == 1.1),
                     expression(gamma^x == 1.2),expression(gamma^x == 1.3))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "randNumPOComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")

# Experiment 5 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-SpatialGearBias-nPO10-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 6 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-SpatialGearBias-nPO11-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.1
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                     gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                     numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                     gammaMultiplier = gammaMultiplier,
                     zetaMultiplier = zetaMultiplier,
                     scenariosPrettyNames = scenariosPrettyNames,
                     xAxisTitle = xAxisTitle)

# Experiment 7 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-SpatialGearBias-nPO12-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.2
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 8 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppNumPOObs/Output-SpatialGearBias-nPO13-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 0
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppNumPOObs/Output-SpatialGearBias-nPO10-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-SpatialGearBias-nPO11-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-SpatialGearBias-nPO12-runs5000"),
                    paste0(getwd(),"/Paper/AppNumPOObs/Output-SpatialGearBias-nPO13-runs5000"))
experimentNames <- c(expression(delta^x == 0),expression(delta^x == 1),expression(delta^x == 2))
experimentNames <- c(expression(gamma^x == 1.0),expression(gamma^x == 1.1),
                     expression(gamma^x == 1.2),expression(gamma^x == 1.3))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "spatialNumPOComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")

#---------------------------
# Ratio of PA to PO (with 2 x bias & 10 surveys) 
#---------------------------

# Experiment 1 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-RandomGearBias-hPAhPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 2 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-RandomGearBias-lPAhPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 3 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-RandomGearBias-hPAlPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)


# Experiment 4 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-RandomGearBias-lPAlPO-runs", format(numRuns))
gearUseStrategyPA <- "rand"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppRatio/Output-RandomGearBias-hPAhPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-RandomGearBias-lPAhPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-RandomGearBias-hPAlPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-RandomGearBias-lPAlPO-runs5000"))
experimentNames <- c(expression(list(n[samp] == 3500,gamma^x == 1.0)),
                     expression(list(n[samp] == 500,gamma^x == 1.0)),
                     expression(list(n[samp] == 3500,gamma^x == 1.3)),
                     expression(list(n[samp] == 500,gamma^x == 1.3)))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "randRatioComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")

# Experiment 5 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-SpatialGearBias-hPAhPO-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 6 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-SpatialGearBias-lPAhPO-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.0
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 7 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-SpatialGearBias-hPAlPO-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 3500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Experiment 8 Setup ...
numRuns <- 5000
scenariosDir <- paste0(getwd(), "/Paper/AppRatio/Output-SpatialGearBias-lPAlPO-runs", format(numRuns))
gearUseStrategyPA <- "covar"
numPA <- 500
numClusters <- 10
deltaMultiplier <- 2.0
gammaMultiplier <- 1.3
zetaMultiplier <- seq(from=0.0, to=2.0, by=0.5)
scenariosPrettyNames <- as.character(format(zetaMultiplier))
xAxisTitle <- expression("gear difference " * (zeta^x)) 

# Run experiment 
expLst <- runExperiment(scenariosDir, numRuns, 
                        gearUseStrategyPA = gearUseStrategyPA, numPA = numPA, 
                        numClusters = numClusters, deltaMultiplier = deltaMultiplier,
                        gammaMultiplier = gammaMultiplier,
                        zetaMultiplier = zetaMultiplier,
                        scenariosPrettyNames = scenariosPrettyNames,
                        xAxisTitle = xAxisTitle)

# Make combined plot.
experimentDirs <- c(paste0(getwd(),"/Paper/AppRatio/Output-SpatialGearBias-hPAhPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-SpatialGearBias-lPAhPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-SpatialGearBias-hPAlPO-runs5000"),
                    paste0(getwd(),"/Paper/AppRatio/Output-SpatialGearBias-lPAlPO-runs5000"))
experimentNames <- c(expression(list(n[samp] == 3500,gamma^x == 1.0)),
                     expression(list(n[samp] == 500,gamma^x == 1.0)),
                     expression(list(n[samp] == 3500,gamma^x == 1.3)),
                     expression(list(n[samp] == 500,gamma^x == 1.3)))
columnHeadings <- experimentNames
plotDir <- paste0(getwd(),"/Paper/Main")
xAxisTitle <- expression("gear difference " * (zeta^x)) 
plotTheseSDMs = c("PA","MsPP","Gear")

# Statistic 1b: difference between true and estimated beta_1 coefficient.
plotExperimentComparisons(experimentDirs, 1, whichCoeffs = c(2), plotSDMs = plotTheseSDMs,
                          plotDevice = "png", plotDir = plotDir, #ylimVals = c(0,0.5,0,5), 
                          plotWidth = 15.8, plotHeight = 12, # plotWidth = 12, plotHeight = 12, 
                          fileName = "spatialRatioComparisonBeta1", xAxisTitle = xAxisTitle,  
                          columnHeadings = columnHeadings, accuracyPrecision = "accuracy")

#--------------------------------------------------------------
# Map plots to demonstrate gear effect for one species example!
#--------------------------------------------------------------

### Plot to show difference between true species dist and MsPP SDM estimated (but without scale).

# Look at just sample bias experiment (with no other degradations).
spatialGearExpDir <- paste0(getwd(), "/Output-SpatialGear-runs5000")

# Load results and plotting object (just simpler to use functions already written for this object)
tmp <- load(paste0(spatialGearExpDir,"/Data/Results.RData"))

# Worst species for beta1 coefficient
worstSpecies <- names(sort(abs(statsObj$stats[[1]]$avg["2.0", ,"MsPP",2]),decreasing = TRUE))[1]

# Get coastline polygon.
load(paste0("../Data/Coastlines/southOceanPoly.RData"))
coastPoly <- cm1
rm(cm1)
if ( ! compareCRS(coastPoly@proj4string, simObj$proj) ) {
  # Project first, then set.
  coastPoly <- spTransform(coastPoly, simObj$proj)
}

# Make figure. NB: this is different to figure in draft 2 as have added exp(meanEstAlpha) to cell values.
plotNames <- c(expression(zeta^x == 0.0), expression(zeta^x == 1.0), expression(zeta^x == 2.0))
rlPlot <- plotScenarioMaps(spatialGearExpDir, worstSpecies, "MsPP", c("ZMult0.0","ZMult1.0","ZMult2.0"), 
                 plotDevice = "png", plotDir = paste0(spatialGearExpDir,"/Plots"), 
                 plotWidth = 10.7, coastLine = coastPoly,  colLand="lightgrey", 
                 nameLand = "Antarctica", posNameLand = c(-700,-100), useValues="intensity", 
                 fileName = "mapDiffTrueSpatialGearMsPP", plotNames = plotNames, 
                 plotDiff=TRUE, useAlpha=FALSE, colNA = "slategray1", zlim=c(-1.0,1.0))



#--------------------------------------------------------------
# Random and biased gear assignment example plot.
#--------------------------------------------------------------

# Data files to use for example gear assignment plots.
dataFile <- paste0("/ZMult0.0/DataDumps/DataExample-run1.RData")
expDirs <- c(paste0(getwd(),"/Paper/AppNumSamps/Output-RandomGear-nPA3500-runs5000"),
             paste0(getwd(),"/Paper/AppNumSamps/Output-SpatialGear-nPA3500-runs5000"),
             paste0(getwd(),"/Paper/AppSurveys/Output-RandomGearClust-nPA3500-runs5000"),
             paste0(getwd(),"/Paper/AppSurveys/Output-SpatialGearClust-nPA3500-runs5000"))

plotTitles <- c("Random gear assignment (Exp No. 1.R.1)",
                "Biased gear assignment (Exp No. 1.B.1)",
                "Random gear assignment (Exp No. 3.R.1)",
                "Biased gear assignment (Exp No. 3.B.1)")

# Plot settings
nPlots <- length(expDirs)
colours <- colourBlindRGB()
if ( length(expDirs) >= length(colours) ) stop("Not enough colours defined to plot this many lines on one plot.")
pchSymbols <- c(NA,1:(length(expDirs)-1))
ltyLines <- rep("solid", times=length(colours))
yRange <- NULL
nGears <- 3              
densities <- array(data=vector(mode="list",length=(nGears+1)*nPlots), 
                   dim = c(nGears+1, nPlots), 
                   dimnames = list(c("All", paste0("g",1:3)),paste0("plot",1:nPlots)))
maxDensityGear <- 0
nGearPlot <- matrix(0, nrow=nGears, ncol=nPlots)

for ( p in 1:nPlots ) {
  # Load and format sample info ...
  tmp <- load(paste0(expDirs[p],dataFile))
  xy <- surveyObj$xy
  gears <- surveyObj$gears
  nGears <- surveyObj$numGears
  surveyBath <- cellObj$covars[surveyObj$rowsInCells,1]

  # Get the density function for all the environment values.
  allBath <- cellObj$covars[ ,1]
  densities[1,p][[1]] <- density(allBath) 
  
  # Get the densities for the survey environment values per gear.
  for ( g in 1:nGears ) {
    indGear <- which(gears == g)
    nGearPlot[g,p] <- length(indGear)
    densities[g+1,p][[1]] <- density(surveyBath[indGear])
    maxDensityGear <- max(densities[g+1,p][[1]]$y, maxDensityGear)
  }

  # Get the required range for all plots.
  xRange <- range(densities[1,p][[1]]$x)
  yRange <- range(c(densities[1,p][[1]]$y, maxDensityGear, yRange))
}

# Which device are we printing to?
plotDevice = "png"
if ( plotDevice == "RStudioGD" ) {
  # plot to the R studio plot window.
  plotToFile <- FALSE
} else {
  fileName <- makeFileName("GearTypeDensities", paste0(getwd(),"/Paper/Main"), plotDevice)
  argList <- list(filename = fileName, width = 15.75, height = 15.75, units = "cm", res = 600)
  do.call(plotDevice, argList)
  plotToFile <- TRUE
}

# Plots ...
opar <- par(mfrow = c(2,2), mar=c(3.6,3.6,2.1,0.5), cex=0.66, cex.main=0.9)
for ( p in 1:nPlots ) {
  # Plot bathymetry density.
  plot(densities[1,p][[1]]$x, densities[1,p][[1]]$y, 
       ylim=yRange, type="l", col=colours[1], 
       lty=ltyLines[1], main=plotTitles[p], xlab="", ylab="")
  title(xlab="centred bathymetry", ylab = "density", line=2.5)
  abline(list(h=0.0), lty="dotted")

  # Add gear densities.
  for ( g in 1:nGears ) {
    x <- densities[g+1,p][[1]]$x
    y <- densities[g+1,p][[1]]$y
    lines(x, y, col=colours[g+1], lty=ltyLines[g+1])
    
    # Plot symbols on lines (not all as too messy)
    numX <- length(x)
    pts <- seq(from=1, to=numX, by=50) + g*10  # offset so symbols aren't all in same place.
    numPts <- length(pts)
    pts[numPts] <- min(pts[numPts], numX)      # make sure the last value isn't too large.
    points(x[pts], y[pts], col=colours[g+1], pch=pchSymbols[g+1])
  }

  # Add legend.
  legend("topleft", legend=c("all cells", paste0("gear ", 1:nGears)), 
         col=colours[1:(nGears+1)], lty=ltyLines[1:(nGears+1)], pch=pchSymbols[1:(nGears+1)])
  
  # Add second legend.
  # legend("topright", legend=paste0("gear ", 1:nGears, " = ", nGearPlot[ ,p]), bty="n",
  #        title = "No. samples using:")
}  
par(opar)  
if ( plotToFile ) dev.off()

