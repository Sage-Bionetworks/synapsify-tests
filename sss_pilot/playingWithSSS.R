# playingWithSSS.R

require(randomForest)
require(sss)
require(illuminaHumanv4.db)

## PULL IN THE HER2 RANDOM FOREST MODEL OBJECT
her2ModEnt <- loadEntity(138541)
load(file.path(her2ModEnt$cacheDir, her2ModEnt$files))

## PULL OUT THE 'IMPORTANCE' RANKING OF FEATURES
her2Features <- importance(her2RF)

## SORT AND PULL OUT ALL FEATURES WITH NON-ZERO WEIGHTS
sortHER2Features <- sort(abs(her2Features[ , 1]), decreasing = TRUE)
topFeatLog <- sortHER2Features > 0
topInd <- grep("TRUE", topFeatLog)
topHER2Features <- sortHER2Features[topInd]

## MAPPING ILLUMINA FEATURES TO GENE SYMBOLS
her2GSym <- as.character(mget(names(topHER2Features), illuminaHumanv4SYMBOL, ifnotfound = NA))

## FOR DEMO PURPOSES, LET'S SELECT THE TOP 500 FEATURES
topFH <- names(topHER2Features)[1:500]

# topFH saved to Synapse Entity #156229

## LOAD IN PREVIOUSLY SAVED METABRIC WORKSPACE
metaEnt <- loadEntity(138635)
load(file.path(metaEnt$cacheDir, metaEnt$files))

## ALIGN THE EXPRESSION AND PHENOTYPIC DATA
sortExprsDat <- exprsDat[ , exprsDatInd]
sortClinDat <- clinDat[clinDatInd, ]

## CREATE TRAINING AND VALIDATION COHORTS
set.seed(022712)
randVec <- rbinom(997, size = 1, prob = 0.5)
trainExpress <- exprs(sortExprsDat)[ , randVec == 0]
validExpress <- exprs(sortExprsDat)[ , randVec == 1]
trainPhen <- sortClinDat@data[randVec == 0, ]
validPhe <- sortClinDat@data[randVec == 1, ]

## LET'S REDUCE THE FEATURESET TO THE 'TOP' HER2 FEATURES
trainExpress <- trainExpress[topFH, ]
validExpress <- validExpress[topFH, ]

## TRANSPOSE
tTrain <- t(trainExpress)
tValid <- t(validExpress)
# 
# ## HER2 BINARY INDICATOR
her2Assign <- ifelse(sortClinDat$HER2_SNP6_state == "GAIN", 1, 0)
trainHER2 <- her2Assign[randVec == 0]
validHER2 <- her2Assign[randVec == 1]


