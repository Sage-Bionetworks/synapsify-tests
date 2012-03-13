# .########......######.....###....########...######..##.....##.##.......########
# .##.....##....##....##...##.##...##.....##.##....##.##.....##.##.......##......
# .##.....##....##........##...##..##.....##.##.......##.....##.##.......##......
# .########.....##.......##.....##.########...######..##.....##.##.......######..
# .##...##......##.......#########.##..............##.##.....##.##.......##......
# .##....##.....##....##.##.....##.##........##....##.##.....##.##.......##......
# .##.....##.....######..##.....##.##.........######...#######..########.########


# playingWithSSS.R

require(randomForest)
require(sss)
require(illuminaHumanv4.db)
require(synapseClient)
require(Biobase)
require(ROCR)
require(ggplot2)

synapseLogin()

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
her2GSym <- as.character(mget(names(topHER2Features), illuminaHumanv4SYMBOL, 
                              ifnotfound = NA))

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

## LET'S MODEL WITH JUST TRAINING AND USE THE PREDICT METHOD
## ON A HELD OUT VALIDATION SET

sssFit <- sss(trainHER2 ~ tTrain)
validVec <- predict(sssFit, newdata = tValid)

## EVALUATE MODEL PERFORMANCE
her2Pred <- prediction(validVec, validHER2)
her2Perf <- performance(her2Pred, "tpr", "fpr")
her2AUC <- performance(her2Pred, "auc")

## FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
her2SSPerf <- performance(her2Pred, "sens", "spec")
youdensJ <- her2SSPerf@x.values[[1]] + her2SSPerf@y.values[[1]] - 1
jMax <- which.max(youdensJ)
optCut <- her2Perf@alpha.values[[1]][jMax]

optSens <- unlist(her2SSPerf@x.values)[jMax]
optSpec <- unlist(her2SSPerf@y.values)[jMax]

# Sensitivity at Youden's J-point = 0.89
# Specificity at Youden's J-point = 0.85

## CREATE A BOXPLOT USING GGPLOT
dfHER2 <- cbind(validVec, validHER2)
colnames(dfHER2) <- c("predictions", "trueStat")
dfHER2 <- as.data.frame(dfHER2)

her2Box <- ggplot(dfHER2, aes(as.factor(trueStat), predictions)) + 
  geom_boxplot() + geom_jitter(aes(colour = as.factor(trueStat)))

her2Box <- her2Box + geom_hline(yintercept = optCut, colour = "red",
                                linetype = 2)
                                  

rankSum <- wilcox.test(validVec, validHER2)
# p-value < 2.2e-16

## CREATE A ROC CURVE USING GGPLOT
dfPerf <- as.data.frame(cbind(unlist(her2Perf@x.values),
                              unlist(her2Perf@y.values)))
colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")

her2ROC <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
  geom_line()
her2ROC <- her2ROC + geom_abline(slope = 1, colour = "red")

# AUROC = 0.944

# ## LET'S PUSH IN THE FULL MATRIX AND SEE WHAT WE GET
# sssFullFit <- sss(her2Assign ~ t(exprs(sortExprsDat)[topFH, ]), 
# training = randVec)


