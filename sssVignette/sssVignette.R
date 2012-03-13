# sssVignette.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# This is a vignette for learning how to use our R package implementation of
# 'Shotgun Stochastic Search' by Chris Hans, Quanli Wang, and Mike West. The 
# package is authored by Brian Bot of Sage Bionetworks and conveniently wraps
# and integrates SSS with R

require(Biobase)
require(factDesign) # A Bioconductor package
require(sss)
require(hgu95av2.db)

## Bring in the 'Estrogen' dataset from the 'factDesign' package
data(estrogen)
exDesign <- pData(estrogen)

# ES TIME
# et1.CEL  A  10h
# et2.CEL  A  10h
# Et1.CEL  P  10h
# Et2.CEL  P  10h
# eT1.CEL  A  48h
# eT2.CEL  A  48h
# ET1.CEL  P  48h
# ET2.CEL  P  48h

# This sample dataset comprises estrogen receptor-positive cell lines that are 
# exposed to estrogen (versus control) over two timepoints. It is a restricted 
# 500 probe features pulled from Affymetrix HGU95av2 arrays for the above
# samples

## Let's create a simple SSS logistic regression model for estrogen exposure
exEstrogen <- exprs(estrogen)

## Let's build a model on the 10h samples and validate on the 48h samples
validSet <- t(exEstrogen[ , 1:4]) # Note that the matrix is transposed
trainSet <- t(exEstrogen[ , 5:8])
validPheno <- ifelse(exDesign$ES[1:4] == "P", 1, 0)
trainPheno <- ifelse(exDesign$ES[5:8] == "P", 1, 0)

# Fit the model
sssFit <- sss(trainPheno ~ trainSet, iters = 5000, nbest = 200)

# Run the training set back through the model object
trainPhenoHat <- predict(sssFit, newdata = trainSet)

# Look at the correlation of 'trainPheno' and 'trainPhenoHat'
trainRho <- cor(trainPheno, trainPhenoHat)
# > trainRho^2 = 0.9928368

# Look at the validation set
validPhenoHat <- predict(sssFit, newdata = validSet)

# Look at the correlation of 'validPheno' and 'validPhenoHat'
validRho <- cor(validPheno, validPhenoHat)
# > validRho^2 = 0.9664588

## Let's look at the features and their posterior marginal probabilities
sPMP <- sort(sssFit@postMargProb, decreasing = TRUE)

## MAPPING ILLUMINA FEATURES TO GENE SYMBOLS
topEstrogenSYMs <- as.character(mget(names(sPMP[1:10]), hgu95av2SYMBOL, 
                              ifnotfound = NA))