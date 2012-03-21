# bfrmTesting.R

require(synapseClient)
require(affy)
require(bfrm)

## BRING IN THE EGFR EXPERIMENT AND CONTROL CEL FILES FROM SYNAPSE
egfrRawEnt <- loadEntity(138511)
egfrBatch <- ReadAffy(filenames=file.path(egfrRawEnt$cacheDir, egfrRawEnt$files))

## NORMALIZE AND SUMMARIZE USING RMA
egfrNorm <- rma(egfrBatch, normalize = T, background = F)
egfrMat <- exprs(egfrNorm)

## IDENTIFY EXPERIMENT VS CONTROL SAMPLES
findControls <- grep("GFP", colnames(egfrMat), ignore.case = F)
experimentInd <- vector()
experimentInd[1:9] <- 1
experimentInd[findControls] <- 0

## DEFINE NORMALIZTION CONTROL PROBESETS
findAFFX <- grep("AFFX", rownames(egfrMat))
covMat <- egfrMat[findAFFX, ]
covComp <- svd(covMat)
batchAdj <- t(covComp$v[ , 1:5])

## RUN BFRM IN SPARSE ANOVA MODE
egfrSparseANOVA <- bfrm(egfrMat, design = experimentInd, 
                        control = batchAdj, burnin = 100, nmcsamples = 500)


