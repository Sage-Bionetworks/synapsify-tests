# bfrmTesting.R

require(synapseClient)
require(affy)
require(bfrm)
require(ggplot2)

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
batchAdj <- t(covComp$v[ , 1:9])

## RUN BFRM IN SPARSE ANOVA MODE
egfrSparseANOVA <- bfrm(egfrMat, design = experimentInd, 
                        control = batchAdj, burnin = 2000, nmcsamples = 5000)
# qplot(egfrSparseANOVA@results$mPostPib[ , 2], geom = "histogram")

egfrSparseANOVA2 <- bfrm(egfrMat, design = experimentInd,
                         burnin = 2000, nmcsamples = 5000)
# qplot(egfrSparseANOVA2@results$mPostPib[ , 2], geom = "histogram")

postProbDF <- as.data.frame(t(rbind(egfrSparseANOVA2@results$mPostPib[ , 2], 
                                    egfrSparseANOVA@results$mPostPib[ , 2])))
colnames(postProbDF) <- c("withoutAdjustment", "withAdjustment")

## VISUALIZE IN GGPLOT AND VISUALLY ASSESS THE DIFFERENCE BETWEEN BATCH 
## ADJUSTED AND UNADJUSTED
meltPPDF <- melt(postProbDF) # For ggplot2
ggplot(meltPPDF, aes(value, fill = factor(variable))) + 
  geom_density(alpha = 0.3)

## GENERATE INDICES FOR THE "TOP" TRANSCRIPTS (POST. PROB ≥ 0.99)
topTxLogical <- egfrSparseANOVA@results$mPostPib[ , 2] >= 0.99
topTxInd <- grep("TRUE", topTxLogical)

## RUN BFRM IN THE EVOLVING MODE WITH THE TOPTXIND
egfrSparseFactors <- bfrm(egfrMat,
                          design = experimentInd,
                          control = batchAdj,
                          burnin = 100,
                          nmcsamples = 500,
                          evol = 1,
                          evolvarin = "topTxInd.txt",
                          evolmaximumfactors = 50,
                          evolmaximumvariables = 5091,
                          evolmaximumvariablesperiteration = 30)



