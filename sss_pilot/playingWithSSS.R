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
# trainExpress <- exprs(sortExprsDat)[ , randVec == 0]
# validExpress <- exprs(sortExprsDat)[ , randVec == 1]
# trainPhen <- sortClinDat@data[randVec == 0, ]
# validPhe <- sortClinDat@data[randVec == 1, ]
# 
# ## LET'S REDUCE THE FEATURESET TO THE 'TOP' HER2 FEATURES
# trainExpress <- trainExpress[topFH, ]
# validExpress <- validExpress[topFH, ]
# 
# ## TRANSPOSE
# tTrain <- t(trainExpress)
# tValid <- t(validExpress)
# 
# ## HER2 BINARY INDICATOR
her2Assign <- ifelse(sortClinDat$HER2_SNP6_state == "GAIN", 1, 0)
# trainHER2 <- her2Assign[randVec == 0]
# validHER2 <- her2Assign[randVec == 1]

## RUN SSS ON THE FULL DATA MATRIX WITH RANDVEC SPECIFYING TRAINING AND VALIDATION
exprMod <- t(exprs(sortExprsDat)[topFH, ])
sssFit <- sss(her2Assign ~ t(exprs(sortExprsDat)[topFH, ]), weights = randVec)

pm <- unlist(sssFit@score)
pm <- exp(pm - max(pm))
pm <- pm/sum(pm)
k <- max(unlist(sssFit@p))
pmax <- length(unlist(sssFit@score))

itrain <- randVec == 0
ivalid <- randVec == 1
trainMean <- apply(exprs(sortExprsDat)[ , itrain], 1, mean)
trainSD <- sqrt(apply(exprs(sortExprsDat)[ , itrain], 1, var))

pmk <- pm[1:k]/sum(pm[1:k]);                        # condition on only the top k models chosen

train.m <- apply(x[,itrain],1,mean)                 # must only standardize based on the observations
train.sd <- sqrt(apply(x[,itrain],1,var))           # used to fit the models
X <- (t(x) - matrix(train.m,ncol=N,nrow=n,byrow=T)) # subtract the mean
X <- X/matrix(train.sd,ncol=N,nrow=n,byrow=T)       # divide the the sd                    

Fit <- matrix(0,nrow=n,ncol=k);                     # to save fits and predictions
pFit <- matrix(0,nrow=n,ncol=k);                     
pmp <- matrix(0,nrow=N,ncol=1);                     # to save marginal inclusive probs

for(j in 1:k)
{
  p <- models[j,1]                              # dim of this model
  if (p>0)
  {
    ig <- models[j,2+(1:p)]                   # predictors in model    
    pmp[ig] <- pmp[ig] + pmk[j];           # posterior probs on predictors 
    
    b <- matrix(as.numeric(models[j,(p+3):(2*p+3)]),ncol=1,nrow=p+1)
    # post mode of regn parameters
    
    A <- matrix(0,nrow=n,ncol=(p+1)); A[,1] <- 1;
    A[,-1] <- X[,ig]                          # design matrix
  }
  if(p==0)
  {
    b <- models[j,3]; A <- matrix(1,ncol=n,nrow=1); 
  }   
  Fit[,j] <- A%*%b;                       # fitted & predicted linear regn
  pFit[,j] <- 1 / (1+exp(-Fit[,j]));      #    ... and plug-in probabilities
  
  # predictions for hold-out validation cases are already in the above
  plot(1:n, pFit[,j], xlab="case id",ylab="Fitted Probabilities",ylim=c(0,1),type="n",main=paste("Model",j,": p =",p," : Prob =",round(pmk[j],4)),xlim=c(1,n + n*0.2),axes=F)
  box(); axis(2); tt <- axTicks(1); tt <- tt[tt<=n]; axis(1,at=tt)
  legend("bottomright", c("obs 0","obs 1","pred 0","pred 1","base"),pch=c(16,16,16,16,-1),lty=c(-1,-1,-1,-1,2),col=c("blue","red","cyan","magenta","black"),merge=T,bg="gray90")
  
  showtv(t(pFit), y, itrain, ivalid, j);
  abline(h=sum(y)/n,lty=2);                         
  par(ask=T)
}

avepFit <- pFit%*%pmk
plot(1:n, pFit[,j], xlab="case id",ylab="", ylim=c(0,1),type="n",main=paste("Model Average"),xlim=c(1,n + n*0.2),axes=F)
box(); axis(2); tt <- axTicks(1); tt <- tt[tt<=n]; axis(1,at=tt)
legend("bottomright", c("obs 0","obs 1","pred 0","pred 1","base"),pch=c(16,16,16,16,-1),lty=c(-1,-1,-1,-1,2),col=c("blue","red","cyan","magenta","black"),merge=T,bg="gray90")
showtv(t(avepFit),y,itrain,ivalid,1)
abline(h=sum(y)/n,lty=2)
