rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)

exampledata <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(exampledata)
#data <- data %>% rename(negbez_t1=negbez_1) %>% rename(negbez_t2=negbez_2)  %>% rename(posbez_t1=posbez_1) %>% rename(posbez_t2=posbez_2)

# create binary
exampledata <- exampledata %>% 
  mutate(schoolbin_1 = ifelse(schoolhigh_1 %in% c(1,2), 1,
                       ifelse(schoolhigh_1 %in% c(3,4), 2, NA))) %>% 
  mutate(schoolbin_2 = ifelse(schoolhigh_2 %in% c(1,2), 1,
                       ifelse(schoolhigh_2 %in% c(3,4), 2, NA))) %>% rename(zy = zyg)
acevars <- c("negbez","kultkapjahre")
acevars1 <- c("negbez1","kultkapjahre1")
acevars2 <- c("negbez2","kultkapjahre2")
acevarswide <- c("negbez_1","kultkapjahre_1","negbez_2","kultkapjahre_2")
covvarsall <- c("age","sex","posbez_1","posbez_2")
covvarswide <- c("posbez_1","posbez_2")
covvarslong_checked <- c("age","sex")
variables <- c(acevarswide,covvarsall)

sep = "_"
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvarsall) # Control variables 
l <- 3*nv*2
t <- m+l+c
# hiernach alle mxMatrix-Objekte checken, ob c verwendet wird! Diese müssen blockiert werden

pathCov_label_variance <- function(string) {
stringend <- substring(string, nchar(string)) == "1"
  if  (stringend == TRUE) {
c(paste0("b",string,1:nv),rep(NA,nv))
  }
else {
  c(rep(NA,nv),paste0("b",string,1:nv))
}
}

# label matrix of effects of wide formatted covariates --> has to be transposed if covariates used in mean matrix!!!
pathCovlabelvariance <- as.matrix(sapply(covvarswide,pathCov_label_variance))
if (length(pathCovlabelvariance)==0) {
  pathCovlabelvariance <- NULL
}
pathCovvaluevariance <- pathCovlabelvariance
pathCovvaluevariance[!is.na(pathCovvaluevariance)] <- .3
pathCovvaluevariance[is.na(pathCovvaluevariance)] <- 0
mode(pathCovvaluevariance) <- "numeric"
pathCovfreevariance <- pathCovvaluevariance==.3

## Cov Vars without variances
pathCov_label_constant <- function(string) {
paste0("b",rep(paste0(string,1:nv),2),rep(c(1,2),each=nv))
}
pathCovlabelconstant <- as.matrix(sapply(covvarslong_checked,pathCov_label_constant))
if (length(pathCovlabelconstant)==0) {
  pathCovlabelconstant <- NULL
}
pathCovvalueconstant <- pathCovlabelconstant
pathCovvalueconstant[!is.na(pathCovvalueconstant)] <- .3
pathCovvalueconstant[is.na(pathCovvalueconstant)] <- 0
mode(pathCovvalueconstant) <- "numeric"
pathCovfreeconstant <- pathCovvalueconstant==.3
pathCovlabel <- cbind(pathCovlabelconstant,pathCovlabelvariance)
pathCovvalue <- cbind(pathCovvalueconstant,pathCovvaluevariance)
pathCovfree <- cbind(pathCovfreeconstant,pathCovfreevariance)

if (!is.null(covvars)) {
pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
}
###############################################################################
# ab hier neuer code
###############################################################################

# before: 
    # covariance vector not part of the variables vector
    # reduced number of variables
    # no effects and no covariances of covariates -> matrix construction without them
    # define definition variables



# Mean Matrix
if (covariance == TRUE) {
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(variables, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = t, ncol = 1, 
                 free = c(rep(TRUE,m+c),rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M")
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")
}
if (covariance == FALSE) {
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(acevarswide, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = m+l, ncol = 1, 
                 free = c(rep(TRUE,m),rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M") 
# matrix with effect sizes of moderation of the means stored in M
pathCovMean <- mxMatrix(type = "Full", nrow = c, ncol = ntv, byrow = FALSE,
                            free = t(pathCovfree),
                            values = t(pathCovvalue),
                            labels = t(pathCovlabel),
                            name = "pCovM")

# Matrix of definition variables for mean moderation
labeldef <- paste0("data",".",covvarsall)
defM      <- mxMatrix( type="Full", nrow=1, ncol=c, free=FALSE, labels=labeldef, name="defM" )

# Matrix effects on means
effMean <- mxAlgebra(expression = defM%*%pCovM, name = "effM")

# Vector of latent variables (all set to zero) just need to concatenate them to the manifests to get the dimensions right
latentmeans <- mxMatrix(type = "Full", nrow = l, ncol = 1, name = "lmeans")

# Concatenate the manifest with the latent means vector
effMeanFull <- mxAlgebra(expression = cbind(effM,lmeans), name = "effMFull")

# Matrix of moderated means
modMean <- mxAlgebra(expression =M+effMFull, name = "modM")

# Matrix of expected means
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%modM), name = "expMean")
}



