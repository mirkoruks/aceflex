# brain storming for dynamic start values of path coefficients (-> less problems when there are mayor differences in variances because then a starting value of e.g. 0.3 is not appropriate!)

rm(list = ls())
library(OpenMx)
library(dplyr)
#library(haven)
library(foreign)
library(psych)
# simulate the data
set.seed(1)
cov1long <- rnorm(n = 1000, mean = 20, sd = 16.2)
cov2long <- rnorm(n = 1000, mean = 10, sd = 9.7)
cov3wide_1 <- rnorm(n = 1000, mean = 7.8, sd = 3.8) 
cov3wide_2 <- rnorm(n = 1000, mean = 7.8, sd = 3.8)
cov4wide_1 <- rnorm(n = 1000, mean = 5.6, sd = 1.9)
cov4wide_2 <- rnorm(n = 1000, mean = 5.6, sd = 1.9)
ey_1 <- rnorm(1000, mean = 0, sd = 2)
ey_2 <- rnorm(1000, mean = 0, sd = 2)
ez_1 <- rnorm(1000, mean = 0, sd = 3)
ez_2 <- rnorm(1000, mean = 0, sd = 3)
y_1 <- cov1long*.7 + cov2long*0.8 + cov3wide_1*1.4 + cov4wide_1*2  + ey_1
y_2 <- cov1long*.7 + cov2long*0.8  + cov3wide_2*1.4 + cov4wide_2*2 + ey_2
z_1 <- cov1long*0.01 + cov2long*2 + cov3wide_1*0.9 + cov4wide_1*3  + ez_1
z_2 <- cov1long*0.01 + cov2long*2  + cov3wide_2*0.9 + cov4wide_2*3 + ez_2
usedata <- data.frame(cbind(y_1,z_1,y_2,z_2,cov1long,cov2long,cov3wide_1,cov4wide_1,cov3wide_2,cov4wide_2))
summary(usedata)


# preliminaries
covvarslong_checked <- c("cov1long","cov2long")
covvarswide_checked <- c("cov3wide_1","cov4wide_1","cov3wide_2","cov4wide_2")
covvarswide_checked <- c("cov3wide_1","cov3wide_2")
covvarsall <- c(covvarslong_checked,covvarswide_checked)
acevars1 <- c("y_1","z_1")
acevars <- c("y","z")
#acevars1 <- c("y_1")
#acevars <- c("y")
nv <- length(acevars)
covvars1 <- grep("_1$", covvarswide_checked, value = TRUE)
covvars2 <- grep("_2$", covvarswide_checked, value = TRUE)
covvarswide <- c(covvars1, covvars2)

################################################################################
################################################################################

Ycov <- as.matrix(subset(usedata, select = acevars1))
Xcov <- cbind(1,as.matrix(subset(usedata, select = c(covvarslong_checked,covvars1))))
pathCovstart <- (t(solve(t(Xcov)%*%Xcov)%*%t(Xcov)%*%Ycov))[,2:ncol(Xcov)]
if (is.null(dim(pathCovstart))) {
pathCovstart <- matrix(t(pathCovstart), ncol = length(pathCovstart))
} 
pathCovstart <- pathCovstart[rep(seq_len(nrow(pathCovstart)), 2), ] # complete matrix 

rownames(pathCovstart) <- NULL
colnames(pathCovstart) <- NULL
pathCovstartlong <- pathCovstart[,1:length(covvarslong_checked)]
pathCovstartwide <- pathCovstart[,(length(covvarslong_checked)+1):ncol(pathCovstart)]
if (is.null(ncol(pathCovstartwide))) {
pathCovstartwide <- matrix(pathCovstartwide, nrow = length(pathCovstartwide))
}
pathCovstartwide <- pathCovstartwide[, rep(1:ncol(pathCovstartwide), each=2)]


# Some helper functions for the labeling process etc.
pathCov_label_variance <- function(string) {
stringend <- substring(string, nchar(string)) == "1"
  if  (stringend == TRUE) {
c(paste0("b",string,1:nv),rep(NA,nv))
  }
else {
  c(rep(NA,nv),paste0("b",string,1:nv))
}
}

pathCovlabelvariance <- as.matrix(sapply(covvarswide,pathCov_label_variance))
if (length(pathCovlabelvariance)==0) {
  pathCovlabelvariance <- NULL
}
pathCovvaluevariance <- pathCovlabelvariance
pathCovvaluevariance[!is.na(pathCovvaluevariance)] <- .3
pathCovvaluevariance[is.na(pathCovvaluevariance)] <- 0
mode(pathCovvaluevariance) <- "numeric"
pathCovfreevariance <- pathCovvaluevariance==.3
pathCovvaluevariance <- (pathCovvaluevariance*10/3)*pathCovstartwide

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
pathCovvalueconstant <- (pathCovvalueconstant*10/3)*pathCovstartlong

pathCovlabel <- cbind(pathCovlabelconstant,pathCovlabelvariance)
pathCovvalue <- cbind(pathCovvalueconstant,pathCovvaluevariance)
pathCovfree <- cbind(pathCovfreeconstant,pathCovfreevariance)


################################################################################



# Calculate start values of effects of latent ACE factors on acevars
    # assumption: correlation = 0.3 for univariate effects for all ACE factors 
        # new argument: coruniv = c(0.3,0.3,0.3)
            # explication: assumed correlation of unique A factor with manifest = 0.3, of unique C factor = 0.3 of unique E factor = 0.3, the order is A, C, E
    # assumption: correlation = 0.2 for bivariate effects for all ACE factors
        # new argument: corbiv = c(0.2,0.2,0.2) 
            # explication: same as above

# formula of bivariate beta (regression): beta = cov(x,y)/var(x)    
    # we know that cor(x,y)=cov(x,y)/sd(x)*sd(y) and we assume that cor(x,y) = 0.3/0.2 while we constrain the variance of the latent factors always to 1 -> var(x)=sd(x)=1
    # so we get 0.3=cov(x,y)/sd(y) <=> cov(x,y)=0.3*sd(y) for the univariate effects
    # and 0.2=cov(x,y)/sd(y) <=> cov(x,y)=0.2*sd(y) for the bivariate effects
    # substituting the terms in the formula for the beta, leads us to:
        # beta = cov(x,y)/var(x) <=> beta = 0.3*sd(y) for the univariate effects
        # beta = cov(x,y)/var(x) <=> beta = 0.2*sd(y) for the bivariate effects

# additional argument: coracerandom = FALSE -> if TRUE the values provided in coruniv and corbiv are used as start values for a random draw of assumed correlations
    # coracerandomtrials = 10
    # coracerandomdist = "normal" # alternative values = "uniform" or "cautchy" or "log" -> values provided in coruniv and corbiv are used as mean of distribution as well, sd of distribution is set to 1

describe(usedata)
varacevars <- diag(var(na.omit(usedata)[,acevars1]))

mat <- matrix(varacevars,nrow = nv,ncol = nv, byrow = FALSE) # first row = Var of first acevar; second row = Var of second acevar, ...
mat[upper.tri(mat, diag = FALSE)] <- 0
mat <- sqrt(mat)
corunivA <- coruniv[1]
corbivA <- corbiv[1]
valuesA <- mat
diag(valuesA) <-diag(valuesA)*corunivA
valuesA[lower.tri(valuesA, diag = FALSE)] <- valuesA[lower.tri(valuesA, diag = FALSE)]*corbivA
corunivC <- coruniv[2]
corbivC <- corbiv[2]
valuesC <- mat
diag(valuesC) <-diag(valuesC)*corunivC
valuesC[lower.tri(valuesC, diag = FALSE)] <- valuesC[lower.tri(valuesC, diag = FALSE)]*corbivC
corunivE <- coruniv[2]
corbivE <- corbiv[2]
valuesE <- mat
diag(valuesE) <-diag(valuesE)*corunivE
valuesE[lower.tri(valuesE, diag = FALSE)] <- valuesE[lower.tri(valuesE, diag = FALSE)]*corbivE
valuesEbeta <- valuesE
valuesEbeta[lower.tri(valuesEbeta, diag = FALSE)] <- 0
