rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)

ordinal <- c("schoolhigh","ado")
sep = "_"
data <- data_orig
existence <- function(variable) {
  result <- NULL
if(variable %in% colnames(data))
{
  result <- NULL
}
else if (!(variable %in% colnames(data))) {
  result <- variable
}
}

existenceerror <- function(result) {
  if (!is.null(result)) {
  stop(c("I could not find a variable(s) in the data frame corresponding to the following variable strings you gave me: ",paste(result, sep = " ", collapse = ", ")))
  
  }
}

#############################################################################
# ordinal part
#############################################################################

if (!is.null(ordinal)) {
# Check for categorical variables in acevars 
if (!is.null(ordinal) & class(ordinal)== "character") {
ordinal1 <-    paste0(ordinal,sep,"1") # Covariates twin 1
ordinal2 <-    paste0(ordinal,sep,"2") # Covariates twin 2
}
ordinalwide <- c(ordinal1,ordinal2)
existence_check_ordvars <- unlist(lapply(ordinalwide, existence))
existenceerror(existence_check_ordvars)
}

ordinalwide
#
#
#
#
usedata <- data

# check number of levels of ordinal variable (crucial question is: ordinal or binary?) -> necessary for later conditional statements 
nlevels <- function(variable) { 
nolevels <- sort(unique(usedata[,variable]))
}

llevels <- function(variable) { 
length <- length(sort(unique(usedata[,variable])))
result <- NULL
if(length>1)
{
  result <- NULL
}
else if (length == 1) {
  result <- variable
}

}

levelserror <- function(result) {
  if (!is.null(result)) {
  stop(c("Check your ordinal variables.. The following variables have only one category: ",paste(result, sep = " ", collapse = ", ")))
  
  }
  }
levelslist <- lapply(ordinalwide, nlevels)
levelserror(unlist(lapply(ordinalwide, llevels)))

usedata[,ordinalwide] <- mxFactor(usedata[,ordinalwide], levels = levelslist)
print(summary(usedata))
ordinallength <- unlist(lapply(levelslist,length))
if (2 %in% ordinallength) {
    stop("Sorry, at the moment the function does not support binary variables")
}
# Define objects for threshold matrix
nTh       <- ordinallength-1 # No of thresholds
print(nTh)
ntvo      <- length(ordinalwide) # Total No of ordinal vars

freeThresholds <- function(nt) { # assumes mxMatrix(byrow = FALSE)
values <- NULL
    if(nt == 1) {
    values <- c(TRUE,rep(FALSE,(max(nTh)-1))) 
    }
    else if (nt > 1 & nt < max(nTh)) {
       values <- c(rep(FALSE,2),rep(TRUE,(nt-2)),rep(FALSE,(max(nTh)-nt)))        
    }
    else if (nt > 1 & nt == max(nTh)) {
       values <- c(rep(FALSE,2),rep(TRUE,(nt-2))) 
    }
}
valThresholds <- function(nt) { # assumes mxMatrix(byrow = FALSE)
values <- NULL
    if(nt == 1) {
    values <- c(.5,rep(0,(max(nTh)-1))) # unused thresholds set to zero (check: post #20 by T. Bates in: https://openmx.ssri.psu.edu/node/4538)
    }
    else if (nt > 1 & nt < max(nTh)) {
      values <- c(0,rep(1,(nt-1)),rep(0,(max(nTh)-nt)))           
    }
    else if (nt > 1 & nt == max(nTh)) {
       values <- c(0,rep(1,(nt-1))) 
    }
}

lbThresholds <- function(nt) { # assumes mxMatrix(byrow = FALSE)
lb <- NULL
    if(nt == 1) {
    lb <- rep(NA,max(nTh))
    }
    else if (nt > 1 & nt < max(nTh)) {
      lb <- c(-3,rep(.0001,(nt-1)),rep(NA,(max(nTh)-nt)))         
    }
    else if (nt > 1 & nt == max(nTh)) {
       lb <- c(-3,rep(.0001,(nt-1)))  
    }
}

labelThresholds <- function(vars) { # assumes mxMatrix(byrow = FALSE)
labels <- NULL
noTh <- (length(levels(usedata[,vars])))-1
varsnew <- unlist(sapply(strsplit(vars, split=sep, fixed=TRUE), function(x) (x[1])))
    if (noTh == 1) { # binary
        labels <- c(paste0("th",varsnew,1),rep(NA,(max(nTh)-1)))
    }
    else if (noTh > 1 & noTh < max(nTh)) { # ordinal
        labels <- paste0("th",varsnew,(1:noTh),rep(NA,(max(nTh)-noTh)))
    }
    else if (noTh > 1 & noTh == max(nTh)) { # ordinal
        labels <- paste0("th",varsnew,(1:noTh))
    }
}

frTh <- unlist(lapply(nTh, freeThresholds))
svTh <- unlist(lapply(nTh, valThresholds))
lbTh <- unlist(lapply(nTh, lbThresholds))
labTh <- unlist(lapply(ordinalwide, labelThresholds))

# Thresholds definieren
thinG     <- mxMatrix(type="Full", nrow=max(nTh), ncol=ntvo, free=frTh, byrow = FALSE, values=svTh, lbound=lbTh, labels=labTh, name="thinG") # matrix of threshold increments
inc       <- mxMatrix(type="Lower", nrow=max(nTh), ncol=max(nTh), free=FALSE, values=1, name="inc") # matrix of lower 1
# Example how the premultiplication with a matrix of lower 1 ensures the ordering of the thresholds
#mat1 <- matrix(svTh, nrow=max(nTh), ncol=ntvo)
#lower <- matrix(1,nrow=max(nTh), ncol=max(nTh))
#lower[upper.tri(lower)] <- 0
threG     <- mxAlgebra(expression= inc %*% thinG, name="threG") # Multiplikation stellt sicher, dass Th3>Th2>Th1, wobei Th2>Th1 schon durch die Fixierung auf 0, bzw. 1 gesichert ist!
#umxThresholdMatrix(df = data, selDVs = ordinalwide, sep = "_", method = "Mehta")



# flag binary variables
acevars <- c("ado_1","ado_2")
covvars <- NULL
variables <- c(acevars,covvars)
binarytrue <- nTh == 1
binaryvar <- ordinalwide[binarytrue]
unlist(sapply(strsplit(binaryvar, split=sep, fixed=TRUE), function(x) (x[1])))

checkcorrespondence <- function(check,comparison) {
    checkresult <- NULL
    checkresult <- check %in% comparison
    checkresult <- TRUE %in% checkresult
}

if (2 %in% ordinallength) {
binaryace <- unlist(lapply(acevars,checkcorrespondence, check = binaryvar))
binaryrest <- rep(FALSE,(length(variables)-length(acevars)))
binaryrest
binaryflag <- c(binaryace,binaryrest)
binaryflag
# fixed means
binaryflag # FALSE = Non binary; TRUE = binary -> for the mean estimation we can reverse it 
 # if there are any binary variables
svmean_manifests <- binaryflag == FALSE
svmean_manifests # insert this vector into the "free"-argument if there are binary variables in the model
}


#if (2 %in% ordinallength) {
# fixed variances
# 2 are binary (1st and 3rd)  
  # No of rows = no of binary vars
nrowfilterbinary <- sum(binaryflag)
  # No of cols = vars in total
ncolsfilterbinary <- length(variables)
ncolsfilterbinary
  # 1 if var = binary and 0 if not
# function: while row
bfilter <- function(x,vec) {
  result <- list()
  vec[-x] <- 0
  result <- vec
  }

valfilterbinary <- binaryflag
valfilterbinary[valfilterbinary== TRUE] <- 1
flag <- which(valfilterbinary == 1)
filtermatvalues <- matrix(unlist(lapply(flag, bfilter, vec = valfilterbinary)),nrow = nrowfilterbinary, ncol = length(valfilterbinary), byrow = TRUE)

filtermatbin <- mxMatrix(type = "Full", values = filtermatvalues, name = "fmatbin")
binarycov <- mxAlgebra(expression = fmatbin %*%expCovMZ %*% t(fmatbin), name = "binCov")
one <- mxMatrix(type = "Unit", nrow = length(valfilterbinary), ncol = 1, name = "Unit")
var1 <- mxConstraint(expression = diag2vec(binCov)==Unit , name = "VConstraint1")

#}
# here are the variances of the binaries that we have to constrain to zero! so we have to apply the mxConstraint to this matrix!
diag(onlybinary) 

rsum.cumdiff <- function(x, n = 3L) (cs <- cumsum(x))[-(1:(n-1))] - c(0,cs[1:(length(x)-n)])
cumsum(valfilterbinary)
