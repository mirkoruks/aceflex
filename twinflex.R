twinflex <- function(acevars, data, zyg, sep, covvars=NULL, covariance = TRUE, ordinal = NULL, optimizer = NULL, tryHard = FALSE, type = "chol") {

if ("OpenMx" %in% (.packages()) == FALSE) {
  stop("You need to load the OpenMx library")
}
  
if (is.null(acevars)) {
  stop("You need to specify a string vector with variables for the ACE decomposition")
}
if (is.null(data)) {
  stop("You need to specify a data set")
}
  
if (is.null(zyg)) {
  stop("You need to specify a zygosity variable")
}
  
if (is.null(sep)) {
  stop("You need to specify a separator. Example: The variable names are 'iq_t1' and 'iq_t2', then your seperator is '_t'")
}
  
if (length(sep) > 1) {
  stop("Please specify the separator with a string vector with one element only.")
}
  
if (!is.null(covvars) & class(covvars)!= "character") {
  stop("Please, use a string vector to specify the covariates")
}

if (!is.null(ordinal)) {
if (!class(ordinal)== "character") {
  stop("The ordinal-argument needs a character vector")
}
ocheck <- ordinal %in% acevars
if (FALSE %in% ocheck) {
  stop("All the ordinal variables need to be added in the acevars argument as well")
}
}
# Output-List -> Print results
output <- list()   


# Create acevars in wide format using sep
acevars1 <-    paste0(acevars,sep,"1") # ACE vars twin 1
acevars2 <-    paste0(acevars,sep,"2") # ACE vars twin 2

############ CHECK IF ACEVARS AND/OR COVVARS 1) EXIST AND 2) FOR THE WIDE FORMATTED ONES: IF THEY HAVE WITHIN-VARIANCE

#######
  # 1. check if variables exist --> for acevars in wide format and for covvars in both long (without within-pair variance) and wide (with within-pair variance)
#######

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


  # 1a: Check for acevars (only in wide since they have to have within-pair-variance)
acevars1 <-    paste0(acevars,sep,"1") # Covariates twin 1
acevars2 <-    paste0(acevars,sep,"2") # Covariates twin 2
acevarswide <- c(acevars1, acevars2)
existence_check_acevars <- unlist(lapply(acevarswide, existence))
if (!is.null(existence_check_acevars)) {
acevars_not_found <- unlist(sapply(strsplit(existence_check_acevars, split=sep, fixed=TRUE), function(x) (x[1])))
acevars_not_found <- unique(acevars_not_found)
#existenceerror(existence_check_acevars)
} else {
  acevars_not_found <- NULL
}

  # 1b: Check for covvars (can be in wide and long) -> 1. Check for wide -> 2. Check for Long if something does not appear as wide
if (!is.null(covvars) & class(covvars)== "character") {
covvars1 <-    paste0(covvars,sep,"1") # Covariates twin 1
covvars2 <-    paste0(covvars,sep,"2") # Covariates twin 2
covvarswide <- c(covvars1, covvars2)

existence_check_covvars <- unlist(lapply(covvarswide, existence))
covvarswide_checked <- covvarswide[!covvarswide %in% existence_check_covvars] # object with wide-formatted covariates 
existence_check_covvars # object with possibly long-formatted covariates -> check if really long or if typo

if (!is.null(existence_check_covvars)) { # if the object with possibly long-formatted covariates is not empty -> check if there are long-formatted covariates
covvars_possibly_long_format <- covvarswide[covvarswide %in% existence_check_covvars] 
covvars_possibly_long_format <- unlist(sapply(strsplit(covvars_possibly_long_format, split=sep, fixed=TRUE), function(x) (x[1])))
covvars_possibly_long_format <- unique(covvars_possibly_long_format)
existence_check_covvars2 <- unlist(lapply(covvars_possibly_long_format, existence))
existenceerror(c(acevars_not_found,existence_check_covvars2))
covvarslong_checked <- covvars_possibly_long_format
} else {
covvarslong_checked <- NULL # if the condition is fulfilled, then all the covariates are in wide format and none in long
}
if (length(covvarswide_checked)==0) {
  covvarswide_checked <- NULL
}
covvarsall <- c(covvarslong_checked,covvarswide_checked) # from now on: if a covariate ends with _1 or _2 -> it is wide!
suf1 <- paste0(sep,"1")
suf2 <- paste0(sep,"2")

covvars1 <- covvarswide_checked[grepl(suf1, covvarswide_checked)]
covvars2 <- covvarswide_checked[grepl(suf2, covvarswide_checked)]
covvarswide <- c(covvars1, covvars2)
} else if (is.null(covvars)) {
covvarsall <- NULL
covvars1 <- NULL
covvars2 <- NULL
covvarswide <- NULL
covvarslong_checked <- NULL
}

varswide1 <- c(acevars1,covvars1)
varswide2 <- c(acevars2,covvars2)



cat("\n\n\nACE Variables: \n\n")
print(acevarswide)

if (!is.null(covvarsall)) {
if (!is.null(covvarswide)) {
  cat("\n\n\nThere are ",length(covvarswide)," wide-formatted covariates: \n\n")
print(covvarswide)
} 
if (is.null(covvarswide)) {
    cat("\n\n\nThere are not wide-formatted covariates\n\n")
} 
if (!is.null(covvarslong_checked)) {
  cat("\n\n\nThere are ",length(covvarslong_checked)," long-formatted covariates: \n\n")
print(covvarslong_checked)
} 
if (is.null(covvarslong_checked)) {
    cat("\n\n\nThere are not long-formatted covariates\n\n")
}
cat("\n\n\nThere are ",length(covvarsall),"covariates in total: \n\n")
print(covvarsall)
} else {
cat("\n\n\nThere are no covariates at all \n\n")
}

varswideonly <- c(acevarswide,covvarswide)
cat("\n\n\nThere are ", length(varswideonly), "wide-formatted Variables: \n\n")
print(varswideonly)

variables <- c(acevarswide,covvarsall)
cat(paste0("\n\n\nAll in all, there are ",length(variables)," variables: \n\n"))
print(variables)

rna <- function(x) replace(x, is.na(x), "")
checkvariance <- function(v1,v2) {
identicalcheck <- as.vector(colSums(ifelse(rna(data[,v1, drop = FALSE])==rna(data[,v2, drop = FALSE]), 0, 1)))
if (0 %in% identicalcheck == TRUE) {
ind <- identicalcheck==0
result1 <- v1[ind]
result2 <- v2[ind]
result <- c(result1,result2)
if (!is.null(result)) {
  stop(c("The following acevars have no within-variance: ",paste(result, sep = " ", collapse = ", ")))
}
}
else {
  cat("\n\n\nAll the wide-formatted variables have within-pair variance")
}
}
checkvariance(varswide1,varswide2)

existence_check_zyg <- existence(zyg)
existenceerror(existence_check_zyg)

usevariables <- c(variables,zyg)
usedata <- subset(data, select = c(usevariables)) # Data set with only variables used for the model
cat("\n\n\nSummary total Data\n\n")
print(summary(usedata))

# Check if zygosity variable is coded correctly
if (min(usedata[,zyg]) != 1 & max(usedata[,zyg]) != 2) {
          stop("Zygosity variable must be coded as follows: 1 = MZ, 2 = DZ. Please, recode the zygosity variable.")
}
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
# flag ordinal variables in acevars
checkcorrespondence <- function(check,comparison) {
    checkresult <- NULL
    checkresult <- check %in% comparison
    checkresult <- TRUE %in% checkresult
}

# Starting Values
  # Means Vector
    # acevars
svmeanacevarswide1 <- colMeans(usedata[,acevars1, drop = FALSE], na.rm=TRUE)
svmeanacevarswide2 <- colMeans(usedata[,acevars2, drop = FALSE], na.rm=TRUE)
svmeanacevars <- rowMeans(cbind(svmeanacevarswide1, svmeanacevarswide2), na.rm=TRUE)
if (!is.null(ordinal)) {
  flagordinal <- unlist(lapply(acevars1,checkcorrespondence, check = ordinalwide))
  svmeanacevars[flagordinal] <- 0
}
cat("\n\n\nStarting Values of the acevars for the mean vector\n\n")
print(svmeanacevars)
    # covvars
      # wide
if (!is.null(covvarswide)) {
svmeancovvarswide1 <- colMeans(usedata[,covvars1, drop = FALSE], na.rm=TRUE)
svmeancovvarswide2 <- colMeans(usedata[,covvars2, drop = FALSE], na.rm=TRUE)
svmeancovvarswide <- rowMeans(cbind(svmeancovvarswide1, svmeancovvarswide2), na.rm=TRUE)
} else {
svmeancovvarswide <- NULL  
}
      # long
if (!is.null(covvarslong_checked)) {
svmeancovvarslong <-  colMeans(usedata[,covvarslong_checked, drop = FALSE], na.rm=TRUE)
} else {
svmeancovvarslong <- NULL  
}
      # all
svmeancovvars <- c(svmeancovvarslong,svmeancovvarswide,svmeancovvarswide)
if (covariance == TRUE) {
svmean <- c(svmeanacevars,svmeanacevars,svmeancovvars)
}

if (covariance == FALSE & !is.null(covvars)) {
svmean <- c(svmeanacevars,svmeanacevars)  
}

cat("\n\n\nStarting Values of the covvars for the mean vector\n\n")
if (!is.null(svmeancovvars)) {
print(svmeancovvars)
} else {
  cat("There are no covariates")
}

cat("\n\n\nStarting Values for the mean vector\n\n")
print(svmean)

if (!is.null(ordinal)) {
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
cat("\n\n\nSummary of used variables\n\n")
print(summary(usedata))
ordinallength <- unlist(lapply(levelslist,length))
#if (2 %in% ordinallength) {
#    stop("Sorry, at the moment the function does not support binary variables")
#}
# Define objects for threshold matrix
nTh       <- ordinallength-1 # No of thresholds
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
}
# define the MZ and DZ data sets
mzData <- usedata %>% filter(!!as.symbol(zyg)==1) %>% select(variables)
dzData <- usedata %>% filter(!!as.symbol(zyg)==2) %>% select(variables)

#mzData    <- subset(usedata, zyg==1, variables)
#dzData    <- subset(usedata, zyg==2, variables)
cat("\n\n\nSummary MZ Data\n\n")
print(summary(mzData))
cat("\n\n\nSummary DZ Data\n\n")
print(summary(dzData))
# Shortcuts for the matrix dimensions
# No. of Variables
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvarsall) # Control variables 
l <- 3*nv*2
t <- m+l+c

# Build elements to construct expected covariance matrix (RAM Notation)

# Matrix A 

# Helper objects for object of manifest paths between acevars f 
mat <- matrix(0.3,nrow = nv,ncol = nv)
mat
freepathB <- lower.tri(mat)
freepathB

mat[upper.tri(mat, diag = TRUE)] <- 0
valuespathB <- mat
valuespathB

nvstring <- as.character(1:nv)
pathBlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("b",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathBlabel[upper.tri(pathBlabel, diag = TRUE)] <- NA
pathBlabel

#pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
#                  free = freepathB,
#                  values = valuespathB,
#                  labels = pathBlabel,
#                  name = "b")
if (type == "aceb") {
pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathB,
                  values = valuespathB,
                  labels = pathBlabel,
                  name = "b")
} else if (type == "chol") {
pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = FALSE,
                  values = 0,
                  labels = pathBlabel,
                  name = "b")  
}
pathZ <- mxMatrix(type = "Zero", nrow = nv, ncol = nv, name = "pZ")

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

if (!is.null(covvars) & covariance == TRUE) {
pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
} 
mat <- matrix(0.3,nrow = nv,ncol = nv)
mat
freepathAC <- lower.tri(mat, diag = TRUE)
freepathAC
mat[upper.tri(mat, diag = FALSE)] <- 0
valuespathAC <- mat
valuespathAC
mat[lower.tri(mat, diag = FALSE)] <- 0
valuespathE <- mat
valuespathE
freepathE <- valuespathE == .3
freepathE


nvstring <- as.character(1:nv)
pathAlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("a",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathAlabel[upper.tri(pathAlabel, diag = FALSE)] <- NA
pathAlabel

pathClabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("c",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathClabel[upper.tri(pathClabel, diag = FALSE)] <- NA
pathClabel

pathElabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("e",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathElabel[upper.tri(pathClabel, diag = FALSE)] <- NA
pathElabel
pathACElb <- diag(0.0001,nv,nv)
pathACElb[pathACElb == 0] <- NA
pathA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  lbound = pathACElb,
                  labels = pathAlabel,
                  name = "a")
pathC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  lbound = pathACElb,
                  labels = pathClabel,
                  name = "c")
if (type == "aceb") {
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathE,
                  values = valuespathE,
                  lbound = pathACElb,
                  labels = pathElabel,
                  name = "e")
} else if (type == "chol") {
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  lbound = pathACElb,
                  labels = pathElabel,
                  name = "e")  
}
pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")
if (!is.null(covvars) & covariance == TRUE) {
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
}
if (is.null(covvars)) {
pathMan <- mxAlgebra(expression = rbind(cbind(b,pZ),
                                              cbind(pZ,b)), name = "pM")    
}
pathACE <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACE")
matA <- mxAlgebra(expression = rbind(cbind(pM,pACE),
                                     Bottom),
                  name = "A")

# MATRIX S
lowerboundcovmat <- function(dimnumber) {
  mat1 <- matrix(NA, dimnumber, dimnumber)
  diag(mat1) <- 0.00001
  mat1[upper.tri(mat1, diag = FALSE)] <- NaN
  mat1 <- as.vector(mat1)
  mat1 <- mat1[!is.nan(mat1)]
  mat1
  return(mat1)
}
labelcovmat <- function(dimlabel) {
#dimlabel <- sapply(strsplit(dimlabel, split = "_", fixed = TRUE), function(x) (x[1]))
labelmat <- matrix("leer",length(dimlabel),length(dimlabel), dimnames = list(dimlabel,dimlabel))
for (r in 1:nrow(labelmat))  {
    for (c in 1:ncol(labelmat)) 
        if (r == c) {
            labelmat[r,c] <- paste("var",dimlabel[r],sep = "_")
        }
        else {
            labelmat[r,c] <- paste("cov",dimlabel[r],dimlabel[c],sep = "_")
    }
}
  labelmat[upper.tri(labelmat, diag = FALSE)] <- NA
  labelmat <- as.vector(labelmat)
  labelmat <- labelmat[!is.na(labelmat)]
  labelmat
return(labelmat)
}
svS <- unname(as.matrix(var(data[,covvarsall], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates

if (!isSymmetric(svS)) {
svS[upper.tri(svS)] <- t(svS)[upper.tri(svS)]
}

if (!is.null(covvars)) {
covCovariates <- mxMatrix(type = "Symm", nrow = c, ncol = c, byrow = FALSE,
                 values = svS,
                 lbound = lowerboundcovmat(c), 
                 labels = labelcovmat(covvarsall),
                 free = TRUE,
                 name = "cCov")
covManCov <- mxMatrix(type = "Zero", nrow = m, ncol = c, name = "cManCov")
covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")
matSMan <- mxAlgebra(expression = rbind(cbind(cMan,cManCov), 
                                        cbind(t(cManCov),cCov),
                                        t(cManCovACE)), name = "matSM")
} 
if (is.null(covvars)) {
covCovariates <- NULL
covManCov <- NULL
covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")
matSMan <- mxAlgebra(expression = rbind(cMan,t(cManCovACE)), name = "matSM")
}

covV <- mxMatrix(type = "Iden", nrow = l/2, ncol = l/2, name = "V")
covCMZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(1,(nv*2)),rep(0,nv)),
                   name = "CMZ")
covCDZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(.5,nv),rep(1,nv),rep(0,nv)),
                   name = "CDZ")

matSMZ <- mxAlgebra(expression = cbind(matSM,rbind(cManCovACE,
                                                   cbind(V,CMZ),
                                                   cbind(CMZ,V))),
                    name = "SMZ")
matSDZ <- mxAlgebra(expression = cbind(matSM,rbind(cManCovACE,
                                                   cbind(V,CDZ),
                                                   cbind(CDZ,V))),
                    name = "SDZ")

filterI <- mxMatrix(type = "Iden", nrow = m+c, ncol = m+c, name = "FI")
filterZ <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "FZ")
matF <- mxAlgebra(expression = cbind(FI,FZ), name = "Filter")

matI <- mxMatrix(type = "Iden", nrow = t, ncol = t, name = "I")
covMZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SMZ%*%t(solve(I-A))%*%t(Filter), name = "expCovMZ")
covDZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SDZ%*%t(solve(I-A))%*%t(Filter), name = "expCovDZ")

###############################################################################
# Mean Matrix
###############################################################################

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

if (covariance == FALSE & !is.null(covvars)) {
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(acevarswide, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = m+l, ncol = 1, 
                 free = c(rep(TRUE,m),rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M") 
}

if (!is.null(ordinal) & (2 %in% ordinallength)) {
binarytrue <- nTh == 1
binaryvar <- ordinalwide[binarytrue]
binaryace <- unlist(lapply(acevarswide,checkcorrespondence, check = binaryvar))
binaryrest <- rep(FALSE,(length(variables)-length(acevarswide)))
if (covariance == TRUE) {
binaryflag <- c(binaryace,binaryrest)
svmean_manifests <- binaryflag == FALSE
matM <- mxMatrix(type = "Full", nrow = t, ncol = 1, 
                 free = c(svmean_manifests,rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M")
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")
}
if (covariance == FALSE & !is.null(covvars)) {
binaryflag <- binaryace
svmean_manifests <- binaryflag == FALSE
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(acevarswide, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = m+l, ncol = 1, 
                 free = c(svmean_manifests,rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M") 
}
}

if (covariance == TRUE) {
# matrix with effect sizes of moderation of the means stored in M
pathCov <- mxMatrix(type = "Full", nrow = c, ncol = ntv, byrow = FALSE,
                            free = t(pathCovfree),
                            values = t(pathCovvalue),
                            labels = t(pathCovlabel),
                            name = "pCov")

# Matrix of definition variables for mean moderation
labeldef <- paste0("data",".",covvarsall)
defM      <- mxMatrix( type="Full", nrow=1, ncol=c, free=FALSE, labels=labeldef, name="defM" )

# Matrix effects on means
effMean <- mxAlgebra(expression = defM%*%pCov, name = "effM")

# Vector of latent variables (all set to zero) just need to concatenate them to the manifests to get the dimensions right
latentmeans <- mxMatrix(type = "Full", nrow = l, ncol = 1, name = "lmeans")

# Concatenate the manifest with the latent means vector
effMeanFull <- mxAlgebra(expression = cbind(effM,lmeans), name = "effMFull")

# Matrix of moderated means
modMean <- mxAlgebra(expression =M+effMFull, name = "modM")

# Matrix of expected means
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%modM), name = "expMean")
}


# Define data object
dataMZ    <- mxData(observed=mzData, type="raw")
dataDZ    <- mxData(observed=dzData, type="raw")

# Variance constraint for binary variables
if (!is.null(ordinal)) {
if (2 %in% ordinallength) {
# fixed variances
# 2 are binary (1st and 3rd)  
  # No of rows = no of binary vars
nrowfilterbinary <- sum(binaryflag)
  # No of cols = vars in total
ncolsfilterbinary <- length(variables)
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

one <- mxMatrix(type = "Unit", nrow = nrowfilterbinary, ncol = 1, name = "Unit")
var1 <- mxConstraint(expression = diag2vec(binCov)==Unit , name = "VConstraint1")
binary <- c(filtermatbin,binarycov,one,var1)
}
}
# Define expectation objects 
if (!is.null(ordinal)) {
expMZ     <- mxExpectationNormal(covariance="expCovMZ", means="expMean", thresholds = "threG", threshnames = ordinalwide,
                                 dimnames=variables)
expDZ     <- mxExpectationNormal(covariance="expCovDZ", means="expMean", thresholds = "threG", threshnames = ordinalwide,
                                 dimnames=variables)
}
else {
expMZ     <- mxExpectationNormal(covariance="expCovMZ", means="expMean",
                                 dimnames=variables)
expDZ     <- mxExpectationNormal(covariance="expCovDZ", means="expMean", 
                                 dimnames=variables)  
}
# Fit function (FIML)
fitfun     <- mxFitFunctionML()

# parameters
pars      <- list(pathB, pathA, pathC, pathE,pathCov, pathZ, matA, pathMan, pathACE,
                  covCovariates, covMan, covManCov, covManCovACE, covV,matSMan,
                  filterI, filterZ, matF, matI,
                  matM, mean, pathBottom)
if (!is.null(ordinal)) {
pars <- c(pars,c(thinG,inc,threG))
}

if (covariance == FALSE & !is.null(covvars)) {
covmeanpars <- c(defM,effMean,latentmeans,effMeanFull,modMean)  
pars <- c(pars,covmeanpars)
}


# group specific model objects
modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun, name="DZ")
if (!is.null(ordinal)) {
if (2 %in% ordinallength) {
  modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, binary,name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun,name="DZ")
}
}
multi     <- mxFitFunctionMultigroup(c("MZ","DZ"))

# overall model object
modelACE  <- mxModel("ACE", pars, modelMZ, modelDZ, multi)
# run model

if (is.null(optimizer) & is.null(ordinal)) {
mxOption(NULL , 'Default optimizer' , 'SLSQP')
} 
else if (is.null(optimizer) & !is.null(ordinal)) {
mxOption(NULL , 'Default optimizer' , 'CSOLNP')   
}
else if (optimizer == "SLSQP") {
mxOption(NULL , 'Default optimizer' , 'SLSQP')    
}
else if (optimizer == "CSOLNP") {
mxOption(NULL , 'Default optimizer' , 'CSOLNP')    
}
else if (optimizer == "NPSOL") {
mxOption(NULL , 'Default optimizer' , 'NPSOL')    
}
set.seed(1)

# Fit full model
#modelACE <- omxAssignFirstParameters(modelACE) # randomly select one starting value if one free parameter has been assigned with more than one starting value

if (tryHard == FALSE){
    fitACE    <- mxRun(modelACE)
}

if (tryHard == TRUE){
if (!is.null(ordinal)) {
fitACE    <- mxTryHardOrdinal(modelACE, extraTries = 10, exhaustive = FALSE)
}
if (is.null(ordinal)) {
fitACE    <- mxTryHard(modelACE, extraTries = 10, exhaustive = FALSE)
}
#fitACE <- mxRun(fitACE)
}
return(fitACE)
}