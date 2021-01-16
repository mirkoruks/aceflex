# Trivariate Cholesky Model
rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)
#data_orig <- data_orig %>% rename(negbez_t1=negbez_1) %>% rename(negbez_t2=negbez_2)  %>% rename(posbez_t1=posbez_1) %>% rename(posbez_t2=posbez_2)


## NEED TO SCALE BEFORE! 
## NEED TO CHECK MISSING DATA STRUCTURE BEFORE! 
  #data <- data[rowSums(is.na(data)) != ncol(data),]
# CONSTANT COVARIATES MUST NOT HAVE A TWIN-SPECIFIC-SUFFIX -> AGE and NOT AGE_1!!

############################################################################################################################################################
############################################################################################################################################################
#################################################### BEGIN OF FUNCTION #####################################################################################
############################################################################################################################################################
############################################################################################################################################################   

# twinflex function
twinflex <- function(acevars, data, zyg, sep, covvars=NULL) {

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
  else {
    print("alles gut")
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
identicalcheck <- as.vector(colSums(ifelse(rna(data_orig[,v1, drop = FALSE])==rna(data_orig[,v2, drop = FALSE]), 0, 1)))
if (0 %in% identicalcheck == TRUE) {
ind <- identicalcheck==0
result1 <- v1[ind]
result2 <- v2[ind]
result <- c(result1,result2)
if (!is.null(result)) {
  stop(c("The following acevars are identical and have no within-variance: ",paste(result, sep = " ", collapse = ", ")))
}
}
else {
  cat("\n\n\nAll the wide-formatted variables have within-pair variance")
}
}
checkvariance(varswide1,varswide2)




usevariables <- c(variables,"zyg")
usedata <- subset(data, select = c(usevariables)) # Data set with only variables used for the model
cat("\n\n\nSummary total Data\n\n")
print(summary(usedata))

# Check if zygosity variable is coded correctly
if (min(data$zyg) != 1 & max(data$zyg) != 2) {
          stop("Zygosity variable must be coded as follows: 1 = MZ, 2 = DZ. Please, recode the zygosity variable.")
}



# Starting Values
  # Means Vector
    # acevars
svmeanacevarswide1 <- colMeans(usedata[,acevars1, drop = FALSE], na.rm=TRUE)
svmeanacevarswide2 <- colMeans(usedata[,acevars2, drop = FALSE], na.rm=TRUE)
svmeanacevars <- rowMeans(cbind(svmeanacevarswide1, svmeanacevarswide2), na.rm=TRUE)
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
svmean <- c(svmeanacevars,svmeanacevars,svmeancovvars)

cat("\n\n\nStarting Values of the covvars for the mean vector\n\n")
if (!is.null(svmeancovvars)) {
print(svmeancovvars)
} else {
  cat("There are no covariates")
}

cat("\n\n\nStarting Values for the mean vector\n\n")
print(svmean)

# define the MZ and DZ data sets
mzData    <- subset(usedata, zyg==1, variables)
dzData    <- subset(usedata, zyg==2, variables)
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

pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathB,
                  values = valuespathB,
                  labels = pathBlabel,
                  name = "b")
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

if (!is.null(covvars)) {
pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
} else {
pathCov <- mxMatrix(type = "Full", nrow = 0, ncol = 0, byrow = FALSE,
                            name = "pCov")  
}
print(pathCov)
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

pathA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  labels = pathAlabel,
                  name = "a")
pathC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  labels = pathClabel,
                  name = "c")
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathE,
                  values = valuespathE,
                  labels = pathElabel,
                  name = "e")
pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
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
svS <- unname(as.matrix(var(data_orig[,covvarsall], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates

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
} else {
covCovariates <- mxMatrix(type = "Symm", nrow = 0, ncol = 0, byrow = FALSE,
                 name = "cCov")
}

print(covCovariates)

covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
covManCov <- mxMatrix(type = "Zero", nrow = m, ncol = c, name = "cManCov")
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")

covV <- mxMatrix(type = "Iden", nrow = l/2, ncol = l/2, name = "V")
covCMZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(1,(nv*2)),rep(0,nv)),
                   name = "CMZ")
covCDZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(.5,nv),rep(1,nv),rep(0,nv)),
                   name = "CDZ")
matSMan <- mxAlgebra(expression = rbind(cbind(cMan,cManCov), 
                                        cbind(t(cManCov),cCov),
                                        t(cManCovACE)), name = "matSM")

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
print(variables)
print(svmean)

# Mean Matrix
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(variables, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = t, ncol = 1, 
                 free = c(rep(TRUE,m+c),rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M")
print(matM)
print(length(c(meanmanifestlabel)))
print(length(c(meanacelabel)))
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")

cat("\n\n\n... Bis hierhin läuft alles durch! Weiter geht's! :-)")
}
############################################################################################################################################################
############################################################################################################################################################
#################################################### END OF FUNCTION #######################################################################################
############################################################################################################################################################
############################################################################################################################################################

twinflex(acevars = c("kultkapjahre"),covvars = c("age","posbez"),data = data_orig,sep = "_",zyg = "zyg")

sep <- "_"
acevars <- c("kultkapjahre","negbez")
covvars <- c("age","posbez")

# check if vars are constant

#######
  # 1. check if variables exist --> for acevars in wide format and for covvars in both long (without within-pair variance) and wide (with within-pair variance)
#######

existence <- function(variable) {
  result <- NULL
if(variable %in% colnames(data_orig))
{
  result <- NULL
}
else if (!(variable %in% colnames(data_orig))) {
  result <- variable
}
}

existenceerror <- function(result) {
  if (!is.null(result)) {
  stop(c("I could not find a variable(s) in the data frame corresponding to the following variable strings you gave me: ",paste(result, sep = " ", collapse = ", ")))
  
  }
  else {
    print("alles gut")
  }
  }


  # 1a: Check for acevars (only in wide since they have to have within-pair-variance)
acevars1 <-    paste0(acevars,sep,"1") # ACE vars twin 1
acevars2 <-    paste0(acevars,sep,"2") # ACE vars twin 2
acevarswide <- c(acevars1, acevars2) # ACE vars variable vector (input for SEM)#
existence_check_acevars <- unlist(lapply(acevarswide, existence))
existence_check_acevars

existenceerror(existence_check_acevars)

acevars_fake1 <- c(acevars1,"fake_t1")
acevars_fake2 <- c(acevars2,"fake_t2")
acevarswide_fake <- c(acevars_fake1,acevars_fake2)
existence_check_acevars <- unlist(lapply(acevarswide_fake, existence))
acevars_not_found <- unlist(sapply(strsplit(existence_check_acevars, split=sep, fixed=TRUE), function(x) (x[1])))
acevars_not_found <- unique(acevars_not_found)
existenceerror(acevars_not_found)
  
  # acevars mit falscher acevar
    # acevarswide1 <- c(acevarswide,"falsch_1","falsch_2")
    # existence_check_acevars <- unlist(lapply(acevarswide1, existence))
    # existenceerror(existence_check_acevars)

  # 1b: Check for covvars (can be in wide and long) -> 1. Check for wide -> 2. Check for Long if something does not appear as wide
covvars1 <-    paste0(covvars,sep,"1") # Covariates twin 1
covvars2 <-    paste0(covvars,sep,"2") # Covariates twin 2
covvarswide <- c(covvars1, covvars2)

existence_check_covvars1 <- unlist(lapply(covvarswide, existence))
covvarswide_checked <- covvarswide[!covvarswide %in% existence_check_covvars1] # the covariates which exist in a wide format
existence_check_covvars1 # for these variables there are no wide formatted variables in the data frame -> check if they are in long format

if (!is.null(existence_check_covvars1)) { # if the condition is not fulfilled, then I could not find some of the wide variables -> maybe they are long variables!
covvars_possibly_long_format <- covvarswide[covvarswide %in% existence_check_covvars1] 
covvars_possibly_long_format <- unlist(sapply(strsplit(covvars_possibly_long_format, split=sep, fixed=TRUE), function(x) (x[1])))
covvars_possibly_long_format <- unique(covvars_possibly_long_format)
existence_check_covvars2 <- unlist(lapply(covvars_possibly_long_format, existence))
existenceerror(existence_check_covvars2)
covvarslong_checked <- covvars_possibly_long_format
} else {
covvarslong_checked <- NULL # if the condition is fulfilled, then all the covariates are in wide format and none in long
}
if (length(covvarswide_checked)==0) {
  covvarswide_checked <- NULL
}

covvarsall <- c(covvarslong_checked,covvarswide_checked) # from now on: if a covariate ends with _1 or _2 -> it is wide!
covvars1 <- covvarswide_checked[grepl('_1', covvarswide_checked)]
covvars2 <- covvarswide_checked[grepl('_2', covvarswide_checked)]
covvarswide <- c(covvars1, covvars2)

covvars <- c(covvarslong_checked,covvarswide_checked) # from now on: if a covariate ends with _1 or _2 -> it is wide!
covvars1 <- covvars[grepl('_1', covvars)]
covvars2 <- covvars[grepl('_2', covvars)]
covvarswide <- c(covvars1, covvars2)
covvars
covvars1
covvars2
covvarswide

  
# 2. check if the wide formatted variables have within-pair-variance
  # acevars
    # generate some fake variables
data_orig$fakewide_1 <- runif(dim(data_orig)[1]) 
ind <- which(data_orig$fakewide_1 %in% sample(data_orig$fakewide_1, 15))
data_orig$fakewide_1[ind]<-NA
data_orig$fakewide_2 <- data_orig$fakewide_1

acevars1
acevars2
acevars1_fake <- c(acevars1,"fakewide_1")
acevars2_fake <- c(acevars2,"fakewide_2")
acevars1_fake
acevars2_fake

rna <- function(x) replace(x, is.na(x), "")
checkvariance <- function(v1,v2) {
identicalcheck <- as.vector(colSums(ifelse(rna(data_orig[,v1, drop = FALSE])==rna(data_orig[,v2, drop = FALSE]), 0, 1)))
if (0 %in% identicalcheck == TRUE) {
ind <- identicalcheck==0
result1 <- v1[ind]
result2 <- v2[ind]
result <- c(result1,result2)
if (!is.null(result)) {
  stop(c("The following acevars are identical and have no within-variance: ",paste(result, sep = " ", collapse = ", ")))
}
}
else {
  print("All the acevars have within-pair variance")
}
}

# Die richtigen Vars
as.vector(colSums(ifelse(rna(data_orig[,acevars1, drop = FALSE])==rna(data_orig[,acevars2, drop = FALSE]), 0, 1)))
checkvariance(acevars1,acevars2)
# Die fake Vars
as.vector(colSums(ifelse(rna(data_orig[,acevars1, drop = FALSE])==rna(data_orig[,acevars2, drop = FALSE]), 0, 1)))
checkvariance(acevars1_fake,acevars2_fake)

variables <- c(acevarswide,covvars)
variables
# Starting Values
  # Means Vector
    # acevars
svmeanacevarswide1 <- colMeans(data_orig[,acevars1, drop = FALSE], na.rm=TRUE)
svmeanacevarswide2 <- colMeans(data_orig[,acevars2, drop = FALSE], na.rm=TRUE)
svmeanacevars <- rowMeans(cbind(svmeanacevarswide1, svmeanacevarswide2), na.rm=TRUE)
cat("\n\n\nStarting Values of the acevars for the mean vector\n\n")
print(svmeanacevars)
    # covvars
      # wide
if (!is.null(covvarswide)) {
svmeancovvarswide1 <- colMeans(data_orig[,covvars1, drop = FALSE], na.rm=TRUE)
svmeancovvarswide2 <- colMeans(data_orig[,covvars2, drop = FALSE], na.rm=TRUE)
svmeancovvarswide <- rowMeans(cbind(svmeancovvarswide1, svmeancovvarswide2), na.rm=TRUE)
} else {
svmeancovvarswide <- NULL  
}
      # long
if (!is.null(covvarslong_checked)) {
svmeancovvarslong <-  colMeans(data_orig[,covvarslong_checked, drop = FALSE], na.rm=TRUE)
} else {
svmeancovvarslong <- NULL  
}
      # all
svmeancovvars <- c(svmeancovvarslong,svmeancovvarswide,svmeancovvarswide)
svmean <- c(svmeanacevars,svmeanacevars,svmeancovvars)
variables
svmean
test <- mxMatrix(type = "Full",nrow = 1, ncol = length(svmean), values = svmean)
test
# Matrix A 
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvars) # Control variables 
l <- 3*nv*2
t <- m+l+c

paste0(c("A_","C_","E_"),"variable")

acedimlabel <- function(var) {
c(paste0("A_",var),paste0("C_",var),paste0("E_",var))
}
as.vector(unlist(sapply(acevarswide,acedimlabel)))


rep(paste0(c("A_","C_"),acevars1),each = nv)
       
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

########################
######### COVARIATES
########################
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvars) # Control variables 
l <- 3*nv*2
t <- m+l+c


### Cov-Vars with Variance

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
pathCovlabelconstant <- as.matrix(sapply(covvars,pathCov_label_constant))
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

pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
print(pathCov)
test <- mxMatrix(type = "Full",nrow = 0, ncol = 0)






pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                    values = c(rep(.5,ntv),
                               rep(.5,nv),rep(0,nv),
                               rep(.5,nv),rep(0,nv),
                               rep(0,nv),rep(.5,nv),
                               rep(0,nv),rep(.5,nv)),
                    labels = c(c("bage1","bage2"),c("bage1","bage2"),
                               c("bdo11","bdo21"),rep(NA,nv),
                               c("bdy11","bdy21"),rep(NA,nv),
                               rep(NA,nv),c("bdo12","bdo22"),
                               rep(NA,nv),c("bdy12","bdy22")),
                    free = c(rep(TRUE,ntv),
                               rep(TRUE,nv),rep(FALSE,nv),
                               rep(TRUE,nv),rep(FALSE,nv),
                               rep(FALSE,nv),rep(TRUE,nv),
                               rep(FALSE,nv),rep(TRUE,nv)),
                    name = "pCov")

covvars <- "age"
acevars <- c("posbez")
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvars) # Control variables 
l <- 3*nv*2
t <- m+l+c

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

pathA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  labels = pathAlabel,
                  name = "a")
pathC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathAC,
                  values = valuespathAC,
                  labels = pathClabel,
                  name = "c")
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathE,
                  values = valuespathE,
                  labels = pathElabel,
                  name = "e")
print(pathA)
print(pathC)
print(pathE)

pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
pathACE <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACE")
matA <- mxAlgebra(expression = rbind(cbind(pM,pACE),
                                     Bottom),
                  name = "A")


covvars = c("age","zyg")
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
svS <- unname(as.matrix(cov(data_orig[,covvars], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates
if (!isSymmetric(svS)) {
svS[upper.tri(svS)] <- t(svS)[upper.tri(svS)]
}

covCovariates <- mxMatrix(type = "Symm", nrow = c, ncol = c, byrow = FALSE,
                 values = svS,
                 lbound = lowerboundcovmat(c), 
                 labels = labelcovmat(covvars),
                 free = TRUE,
                 name = "cCov")
print(covCovariates)
cov(data_orig[,"age"])

##############################################################################
##############################################################################
##############################################################################
# next step: Matrix S 
##############################################################################
##############################################################################
##############################################################################

covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
covManCov <- mxMatrix(type = "Zero", nrow = m, ncol = c, name = "cManCov")
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")

covV <- mxMatrix(type = "Iden", nrow = l/2, ncol = l/2, name = "V")
covCMZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(1,(nv*2)),rep(0,nv)),
                   name = "CMZ")
covCDZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(.5,nv),rep(1,nv),rep(0,nv)),
                   name = "CDZ")
matSMan <- mxAlgebra(expression = rbind(cbind(cMan,cManCov), 
                                        cbind(t(cManCov),cCov),
                                        t(cManCovACE)), name = "matSM")

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

# Mean Matrix
matM <- mxMatrix(type = "Full", nrow = t, ncol = 1, 
                 free = c(rep(TRUE,m+c),rep(FALSE,l)), 
                 labels = c(rep(c("int1","int2"),2),"meanage",c("meando","meandy","meando","meandy"),rep(NA,l)),
                 values = c(rep(0,(m)),svM,rep(0,l)), 
                 name = "M")
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")

# Define data object
dataMZ    <- mxData(observed=mzData, type="raw")
dataDZ    <- mxData(observed=dzData, type="raw")

# Define expectation objects 
expMZ     <- mxExpectationNormal(covariance="expCovMZ", means="expMean",
                                 dimnames=vars_cov)
expDZ     <- mxExpectationNormal(covariance="expCovDZ", means="expMean", 
                                 dimnames=vars_cov)
# Fit function (FIML)
fitfun     <- mxFitFunctionML()

# parameters
pars      <- list(pathB, pathA, pathC, pathE,pathCov, pathZ, matA, pathMan, pathACE,
                  covCovariates, covMan, covManCov, covManCovACE, covV,matSMan,
                  filterI, filterZ, matF, matI,
                  matM, mean, pathBottom)

# group specific model objects
modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun, name="DZ")
multi     <- mxFitFunctionMultigroup(c("MZ","DZ"))

# overall model object
modelACE  <- mxModel("ACE", pars, modelMZ, modelDZ, multi)

# run model
#mxOption(NULL , 'Default optimizer' , 'NPSOL')
#mxOption(NULL , 'Default optimizer' , 'CSOLNP')
mxOption(NULL , 'Default optimizer' , 'SLSQP')


set.seed(1)

# Fit full model
modelACE <- omxAssignFirstParameters(modelACE) # randomly select one starting value if one free parameter has been assigned with more than one starting value
#startACE <- mxAutoStart(modelACE)
fitACE    <- mxTryHard(modelACE, extraTries = 50, exhaustive = TRUE)
fitACE <- mxRun(fitACE)
# Summarize model
sumACE    <- summary(fitACE) 
sumACE

# Check identification status
fitACEIdent <- mxCheckIdentification(fitACE)
fitACEIdent$status
fitACEIdent$non_identified_parameters

