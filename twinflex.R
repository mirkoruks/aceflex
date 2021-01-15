# Trivariate Cholesky Model
rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)


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
existenceerror(existence_check_acevars)

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
svmeancovvars <- c(svmeancovvarslong,svmeancovvarswide)
svmean <- c(svmeanacevars,svmeancovvars)

cat("\n\n\nStarting Values of the covvars for the mean vector\n\n")
if (!is.null(svmeancovvars)) {
print(svmeancovvars)
} else {
  cat("There are no covariates")
}

cat("\n\n\nStarting Values for the mean vector\n\n")
print(svmean)

  # Covariance Matrix of the covariates (at the moment the function uses default values)
if (!is.null(covvars)) { 
# here comes the code ! :-)
  # here something to begin with
    # unname(as.matrix(cor(data[,covariates], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates
    # svS <- round(unname(as.matrix(cov(data[,covariates], use = "na.or.complete"))),3) # S matrix starting values for non-decomposed covariates
    # svS
}

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

pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
print(pathCov)
cat("\n\n\n... Bis hierhin läuft alles durch! Weiter geht's! :-)")
}
############################################################################################################################################################
############################################################################################################################################################
#################################################### END OF FUNCTION #######################################################################################
############################################################################################################################################################
############################################################################################################################################################

twinflex(acevars = c("posbez","kultkapjahre"), covvars = c("age","negbez"),data = data_orig,sep = "_",zyg = "zyg")

varswide1 <- c(acevars1,NULL)
svmeanacevars1 <- colMeans(data_orig[,acevars1], na.rm=TRUE)
svmeanacevars2 <- colMeans(data_orig[,acevars2], na.rm=TRUE)
acevars_svmean <- rowMeans(cbind(svmeanacevars1,svmeanacevars2), na.rm=TRUE)


sep <- "_"
# Build elements to construct expected covariance matrix (RAM Notation)
acevars <- c("posbez", "schoolhigh", "iseiempmean")

acevars1 <-    paste0(acevars,sep,"1") # ACE vars twin 1
acevars2 <-    paste0(acevars,sep,"2") # ACE vars twin 2
acevarswide <- c(acevars1, acevars2) # ACE vars variable vector (input for SEM)

covvars <- c("age","agediff")
covvars1 <-    paste0(covvars,sep,"1") # Covariates twin 1
covvars2 <-    paste0(covvars,sep,"2") # Covariates twin 2
covvarswide <- c(covvars1, covvars2)

vars1 <- c(acevars1,covvars1) # All variables for twin 1 
vars2 <- c(acevars2,covvars2) # All variables for twin 2



varswide <- c(acevarswide,covvarswide) # all variables wide format

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
existence_check_acevars <- unlist(lapply(acevarswide, existence))
existenceerror(existence_check_acevars)

acevars_fake1 <- c(acevars1,"fake_1")
acevars_fake2 <- c(acevars2,"fake_2")
acevarswide_fake <- c(acevars_fake1,acevars_fake2)
existence_check_acevars <- unlist(lapply(acevarswide_fake, existence))
existenceerror(existence_check_acevars)
  
  # acevars mit falscher acevar
    # acevarswide1 <- c(acevarswide,"falsch_1","falsch_2")
    # existence_check_acevars <- unlist(lapply(acevarswide1, existence))
    # existenceerror(existence_check_acevars)

  # 1b: Check for covvars (can be in wide and long) -> 1. Check for wide -> 2. Check for Long if something does not appear as wide
covvars <- c("age")
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

 as.vector(colSums(ifelse(rna(data_orig[,acevars1[1], drop = FALSE])==rna(data_orig[,acevars2[1], drop = FALSE]), 0, 1)))
 
a1 <- ifelse(rna(data_orig[,acevars1_fake])==rna(data_orig[,acevars2_fake]), 0, 1)
a2 <- colSums(ifelse(rna(data_orig[,acevars1_fake])==rna(data_orig[,acevars2_fake]), 0, 1))

colSums(a1)
a2


# Matrix A 
nv <- length(acevars) # Vars per twin
ntv <- nv*2 # Vars per twin pair
m <- (nv*2) # Decomposed manifest variables
c <- length(covvars) # Control variables 
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

pathA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = TRUE,
                  values = .8,
                  labels = c("a11","a21","a22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "a")
pathC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = TRUE,
                  values = .8,
                  labels = c("c11","c21","c22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "c")
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = c(TRUE,FALSE,TRUE),
                  values = c(.8,0,.8),
                  labels = c("e11","e21","e22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "e")
pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
pathACE <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACE")
matA <- mxAlgebra(expression = rbind(cbind(pM,pACE),
                                     Bottom),
                  name = "A")

# Matrix S 

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
            labelmat[r,c] <- paste0("var",dimlabel[r])
        }
        else {
            labelmat[r,c] <- paste0("cov",dimlabel[r],dimlabel[c])
    }
}
  labelmat[upper.tri(labelmat, diag = FALSE)] <- NA
  labelmat <- as.vector(labelmat)
  labelmat <- labelmat[!is.na(labelmat)]
  labelmat
return(labelmat)
}

labelcov <- c("vAge","cAgeDo1","cAgeDy1","cAgeDo2","cAgeDy2",
                  "vDo1","cDo1Dy1","cDo1Do2","cDo1Dy2",
                  "vDy1","cDo2Dy1","cDy1Dy2",
                  "vDo2","cDo2Dy2",
                  "vDy2")


covCovariates <- mxMatrix(type = "Symm", nrow = c, ncol = c, byrow = FALSE,
                 values = svS,
                 lbound = lowerboundcovmat(c), 
                 labels = labelcov,
                 free = TRUE,
                 name = "cCov")
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



# Parameter constraints 

# Covariates
modelACEred_11   <- mxModel(fitACE, name="ACEreduced_11")
ACEred_11  <- omxSetParameters(modelACEred_11, labels=c("bdo11","bdo12"), free=TRUE, values=--.14, newlabels='bdo1')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdo21","bdo22"), free=TRUE, values=-.025, newlabels='bdo2')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdy11","bdy12"), free=TRUE, values=--.15, newlabels='bdy1')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdy21","bdy22"), free=TRUE, values=-.04, newlabels='bdy2')

ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cAgeDo1","cAgeDo2"), free=TRUE, values=-0.028, newlabels='cAgeDo')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cAgeDy1","cAgeDy2"), free=TRUE, values=-.009, newlabels='cAgeDy')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cDo1Dy1","cDo2Dy1","cDo1Dy2","cDo2Dy2"), free=TRUE, values=-0.25, newlabels="cDoDy")
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("vDo1","vDo2"), free=TRUE, values=.99, newlabels='vDo')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("vDy1","vDy2"), free=TRUE, values=.99, newlabels='vDy')
fitACEred_11     <- mxTryHard(ACEred_11, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_11) 
mxCompare(fitACE, fitACEred_11) # p = 0.93 -> OK!

modelACEred_12   <- mxModel(fitACEred_11, name="ACEreduced_12")
ACEred_12   <- omxSetParameters(modelACEred_12, labels=c("bdo2","bdy2","cAgeDo","cAgeDy"), free=FALSE, values=0)
fitACEred_12     <- mxTryHard(ACEred_12, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_12) 
mxCompare(fitACEred_11, fitACEred_12) # p = 0.64 -> OK!

# c22
modelACEred_21   <- mxModel(fitACEred_12, name="ACEreduced_21")
ACEred_21   <- omxSetParameters(modelACEred_21, labels=c("c22"), free=FALSE, values=0)
fitACEred_21     <- mxTryHard(ACEred_21, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_21) 
mxCompare(fitACEred_12, fitACEred_21) # p = 0.63 -> OK!

# b21
modelACEred_22   <- mxModel(fitACEred_21, name="ACEreduced_22")
ACEred_22   <- omxSetParameters(modelACEred_22, labels=c("b21"), free=FALSE, values=0)
fitACEred_22     <- mxTryHard(ACEred_22, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_22) 
mxCompare(fitACEred_21, fitACEred_22) # p = 0.20 -> OK!

# a21
modelACEred_31   <- mxModel(fitACEred_22, name="ACEreduced_31")
ACEred_31   <- omxSetParameters(modelACEred_31, labels=c("a21"), free=FALSE, values=0)
fitACEred_31     <- mxTryHard(ACEred_31, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_31) 
mxCompare(fitACEred_22, fitACEred_31) # p < 0.05 -> Not OK!

# Best fitting model
fitBestModel <- fitACEred_22
sumBestModel <- summary(fitBestModel) 
sumBestModel

# Generate Output
sig.level <- function(coef,std.error){
  ifelse(abs(coef/std.error)>= qnorm(.0005,lower.tail=FALSE),"$^{***}$",
         ifelse(abs(coef/std.error) >= qnorm(.005,lower.tail=FALSE), "$^{**}$",
                ifelse(abs(coef/std.error) >= qnorm(.025,lower.tail=FALSE), "$^{*}$","")
        )
 )
}
nobs <- c("N",as.character(dim(data)[1]))
loli <- c("-2LL",as.character(round(sumBestModel[["Minus2LogLikelihood"]], digits = 3)))
aic <- c("AIC",as.character(round(sumBestModel[["informationCriteria"]][1,3], digits = 3)))
bic <- c("BIC",as.character(round(sumBestModel[["informationCriteria"]][2,3], digits = 3)))
dvname <- c("IV",iv)

# How many parameters are there for the table?
np <- 7
ACE.printcoefs <-
  function(modelfit){
    namelist <- rbind(names(modelfit@output[["estimate"]][1:np]),
                      rep("",length(names(modelfit@output[["estimate"]][1:np]))))
    coeflist <- rbind(paste(round(modelfit@output[["estimate"]][1:np],3),
                            sig.level(modelfit@output[["estimate"]][1:np],
                                      modelfit@output[["standardErrors"]][1:np]),sep=""),
                      paste("(",round(modelfit@output[["standardErrors"]][1:np],3),")",sep="")
   )
    list(coefs = rbind(as.vector(dvname),cbind(as.vector(namelist),as.vector(coeflist)),as.vector(nobs),as.vector(loli),as.vector(aic),as.vector(bic)))
  }

output[[iv]] <- ACE.printcoefs(fitBestModel)$coefs

output_full <- do.call(cbind, output)

print.xtable(xtable(output_full,
                    caption = paste0("Bivariate ACE-$\\beta$ model ",iv)),type='latex', sanitize.text.function=identity,
             caption.placement = "top", size="\\tiny", table.placement = "H",
             file = paste0("C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/Output/","biv_",iv,"2dummy",".tex"), 
             fileEncoding = "UTF-8",include.rownames = FALSE, include.colnames = FALSE)

=======
# Trivariate Cholesky Model
rm(list = ls())

library(dplyr)
#library(OpenMx)
library(xtable)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)

acevars <- c("posbez", "schoolhigh", "iseiempmean")
acevars1 <-    paste0(acevars,"_1") # ACE vars twin 1
vars1 <- c(acevars1) # All variables for twin 1 
acevars2 <-    paste0(acevars,"_2") # ACE vars twin 2
vars2 <- c(acevars2) # All variables for twin 2
acevarswide <- c(acevars1, acevars2) # ACE vars variable vector (input for SEM)
varswide <- c(vars1,vars2)

## NEED TO SCALE BEFORE! 
## NEED TO CHECK MISSING DATA STRUCTURE BEFORE! 
  #data <- data[rowSums(is.na(data)) != ncol(data),]

#####################################################################################################     
#####################################################################################################     
#####################################################################################################     
#####################################################################################################     

# twinflex function
twinflex <- function(acevars, data, zyg, covvars=NULL) {

if ("OpenMx" %in% (.packages()) == FALSE) {
  stop("You need to load the OpenMx library: \nuse library(OpenMx)\n...Science is standing on the shoulders of giants and so does this function... :-)")
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
  
if (!is.null(covvars) & class(covvars)!= "character") {
  stop("Please, use a string vector to specify the covariates")
}

# Output-List -> Print results
output <- list()   

## check here if covariates are constants or not! 

# Transform vars from long to wide
acevars1 <-    paste0(acevars,"_1") # ACE vars twin 1
acevars2 <-    paste0(acevars,"_2") # ACE vars twin 2

if (is.null(covvars)) {
  covvars1 <- covvars
  covvars2 <- covvars
}
# here comes the else if option for constant covariates!
else {
covvars1 <-    paste0(covvars,"_1") # Covariates twin 1
covvars2 <-    paste0(covvars,"_2") # Covariates twin 2 
}


vars1 <- c(acevars1,covvars1) # All variables for twin 1 
vars2 <- c(acevars2,covvars2) # All variables for twin 2


variables <- c(acevars1, acevars2, covvars1, covvars2) # ACE vars variable vector (input for SEM)

# Check if acevars have within-twin-pair-variance -> give warning if correlation > .9
#####################################################################################################     
#### When function allows for covariates in the covariance matrix we can add them here as well ! #### 
#####################################################################################################     

acevarscor <- mapply(cor,data[,acevars1],data[,acevars2], use = "pairwise.complete.obs")
#acevarscor[3] <- .98
s = attr(acevarscor, "names")
s1 = unlist(sapply(strsplit(s, split='_', fixed=TRUE), function(x) (x[1])))
attr(acevarscor, "names") <- s1

lowcorrelation <- acevarscor < 0.9 
alllowcorrelation <- all(lowcorrelation)
if (alllowcorrelation == FALSE) {
    highcorrelation <- acevarscor[acevarscor >=.9]
    highcorrelationnames <- attributes(highcorrelation)
    print("Oups! For the following variables the within-pair correlation is >= 0.9 and < 1. There might be estimation problems due to (multi-)collinearity")
    print(highcorrelation)
    proceedcollinearity <- readline("Do you want to proceed? You still want to proceed? Then type 'yes' You want to stop? Then type 'no' (Without quotes)") 
    if (proceedcollinearity == "no") {
        stop("User did not want to proceed due to (multi-)collinearity in one of the variables. Function twinflex has stopped! See you soon!")
    }
}

# Check if covariates are constants or have within-pair-variance
#
#
#
#
#


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
svmean1 <- colMeans(usedata[,vars1], na.rm=TRUE)
svmean2 <- colMeans(usedata[,vars2], na.rm=TRUE)
svmean <- rowMeans(cbind(svmean1,svmean2), na.rm=TRUE)
cat("\n\n\nStarting Values for the mean vector\n\n")
print(svmean)

  # Covariance Matrix of the covariates (at the moment the function uses default values)
if (!is.null(covvars)) { 
# here comes the code ! :-)
  # here something to begin with
    # unname(as.matrix(cor(data[,covariates], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates
    # svS <- round(unname(as.matrix(cov(data[,covariates], use = "na.or.complete"))),3) # S matrix starting values for non-decomposed covariates
    # svS
}

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
c <- length(acevars) # Control variables 
l <- 3*nv*2
t <- m+l+c

# Build elements to construct expected covariance matrix (RAM Notation)

# Matrix A 

# Helper objects for object of manifest paths between acevars
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
print(pathB)
cat("\n\n\n... Bis hierhin läuft alles durch! Weiter geht's! :-)")
}
 
twinflex(acevars = acevars, data = data_orig,zyg = "zyg")




# Build elements to construct expected covariance matrix (RAM Notation)

# Matrix A 
pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathB,
                  values = valuespathB,
                  labels = pathBlabel,
                  name = "b")
pathB
pathZ <- mxMatrix(type = "Zero", nrow = nv, ncol = nv, name = "pZ")

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

pathA <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = TRUE,
                  values = .8,
                  labels = c("a11","a21","a22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "a")
pathC <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = TRUE,
                  values = .8,
                  labels = c("c11","c21","c22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "c")
pathE <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = c(TRUE,FALSE,TRUE),
                  values = c(.8,0,.8),
                  labels = c("e11","e21","e22"),
                  lbound = c(0.000001,NA,0.000001),
                  name = "e")
pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
pathACE <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACE")
matA <- mxAlgebra(expression = rbind(cbind(pM,pACE),
                                     Bottom),
                  name = "A")

# Matrix S 

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
            labelmat[r,c] <- paste0("var",dimlabel[r])
        }
        else {
            labelmat[r,c] <- paste0("cov",dimlabel[r],dimlabel[c])
    }
}
  labelmat[upper.tri(labelmat, diag = FALSE)] <- NA
  labelmat <- as.vector(labelmat)
  labelmat <- labelmat[!is.na(labelmat)]
  labelmat
return(labelmat)
}

labelcov <- c("vAge","cAgeDo1","cAgeDy1","cAgeDo2","cAgeDy2",
                  "vDo1","cDo1Dy1","cDo1Do2","cDo1Dy2",
                  "vDy1","cDo2Dy1","cDy1Dy2",
                  "vDo2","cDo2Dy2",
                  "vDy2")


covCovariates <- mxMatrix(type = "Symm", nrow = c, ncol = c, byrow = FALSE,
                 values = svS,
                 lbound = lowerboundcovmat(c), 
                 labels = labelcov,
                 free = TRUE,
                 name = "cCov")
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



# Parameter constraints 

# Covariates
modelACEred_11   <- mxModel(fitACE, name="ACEreduced_11")
ACEred_11  <- omxSetParameters(modelACEred_11, labels=c("bdo11","bdo12"), free=TRUE, values=--.14, newlabels='bdo1')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdo21","bdo22"), free=TRUE, values=-.025, newlabels='bdo2')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdy11","bdy12"), free=TRUE, values=--.15, newlabels='bdy1')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("bdy21","bdy22"), free=TRUE, values=-.04, newlabels='bdy2')

ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cAgeDo1","cAgeDo2"), free=TRUE, values=-0.028, newlabels='cAgeDo')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cAgeDy1","cAgeDy2"), free=TRUE, values=-.009, newlabels='cAgeDy')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("cDo1Dy1","cDo2Dy1","cDo1Dy2","cDo2Dy2"), free=TRUE, values=-0.25, newlabels="cDoDy")
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("vDo1","vDo2"), free=TRUE, values=.99, newlabels='vDo')
ACEred_11  <- omxSetParameters(ACEred_11, labels=c("vDy1","vDy2"), free=TRUE, values=.99, newlabels='vDy')
fitACEred_11     <- mxTryHard(ACEred_11, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_11) 
mxCompare(fitACE, fitACEred_11) # p = 0.93 -> OK!

modelACEred_12   <- mxModel(fitACEred_11, name="ACEreduced_12")
ACEred_12   <- omxSetParameters(modelACEred_12, labels=c("bdo2","bdy2","cAgeDo","cAgeDy"), free=FALSE, values=0)
fitACEred_12     <- mxTryHard(ACEred_12, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_12) 
mxCompare(fitACEred_11, fitACEred_12) # p = 0.64 -> OK!

# c22
modelACEred_21   <- mxModel(fitACEred_12, name="ACEreduced_21")
ACEred_21   <- omxSetParameters(modelACEred_21, labels=c("c22"), free=FALSE, values=0)
fitACEred_21     <- mxTryHard(ACEred_21, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_21) 
mxCompare(fitACEred_12, fitACEred_21) # p = 0.63 -> OK!

# b21
modelACEred_22   <- mxModel(fitACEred_21, name="ACEreduced_22")
ACEred_22   <- omxSetParameters(modelACEred_22, labels=c("b21"), free=FALSE, values=0)
fitACEred_22     <- mxTryHard(ACEred_22, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_22) 
mxCompare(fitACEred_21, fitACEred_22) # p = 0.20 -> OK!

# a21
modelACEred_31   <- mxModel(fitACEred_22, name="ACEreduced_31")
ACEred_31   <- omxSetParameters(modelACEred_31, labels=c("a21"), free=FALSE, values=0)
fitACEred_31     <- mxTryHard(ACEred_31, extraTries = 50, exhaustive = TRUE)
# summary(fitACEred_31) 
mxCompare(fitACEred_22, fitACEred_31) # p < 0.05 -> Not OK!

# Best fitting model
fitBestModel <- fitACEred_22
sumBestModel <- summary(fitBestModel) 
sumBestModel

# Generate Output
sig.level <- function(coef,std.error){
  ifelse(abs(coef/std.error)>= qnorm(.0005,lower.tail=FALSE),"$^{***}$",
         ifelse(abs(coef/std.error) >= qnorm(.005,lower.tail=FALSE), "$^{**}$",
                ifelse(abs(coef/std.error) >= qnorm(.025,lower.tail=FALSE), "$^{*}$","")
        )
 )
}
nobs <- c("N",as.character(dim(data)[1]))
loli <- c("-2LL",as.character(round(sumBestModel[["Minus2LogLikelihood"]], digits = 3)))
aic <- c("AIC",as.character(round(sumBestModel[["informationCriteria"]][1,3], digits = 3)))
bic <- c("BIC",as.character(round(sumBestModel[["informationCriteria"]][2,3], digits = 3)))
dvname <- c("IV",iv)

# How many parameters are there for the table?
np <- 7
ACE.printcoefs <-
  function(modelfit){
    namelist <- rbind(names(modelfit@output[["estimate"]][1:np]),
                      rep("",length(names(modelfit@output[["estimate"]][1:np]))))
    coeflist <- rbind(paste(round(modelfit@output[["estimate"]][1:np],3),
                            sig.level(modelfit@output[["estimate"]][1:np],
                                      modelfit@output[["standardErrors"]][1:np]),sep=""),
                      paste("(",round(modelfit@output[["standardErrors"]][1:np],3),")",sep="")
   )
    list(coefs = rbind(as.vector(dvname),cbind(as.vector(namelist),as.vector(coeflist)),as.vector(nobs),as.vector(loli),as.vector(aic),as.vector(bic)))
  }

output[[iv]] <- ACE.printcoefs(fitBestModel)$coefs

output_full <- do.call(cbind, output)

print.xtable(xtable(output_full,
                    caption = paste0("Bivariate ACE-$\\beta$ model ",iv)),type='latex', sanitize.text.function=identity,
             caption.placement = "top", size="\\tiny", table.placement = "H",
             file = paste0("C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/Output/","biv_",iv,"2dummy",".tex"), 
             fileEncoding = "UTF-8",include.rownames = FALSE, include.colnames = FALSE)

