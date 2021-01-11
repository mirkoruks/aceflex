# Trivariate Cholesky Model
rm(list = ls())

library(dplyr)
library(OpenMx)
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

# aceflex function
aceflex <- function(acevars, data, zyg, covvars=NULL) {

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
        stop("User did not want to proceed due to (multi-)collinearity in one of the variables. Function aceflex has stopped! See you soon!")
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
print(pathB)
cat("\n\n\n... Bis hierhin läuft alles durch! Weiter geht's! :-)")
}
 
aceflex(acevars = acevars, data = data_orig,zyg = "zyg")




# Build elements to construct expected covariance matrix (RAM Notation)

# Matrix A 
pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathB,
                  values = valuespathB,
                  labels = pathBlabel,
                  name = "b")
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

# aceflex function
aceflex <- function(acevars, data, zyg, covvars=NULL) {

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
        stop("User did not want to proceed due to (multi-)collinearity in one of the variables. Function aceflex has stopped! See you soon!")
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
 
aceflex(acevars = acevars, data = data_orig,zyg = "zyg")




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

>>>>>>> 6b2c1e6ff10920cb8194250c757ae35b4510bad8
