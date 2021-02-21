twinflex <- function(acevars, data, zyg, sep, covvars=NULL, covariance = TRUE, ordinal = NULL, optimizer = NULL, tryHard = FALSE, type = "chol", modACEuniv = NULL, modACEbiv = NULL, modBeta = NULL) {

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
  
# Check if moderation at all
if (is.null(modBeta) & is.null(modACEuniv) & is.null(modACEbiv)) {
  moderation <- FALSE
} else {
  moderation <- TRUE
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

varswideonly <- c(acevarswide,covvarswide)
# Variables that define the dimensions of the covariance matrix of the model 
variables <- c(acevarswide,covvarsall)

###############################################################################
# Prepare moderatorvars vector and some other stuff
###############################################################################
moderatorvars <- NULL
# Check if beta paths are moderated
if (!is.null(modBeta)) {
  Betamoderation <- TRUE
  if (type != "aceb") {
    stop("You cannot moderate a beta path if there are none. You need set <type> to <aceb> to estimate beta paths!")
  }
} else {
  Betamoderation <- FALSE
}
# Check if ACE paths are moderated
if (!is.null(modACEuniv) | !is.null(modACEbiv)) {
  ACEmoderation <- TRUE
  
} else {
  ACEmoderation <- FALSE
}  
print(moderation)
if (moderation == TRUE) {
  
# Split function to get the moderators out of the the strings supplied to the function
splitit <- function(input) {
if (grepl("+",input) == TRUE) {
result <- unlist(strsplit(input, split = "\\->|BY|\\+"))
} else {
result <- unlist(strsplit(input, split = "\\->|BY"))
}
}

# Define moderator vector
if (ACEmoderation == TRUE) {
  if (!is.null(modACEuniv)) {
varsACEuniv <- lapply(modACEuniv,splitit)
moderatedunivACE <- list()
moderatorunivACE <- list()
for (i in 1:length(modACEuniv)) {
varsACEuniv[[i]] <- trimws(varsACEuniv[[i]])
moderatedunivACE[[i]] <- varsACEuniv[[i]][1] # save moderated vars
moderatorunivACE[[i]] <- varsACEuniv[[i]][2:length(varsACEuniv[[i]])] # save moderators
}

modvarsACEuniv <- unique(unlist(moderatorunivACE))

  }
if (is.null(modACEuniv)) {
modvarsACEuniv <- NULL
varsACEuniv <- NULL
moderatedunivACE <- NULL
moderatorunivACE <- NULL
}
if (!is.null(modACEbiv)) {
varsACEbiv <- lapply(modACEbiv,splitit)
moderatedbivACE <- list()
moderatorbivACE <- list()
for (i in 1:length(modACEbiv)) {
varsACEbiv[[i]] <- trimws(varsACEbiv[[i]])
moderatedbivACE[[i]] <- varsACEbiv[[i]][1:2] # save moderated vars
moderatorbivACE[[i]] <- varsACEbiv[[i]][3:length(varsACEbiv[[i]])] # save moderators
}
modvarsACEbiv <- unique(unlist(moderatorbivACE))
    }
  if (is.null(modACEbiv)) {
    modvarsACEbiv <- NULL
    varsACEbiv <- NULL
  moderatedbivACE <- NULL
  moderatorbivACE <- NULL
  }
}
if (Betamoderation == TRUE) {
varsBeta <- lapply(modBeta,splitit)
moderatedBeta <- list()
moderatorBeta <- list()
for (i in 1:length(varsBeta)) {
varsBeta[[i]] <- trimws(varsBeta[[i]])
moderatedBeta[[i]] <- varsBeta[[i]][1:2] # save moderated vars
moderatorBeta[[i]] <- varsBeta[[i]][3:length(varsBeta[[i]])] # save moderators
}
varsBeta
moderatedBeta
moderatorBeta
modvarsBeta <- unique(unlist(moderatorBeta))
modvarsBeta
}
if (Betamoderation == FALSE) {
  modvarsBeta <- NULL
}
modvarsuser <- unique(c(modvarsACEuniv,modvarsACEbiv,modvarsBeta))
if (length(modvarsuser) < 5) {
  modvarsuser[(length(modvarsuser)+1):5] <- "NA"
}

modvarsmachine <- paste0("Mod",1:length(modvarsuser))


modvarsmachine
modvarsBetauser <- unique(modvarsBeta)
modvarsACEuser <- unique(c(modvarsACEuniv,modvarsACEbiv))

modvarsACEmachine <- modvarsmachine[unique(grep(paste(modvarsACEuser,collapse="|"),modvarsuser))]
modvarsBetamachine <- modvarsmachine[unique (grep(paste(modvarsBetauser,collapse="|"),modvarsuser))]


modvarslegend <- data.frame(cbind(modvarsmachine,modvarsuser))

if (ACEmoderation == TRUE) {
 modvarsACElegend <- data.frame(cbind(modvarsACEmachine,modvarsACEuser))
modvarslegend <- merge(modvarslegend, as.data.frame(modvarsACElegend), by.x = "modvarsmachine", by.y = "modvarsACEmachine", all = TRUE)
}

if (Betamoderation == TRUE) {
modvarsBetalegend <- data.frame(cbind(modvarsBetamachine,modvarsBetauser))
modvarslegend <- merge(modvarslegend, as.data.frame(modvarsBetalegend), by.x = "modvarsmachine", by.y = "modvarsBetamachine", all = TRUE)
}


aceunivmat <- NULL
if (!is.null(modACEuniv)) {
for (j in 1:length(moderatorunivACE)) {
for (i in 1:length(modvarsuser)) {
if (modvarsuser[i] %in% moderatorunivACE[[j]]) {
    newpair <- c(modvarsuser[i],paste(moderatedunivACE[[j]]))
    if (newpair[1] %in% aceunivmat[,1]) { # wenn Moderator schon in der Matrix vorhanden
    posmatch <- which(grepl(newpair[1],aceunivmat[,1]))
    new <- paste(c(aceunivmat[posmatch,2],newpair[2]),collapse=" ")
    aceunivmat[posmatch,2] <- new
      } else {
    aceunivmat <- rbind(aceunivmat,newpair)
  }
  }
  }
}
    colnames(aceunivmat) <- c("modvar","avuniv")
    rownames(aceunivmat) <- NULL
    aceunivmat <- as.data.frame(aceunivmat)
  }


acebivmat <- NULL
if (!is.null(modACEbiv)) {
for (j in 1:length(moderatorbivACE)) {
for (i in 1:length(modvarsuser)) {
if (modvarsuser[i] %in% moderatorbivACE[[j]]) {
    newpair <- c(modvarsuser[i],paste(moderatedbivACE[[j]][2]))
    if (newpair[1] %in% acebivmat[,1]) { # wenn Moderator schon in der Matrix vorhanden
    posmatch <- which(grepl(newpair[1],acebivmat[,1]))
    new <- paste(c(acebivmat[posmatch,2],newpair[2]),collapse=" ")
    acebivmat[posmatch,2] <- new
      } else {
    acebivmat <- rbind(acebivmat,newpair)
  }
  }
  }
}
    colnames(acebivmat) <- c("modvar","avbiv")
    rownames(acebivmat) <- NULL
    acebivmat <- as.data.frame(acebivmat)
  }

betamat <- NULL
if (!is.null(modBeta)) {
for (j in 1:length(moderatorBeta)) {
for (i in 1:length(modvarsuser)) {
if (modvarsuser[i] %in% moderatorBeta[[j]]) {
    newpair <- c(modvarsuser[i],paste(moderatedBeta[[j]][2]))
    if (newpair[1] %in% betamat[,1]) { # wenn Moderator schon in der Matrix vorhanden
    posmatch <- which(grepl(newpair[1],betamat[,1]))
    new <- paste(c(betamat[posmatch,2],newpair[2]),collapse=" ")
    betamat[posmatch,2] <- new
      } else {
    betamat <- rbind(betamat,newpair)
  }
  }
  }
}
    colnames(betamat) <- c("modvar","avbeta")
    rownames(betamat) <- NULL
    betamat <- as.data.frame(betamat)
  }



modvarslegend
  if (!is.null(modACEuniv)) {
modvarslegend <- merge(modvarslegend, aceunivmat, by.x = "modvarsuser", by.y = "modvar", all = TRUE)
}
  if (!is.null(modACEbiv)) {
modvarslegend <- merge(modvarslegend, acebivmat, by.x = "modvarsuser", by.y = "modvar", all = TRUE)
}
if (Betamoderation == TRUE) {
modvarslegend <- merge(modvarslegend, betamat, by.x = "modvarsuser", by.y = "modvar", all = TRUE)
}
modvarslegend <- modvarslegend[,c(2,1,3:ncol(modvarslegend))]
modvarslegend <- modvarslegend[order(modvarslegend$modvarsmachine),]
modvarslegend
print(modvarslegend)
# Check if moderator variables are in data frame
moderatorvars <- modvarslegend[modvarslegend$modvarsuser != "NA","modvarsuser"]
moderatorvars1 <-    paste0(moderatorvars,sep,"1") 
moderatorvars2 <-    paste0(moderatorvars,sep,"2") 
moderatorvarswide <- c(moderatorvars1, moderatorvars2)
moderatorvarsnotwide <- unlist(lapply(moderatorvarswide, existence))
print(moderatorvarsnotwide)
moderatorvarswide <- moderatorvarswide[!moderatorvarswide %in% moderatorvarsnotwide]
if (!is.null(moderatorvarsnotwide)) {
moderatorvarsnotwide <- unique(sapply(strsplit(moderatorvarsnotwide, split = c(paste0(sep,1),paste0(sep,2)), fixed = TRUE), function(x) (x[1])))
} else {
moderatorvarsnotwide <- NULL
}
moderatorvarsnotfound <- unlist(lapply(moderatorvarsnotwide, existence))
existenceerror(moderatorvarsnotfound)
moderatorvars <- c(moderatorvarswide,moderatorvarsnotwide)
print(moderatorvars)
# Check if there are matches between moderator variables and other variables
matchwithvars <- grep(paste(variables,collapse="|"),moderatorvars)
if (length(matchwithvars) > 0) {
moderatorvars <- moderatorvars[-matchwithvars] 
}
# the final vector with all unique moderator variables!
moderatorvars
}

###############################################################################
###############################################################################


# Check for within-variance of wide-formatted variables
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

usevariables <- c(variables,moderatorvars,zyg)
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
if (moderation == TRUE) {
mzData <- usedata %>% filter(!!as.symbol(zyg)==1) %>% select(variables,moderatorvars)
dzData <- usedata %>% filter(!!as.symbol(zyg)==2) %>% select(variables,moderatorvars)
}
if (moderation == FALSE) {
mzData <- usedata %>% filter(!!as.symbol(zyg)==1) %>% select(variables)
dzData <- usedata %>% filter(!!as.symbol(zyg)==2) %>% select(variables)
}
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
if (covariance == TRUE) {
c <- length(covvarsall) # Control variables 
}
if (covariance == FALSE) {
  c <- 0
}
l <- 3*nv*2
t <- m+l+c

###############################################################################
###############################################################################
# Matrix A: The matrix of the path coefficients
###############################################################################
###############################################################################

###############################################################################
# Matrix "b": The matrix of the phenotypic effects
###############################################################################
mat <- matrix(0.3,nrow = nv,ncol = nv)
freepathB <- lower.tri(mat)
mat[upper.tri(mat, diag = TRUE)] <- 0
valuespathB <- mat
nvstring <- as.character(1:nv)
pathBlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("b",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathBlabel[upper.tri(pathBlabel, diag = TRUE)] <- NA
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

###############################################################################
# Matrix "pZ": Just some zeros to fill some space in the matrix
###############################################################################
pathZ <- mxMatrix(type = "Zero", nrow = nv, ncol = nv, name = "pZ")

###############################################################################
# Matrix "pCov": The matrix with the effects of the covariates
  # Condition: Covariates specified in covariance matrix (covariance = TRUE)
###############################################################################

# Some helper functions for the labeling process etc.
if (!is.null(covvars)) {
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

if (covariance == TRUE) {
pathCov <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = pathCovfree,
                            values = pathCovvalue,
                            labels = pathCovlabel,
                            name = "pCov")
}
}

###############################################################################
# Matrices "a","c","e": The matrices with the unmoderated ACE effects
###############################################################################
mat <- matrix(0.3,nrow = nv,ncol = nv)
freepathAC <- lower.tri(mat, diag = TRUE)
mat[upper.tri(mat, diag = FALSE)] <- 0
valuespathAC <- mat
mat[lower.tri(mat, diag = FALSE)] <- 0
valuespathE <- mat
freepathE <- valuespathE == .3
nvstring <- as.character(1:nv)
pathAlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("a",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathAlabel[upper.tri(pathAlabel, diag = FALSE)] <- NA
pathClabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("c",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathClabel[upper.tri(pathClabel, diag = FALSE)] <- NA
pathElabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("e",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathElabel[upper.tri(pathClabel, diag = FALSE)] <- NA
pathACElb <- diag(0.0001,nv,nv)
pathACElb[pathACElb == 0] <- NA

# ACE paths with constraints if there is no bivariate moderation of ace and beta paths 
if (moderation == FALSE | (moderation == TRUE & is.null(modACEbiv) & Betamoderation == FALSE)) {
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
}

# ACE paths without constraints if there is a bivariate moderation of ace or beta paths 
if (moderation == TRUE & (!is.null(modACEbiv) | Betamoderation == TRUE)) {
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
}

###############################################################################
# Matrix "Bottom": The matrix to fill the bottom with 0 since it is not used
###############################################################################
pathBottom <- mxMatrix(type = "Zero", nrow = l+c, ncol = t, name = "Bottom")

###############################################################################
# Construct the A Matrix with mxAlgebra
###############################################################################

# Integrate the "pCov"-Matrix only if there are covariates and they should be used
# in the covariance matrix
if (!is.null(covvars) & covariance == TRUE) {
pathMan <- mxAlgebra(expression = cbind(rbind(cbind(b,pZ),
                                              cbind(pZ,b)),pCov), name = "pM")
}

# Do not integrate the "pCov"-Matrix if there are no covariates or when the 
# covariates should enter the mean vector
#if (Betamoderation == FALSE) {
if (is.null(covvars) | covariance == FALSE) {
pathMan <- mxAlgebra(expression = rbind(cbind(b,pZ),
                                              cbind(pZ,b)), name = "pM")    
}
#}
# All the ACE paths
#if (ACEmoderation == FALSE) {
pathACE <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACE")
#}
###############################################################################
# Moderation of the paths
###############################################################################
# The end product of this condition is a moderated 
  # "pM" (if there is a moderation of the betas)
  # "pCov" (if there is a moderation of the ACE vars)
if (moderation == TRUE) {
###############################################################################
# Moderation of ACE paths 
###############################################################################
if (ACEmoderation == TRUE) {
###############################################################################
# Check if moderator long or wide formatted
###############################################################################
modmatACE <- matrix(0,nrow = nv,ncol = nv, dimnames = list(acevars,acevars))
pathModACEfree <- modmatACE != 0

modvarsACEuser1 <-    paste0(modvarsACEuser,sep,"1") 
modvarsACEuser2 <-    paste0(modvarsACEuser,sep,"2") 
modvarsACEuserwide <- c(modvarsACEuser1, modvarsACEuser2)

modvarsACEusernotwide <- unlist(lapply(modvarsACEuserwide, existence))
modACElong <- rep(FALSE,length(modvarsACEuser))
if (!is.null(modvarsACEusernotwide)) {
modvarsACEusernotwide <- unique(sapply(strsplit(modvarsACEusernotwide, split = c(paste0(sep,1),paste0(sep,2)), fixed = TRUE), function(x) (x[1])))
for (i in 1:length(modvarsACEuser)) {
  for (j in 1:length(modvarsACEusernotwide)) {
    if (modvarsACEuser[i]==modvarsACEusernotwide[j]) {
      modACElong[i] <- TRUE
    }
  }
}
} else {
modvarsACEusernotwide <- NULL
}

print(modvarsACEuser)
print(modACElong)



###############################################################################
# Create matrix of interaction effects
###############################################################################

# Change names of Moderators
if (!is.null(modACEuniv)) {
  moderatormachineunivACE <- moderatorunivACE
for (i in 1:length(modvarsACEuser)) {
moderatormachineunivACE <- lapply(moderatormachineunivACE,gsub, pattern = modvarsACEuser[i], replacement = modvarsACEmachine[i])


}
  for (i in 1:length(varsACEuniv)) {
  varsACEuniv[[i]] <- c(varsACEuniv[[i]][1],moderatormachineunivACE[[i]])
}
}
if (is.null(modACEuniv)) {
  varsACEuniv <- NULL

}


# Change names of Moderators
if (!is.null(modACEbiv)) {
  moderatormachinebivACE <- moderatorbivACE
for (i in 1:length(modvarsACEuser)) {
moderatormachinebivACE <- lapply(moderatormachinebivACE,gsub, pattern = modvarsACEuser[i], replacement = modvarsACEmachine[i])


}
  for (i in 1:length(varsACEbiv)) {
  varsACEbiv[[i]] <- c(varsACEbiv[[i]][c(1,2)],moderatormachinebivACE[[i]])
}
}
if (is.null(modACEbiv)) {
  varsACEbiv <- NULL

}
print(varsACEuniv)
print(varsACEbiv)

# create list that indexes free and fixed interaction effects
freeModACE <- list(pathModACEfree,pathModACEfree,pathModACEfree,pathModACEfree,pathModACEfree)
names(freeModACE) <- modvarsmachine
    # Univariate ACE interaction effects
if (!is.null(modACEuniv)) {
for (j in 1:length(modvarsACEmachine)) {
for (i in varsACEuniv) {
if (modvarsACEmachine[j] %in% i) {
index <- modvarsACEmachine[j]
#diag(freeModACE[[index]]) <- TRUE
freeModACE[[index]][as.vector(i)[1],as.vector(i)[1]] <- TRUE
}
}
}
} 

  # Bivariate ACE interaction effects
if (!is.null(modACEbiv)) {
for (j in 1:length(modvarsACEmachine)) {
for (i in varsACEbiv) {
if (modvarsACEmachine[j] %in% i) {
index <- modvarsACEmachine[j]
print(as.vector(i)[2])
if (!is.null(freeModACE[[index]])) {
freeModACE[[index]][as.vector(i)[2],as.vector(i)[1]] <- TRUE
}
if (is.null(freeModACE[[index]])) {
freeModACE[[index]] <-  freevector
freeModACE[[index]][as.vector(i)[2],as.vector(i)[1]] <- TRUE
}
}
}
}
}


# create list that stores labels of interaction effects
labelACE <- list()
nvstring <- as.character(1:nv)
for (j in 1:length(freeModACE)) {
for (i in c("a","c","e")) {
pathModACElabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("Mod",j,"b",i,x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathModACElabel[upper.tri(pathModACElabel, diag = FALSE)] <- NA
labelACE[[paste0("Mod",j,i)]] <-pathModACElabel
}
}

# create list that stores matrices of interaction effects
pathModACEStore <- list()
for (j in names(labelACE)) {
for (i in names(freeModACE)) {
if (grepl(i,j) == TRUE) {
  index <- paste0("path",j)
  name <- paste0("p",j)
matrixstore <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = labelACE[[j]],
                       free = freeModACE[[i]],
                       values = 0,
                       name = name)
pathModACEStore[[index]] <- matrixstore
}
}
}

## differentiate between long and wide formatted moderators
modACElong
# create list with matrices of moderators of ACE paths
defModACE <- list()
  for (i in 1:length(modvarsmachine)) {
name1 <- paste0("d",modvarsmachine[i],"1","ACE")
name2 <- paste0("d",modvarsmachine[i],"2","ACE")
index1 <- paste0("def",modvarsmachine[i],"1","ACE")
index2 <- paste0("def",modvarsmachine[i],"2","ACE")
defModACE[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=NA, name= name1)
defModACE[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=NA, name= name2)  
}

for (i in 1:length(modvarsACEmachine)) { # needed: new vector with same length as modvarsACEmachine, indicating whether moderator is long formatted or not (TRUE/FALSE)
if (modACElong[i]==TRUE) { # for long formatted moderators
name1 <- paste0("d",modvarsACEmachine[i],"1","ACE")
name2 <- paste0("d",modvarsACEmachine[i],"2","ACE")
index1 <- paste0("def",modvarsACEmachine[i],"1","ACE")
index2 <- paste0("def",modvarsACEmachine[i],"2","ACE")
label <- paste0("data.",modvarsACEuser[i]) 
defModACE[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label, name= name1)
defModACE[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label, name= name2)
  }
if (modACElong[i]==FALSE) { # for wide formatted moderators
name1 <- paste0("d",modvarsACEmachine[i],"1","ACE")
name2 <- paste0("d",modvarsACEmachine[i],"2","ACE")
index1 <- paste0("def",modvarsACEmachine[i],"1","ACE")
index2 <- paste0("def",modvarsACEmachine[i],"2","ACE")
label1 <- paste0("data.",modvarsACEuser[i],sep,"1")
label2 <- paste0("data.",modvarsACEuser[i],sep,"2")
defModACE[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label1, name= name1)
defModACE[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label2, name= name2)  
  }
}

# Save matrices with ACE moderators as R objects
def <- NULL
for (i in names(defModACE)) {
assign(i, defModACE[[i]])
}

# Save matrices with ACE interaction effects as R objects
path <- NULL
for (i in names(pathModACEStore)) {
assign(i, pathModACEStore[[i]])
}


aModFull1 <- mxAlgebra(expression = a + pMod1a*dMod11ACE + pMod2a*dMod21ACE + pMod3a*dMod31ACE + pMod4a*dMod41ACE + pMod5a*dMod51ACE, name = "aMod1")
aModFull2 <- mxAlgebra(expression = a + pMod1a*dMod12ACE + pMod2a*dMod22ACE + pMod3a*dMod32ACE + pMod4a*dMod42ACE + pMod5a*dMod52ACE, name = "aMod2")
cModFull1 <- mxAlgebra(expression = c + pMod1c*dMod11ACE + pMod2c*dMod21ACE + pMod3c*dMod31ACE + pMod4c*dMod41ACE + pMod5a*dMod51ACE, name = "cMod1")
cModFull2 <- mxAlgebra(expression = c + pMod1c*dMod12ACE + pMod2c*dMod22ACE + pMod3c*dMod32ACE + pMod4c*dMod42ACE + pMod5a*dMod52ACE, name = "cMod2")
eModFull1 <- mxAlgebra(expression = e + pMod1e*dMod11ACE + pMod2e*dMod21ACE + pMod3e*dMod31ACE + pMod4e*dMod41ACE + pMod5a*dMod51ACE, name = "eMod1")
eModFull2 <- mxAlgebra(expression = e + pMod1e*dMod12ACE + pMod2e*dMod22ACE + pMod3e*dMod32ACE + pMod4e*dMod42ACE + pMod5a*dMod52ACE, name = "eMod2") 

defacepars <- c(defMod11ACE,defMod12ACE,defMod21ACE,defMod22ACE,defMod31ACE,defMod32ACE,defMod41ACE,defMod42ACE,defMod51ACE,defMod52ACE)
modacepars <- c(pathMod1a,pathMod2a,pathMod3a,pathMod4a,pathMod5a,pathMod1c,pathMod2c,pathMod3c,pathMod4c,pathMod5c,pathMod1e,pathMod2e,pathMod3e,pathMod4e,pathMod5e,
                aModFull1,cModFull1,eModFull1,aModFull2,cModFull2,eModFull2)

pathACE <- mxAlgebra(expression = rbind(cbind(aMod1,cMod1,eMod1,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,aMod2,cMod2,eMod2)), name = "pACE")
pathACEnoMod <- mxAlgebra(expression = rbind(cbind(a,c,e,pZ,pZ,pZ),
                                        cbind(pZ,pZ,pZ,a,c,e)), name = "pACEnM")
}


###############################################################################
# Moderation of Beta paths
###############################################################################
if (Betamoderation == TRUE) {

###############################################################################
# Check if moderator long or wide formatted
###############################################################################

modmatBeta <- matrix(0,nrow = nv,ncol = nv, dimnames = list(acevars,acevars))
pathModBetafree <- modmatBeta!=0

modvarsBetauser1 <-    paste0(modvarsBetauser,sep,"1") # Covariates twin 1
modvarsBetauser2 <-    paste0(modvarsBetauser,sep,"2") # Covariates twin 2
modvarsBetauserwide <- c(modvarsBetauser1, modvarsBetauser2)

modvarsBetausernotwide <- unlist(lapply(modvarsBetauserwide, existence))
modvarsBetausernotwide <- unique(sapply(strsplit(modvarsBetausernotwide, split = c(paste0(sep,1),paste0(sep,2)), fixed = TRUE), function(x) (x[1])))

modBetalong <- rep(FALSE,length(modvarsBetauser))

for (i in 1:length(modvarsBetauser)) {
  for (j in 1:length(modvarsBetausernotwide)) {
    if (modvarsBetauser[i]==modvarsBetausernotwide[j]) {
      modBetalong[i] <- TRUE
    }
  }
}
modvarsBetauser
modBetalong


###############################################################################
# Create matrix of interaction effects
###############################################################################

# Change names of Moderators
for (i in 1:length(modvarsBetauser)) {
varsBeta <- lapply(varsBeta,gsub, pattern = modvarsBetauser[i], replacement = modvarsBetamachine[i])
}
print(varsBeta)

# create list that indexes free and fixed interaction effects
freeModBeta <- list(pathModBetafree,pathModBetafree,pathModBetafree,pathModBetafree,pathModBetafree)
names(freeModBeta) <- modvarsmachine
for (j in 1:length(modvarsBetamachine)) {
freevector <- pathModBetafree
for (i in varsBeta) {
if (modvarsBetamachine[j] %in% i) {
index <- modvarsBetamachine[j]
print(index)
freeModBeta[[index]][as.vector(i)[2],as.vector(i)[1]] <- TRUE
}
}
}

# create list that stores labels of interaction effects
labelBeta <- list()
nvstring <- as.character(1:nv)
for (j in 1:length(freeModBeta)) {
pathBlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("Mod",j,"b",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathBlabel[upper.tri(pathBlabel, diag = FALSE)] <- NA
labelBeta[[paste0("Mod",j)]] <-pathBlabel
}

# create list that stores matrices of interaction effects
pathModBetaStore <- list()
for (i in names(freeModBeta)) {
  print(freeModBeta[[i]])
  index <- paste0("path",i,"Beta")
  name <- paste0("p",i,"B")
matrixstore <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = labelBeta[[i]],
                       free = freeModBeta[[i]],
                       values = 0,
                       name = name)
pathModBetaStore[[index]] <- matrixstore
}

## differentiate between long and wide formatted moderators
modBetalong
# create list with matrices of moderators of Beta paths
defModBeta <- list()
  for (i in 1:length(modvarsmachine)) {
name1 <- paste0("d",modvarsmachine[i],"1","Beta")
name2 <- paste0("d",modvarsmachine[i],"2","Beta")
index1 <- paste0("def",modvarsmachine[i],"1","Beta")
index2 <- paste0("def",modvarsmachine[i],"2","Beta")
defModBeta[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=NA, name= name1)
defModBeta[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=NA, name= name2)  
}
for (i in 1:length(modvarsBetamachine)) { # needed: new vector with same length as modvarsBetamachine, indicating whether moderator is long formatted or not (TRUE/FALSE)
  if (modBetalong[i]==TRUE) { # for long formatted moderators
name1 <- paste0("d",modvarsBetamachine[i],"1","B")
name2 <- paste0("d",modvarsBetamachine[i],"2","B")
index1 <- paste0("def",modvarsBetamachine[i],"1","Beta")
index2 <- paste0("def",modvarsBetamachine[i],"2","Beta")
label <- paste0("data.",modvarsBetauser[i]) 
defModBeta[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label, name= name1)
defModBeta[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label, name= name2)
  }
  if (modBetalong[i]==FALSE) { # for wide formatted moderators
name1 <- paste0("d",modvarsBetamachine[i],"1","B")
name2 <- paste0("d",modvarsBetamachine[i],"2","B")
index1 <- paste0("def",modvarsBetamachine[i],"1","Beta")
index2 <- paste0("def",modvarsBetamachine[i],"2","Beta")
label1 <- paste0("data.",modvarsBetauser[i],sep,"1")
label2 <- paste0("data.",modvarsBetauser[i],sep,"2")
defModBeta[[index1]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label1, name= name1)
defModBeta[[index2]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label2, name= name2)  
  }
}
# Save matrices with Beta moderators as R objects
def <- NULL
for (i in names(defModBeta)) {
assign(i, defModBeta[[i]])
}

# Save matrices with Beta interaction effects as R objects
path <- NULL
for (i in names(pathModBetaStore)) {
assign(i, pathModBetaStore[[i]])
}
  
betaModFull1 <- mxAlgebra(expression = b + pMod1B*dMod11B + pMod2B*dMod21B + pMod3B*dMod31B + pMod4B*dMod41B + pMod5B*dMod51B, name = "betaMod1")
betaModFull2 <- mxAlgebra(expression = b + pMod1B*dMod12B + pMod2B*dMod22B + pMod3B*dMod32B + pMod4B*dMod42B + pMod5B*dMod52B, name = "betaMod2")

defbetapars <- c(defMod11Beta,defMod12Beta,defMod21Beta,defMod22Beta,defMod31Beta,defMod32Beta,defMod41Beta,defMod42Beta,defMod51Beta,defMod52Beta)
modbetapars <- c(pathMod1Beta,pathMod2Beta,pathMod3Beta,pathMod4Beta,pathMod5Beta,betaModFull1,betaModFull2)

pathMan <- mxAlgebra(expression = rbind(cbind(betaMod1,pZ),
                                              cbind(pZ,betaMod2)), name = "pM")
pathMannoMod <- mxAlgebra(expression = rbind(cbind(beta,pZ),
                                              cbind(pZ,beta)), name = "pMnM")
}
}

# The A Matrix!
matA <- mxAlgebra(expression = rbind(cbind(pM,pACE),
                                     Bottom),
                  name = "A")

if (moderation == TRUE) {
if (ACEmoderation == TRUE & Betamoderation == TRUE) {
matAnoMod <- mxAlgebra(expression = rbind(cbind(pMnM,pACEnM),
                                          Bottom),name = "AnM")
modbinconstrpars <- c(pathACEnoMod,pathMannoMod,matAnoMod)
} else if (ACEmoderation == TRUE & Betamoderation == FALSE) {
matAnoMod <- mxAlgebra(expression = rbind(cbind(pM,pACEnM),
                                          Bottom),name = "AnM")
modbinconstrpars <- c(pathACEnoMod,matAnoMod)
} else if (ACEmoderation == FALSE & Betamoderation == TRUE) {
matAnoMod <- mxAlgebra(expression = rbind(cbind(pMnM,pACE),
                                          Bottom),name = "AnM")
modbinconstrpars <- c(pathMannoMod,matAnoMod)
}
}

###############################################################################
###############################################################################
# Matrix S: The matrix of the covariances
###############################################################################
###############################################################################

# some helper functions for the labeling process etc.
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
return(labelmat)
}

# the starting values
svS <- unname(as.matrix(var(data[,covvarsall], use = "na.or.complete"))) # S matrix starting values for non-decomposed covariates

# make sure the starting value matrix is symmetric --> same starting value for the same parameter
if (!isSymmetric(svS)) {
svS[upper.tri(svS)] <- t(svS)[upper.tri(svS)]
}

# If there are Covariates and they should go into the covariance matrix
# specify all covariances with covariates
if (!is.null(covvars) & covariance == TRUE) {
  # The Covariances between the covariates
covCovariates <- mxMatrix(type = "Symm", nrow = c, ncol = c, byrow = FALSE,
                 values = svS,
                 lbound = lowerboundcovmat(c), 
                 labels = labelcovmat(covvarsall),
                 free = TRUE,
                 name = "cCov")
# The Covariances between the Covariates and the ACEvars -> all zero
covManCov <- mxMatrix(type = "Zero", nrow = m, ncol = c, name = "cManCov")
# The Covariances between the ACEvars -> all zero
covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
# The Covariances between the manifest, covariates and latent ACE-Vars -> all zero
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")
matSMan <- mxAlgebra(expression = rbind(cbind(cMan,cManCov), 
                                        cbind(t(cManCov),cCov),
                                        t(cManCovACE)), name = "matSM")
}

# If there are no Covariates or they should go into the mean,
# don't specify any covariances since they are not part of the covariance matrix
if (is.null(covvars) | covariance == FALSE) {
covCovariates <- NULL
covManCov <- NULL
# The "cMan" and "cManCovACE" Matrices remain unchanged
covMan <- mxMatrix(type = "Zero", nrow = m, ncol = m, name = "cMan")
covManCovACE <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "cManCovACE")
matSMan <- mxAlgebra(expression = rbind(cMan,t(cManCovACE)), name = "matSM")
}

# The Variances of the latent ACE factors
covV <- mxMatrix(type = "Iden", nrow = l/2, ncol = l/2, name = "V")
# The Covariances of the latent ACE factors for MZ twins
covCMZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(1,(nv*2)),rep(0,nv)),
                   name = "CMZ")
# The Covariances of the latent ACE factors for DZ twins
covCDZ <- mxMatrix(type = "Diag", nrow = l/2, ncol = l/2,
                   values = c(rep(.5,nv),rep(1,nv),rep(0,nv)),
                   name = "CDZ")
# The Covariance Matrix for MZ twins
matSMZ <- mxAlgebra(expression = cbind(matSM,rbind(cManCovACE,
                                                   cbind(V,CMZ),
                                                   cbind(CMZ,V))),
                    name = "SMZ")

# The Covariance Matrix for DZ twins
matSDZ <- mxAlgebra(expression = cbind(matSM,rbind(cManCovACE,
                                                   cbind(V,CDZ),
                                                   cbind(CDZ,V))),
                    name = "SDZ")

# Build the Filter and Identity Matrix (necessary for the Construction of the 
# expected covariance matrix via RAM matrices)
filterI <- mxMatrix(type = "Iden", nrow = m+c, ncol = m+c, name = "FI")
filterZ <- mxMatrix(type = "Zero", nrow = m+c, ncol = l, name = "FZ")
matF <- mxAlgebra(expression = cbind(FI,FZ), name = "Filter")
matI <- mxMatrix(type = "Iden", nrow = t, ncol = t, name = "I")

# The expected covariance matrix for the MZ twins
covMZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SMZ%*%t(solve(I-A))%*%t(Filter), name = "expCovMZ")
# The expected covariance matrix for the DZ twins
covDZ <- mxAlgebra(expression = Filter%*%solve(I-A)%*%SDZ%*%t(solve(I-A))%*%t(Filter), name = "expCovDZ")

if (moderation == TRUE) {
covMZnoMod <- mxAlgebra(expression = Filter%*%solve(I-AnM)%*%SMZ%*%t(solve(I-AnM))%*%t(Filter), name = "expCovMZnMod")
modbinconstrpars <- c(modbinconstrpars,covMZnoMod)
}

###############################################################################
###############################################################################
# Mean Matrix
###############################################################################
###############################################################################

###################################
### MATRIX OF UNMODERATED MEANS ###
###################################

# If covariates in covariance matrix: Define matrix of unmoderated Means ("M") and matrix of expected means ("expMean")
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
# If covariates in means matrix: Define matrix of unmoderated Means ("M")
if (covariance == FALSE) {
meanmanifestlabel <- c(paste0("mean_",unlist(sapply(strsplit(acevarswide, split=sep, fixed=TRUE), function(x) (x[1])))))
meanacelabel <- paste0("mean",c(paste0("A_",acevars1),paste0("C_",acevars1),paste0("E_",acevars1),paste0("A_",acevars2),paste0("C_",acevars2),paste0("E_",acevars2)))
meanlabel <- c(meanmanifestlabel,meanacelabel)
matM <- mxMatrix(type = "Full", nrow = m+l, ncol = 1, 
                 free = c(rep(TRUE,m),rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M") 
}

# If there are binary ACE vars
if (!is.null(ordinal)) {
  if (2 %in% ordinallength) {
binarytrue <- nTh == 1
binaryvar <- ordinalwide[binarytrue]
binaryace <- unlist(lapply(acevarswide,checkcorrespondence, check = binaryvar))
print(binaryace)
binaryrest <- rep(FALSE,(length(variables)-length(acevarswide)))
      # If there are binary vars + covariates in covariance matrix: Define matrix of unmoderated Means ("M") and matrix of expected means ("expMean")
if (covariance == TRUE) { # covariates in model with binary vars with covariance = TRUE
binaryflag <- c(binaryace,binaryrest)
svmean_manifests <- binaryflag == FALSE
matM <- mxMatrix(type = "Full", nrow = t, ncol = 1, 
                 free = c(svmean_manifests,rep(FALSE,l)), 
                 labels = c(meanmanifestlabel,meanacelabel),
                 values = c(svmean,rep(0,l)), 
                 name = "M")
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%M), name = "expMean")
}
      # If there are binary vars + covariates in means matrix: Define matrix of unmoderated Means ("M") 
if (covariance == FALSE) { 
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
}
# If covariates in means matrix: Until now we only defined the matrix of unmoderated Means ("M"), now we define 
      # the effects of the covariates in the mean matrix ("pCov")
      # the covariates ("defMCov")
      # the regressions of the covariates without the latents ("effMCov")
      # the regressions ot the covariates with the latents fixed to zero and unmoderated ("effMCovFull")
      # the moderated means matrix ("modMCov")
if (covariance == FALSE) { 
  c <- length(covvarsall)
# matrix with effect sizes of moderation of the means stored in M
pathCov <- mxMatrix(type = "Full", nrow = c, ncol = ntv, byrow = FALSE,
                            free = t(pathCovfree),
                            values = t(pathCovvalue),
                            labels = t(pathCovlabel),
                            name = "pCov")

# Matrix of definition variables for mean moderation
labeldef <- paste0("data",".",covvarsall)
defCovMean      <- mxMatrix( type="Full", nrow=1, ncol=c, free=FALSE, labels=labeldef, name="defMCov" )

# Matrix effects on means
effMeanCov <- mxAlgebra(expression = defMCov%*%pCov, name = "effMCov")

# Vector of latent variables (all set to zero) just need to concatenate them to the manifests to get the dimensions right
latentmeans <- mxMatrix(type = "Full", nrow = 1, ncol = l, name = "lmeans")

# Concatenate the manifest with the latent means vector
effMeanCovFull <- mxAlgebra(expression = cbind(effMCov,lmeans), name = "effMCovFull") # this matrix will be used if there is a moderation AND covariates in the means matrix -> modM and expMean will be overwritten if there is a moderation

# Matrix of moderated means
modMeanCov <- mxAlgebra(expression =M+t(effMCovFull), name = "modM")

# Matrix of expected means
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%modM), name = "expMean")
}

# add here the construction of the expected means matrix if there is a moderation 
  # it does not matter that there are already expected means defined, here they will
  # be overruled
###############################################################################
# Moderation of the mean
###############################################################################
if (moderation == TRUE) {
# Moderation of the mean matrix
modvarsuser <- modvarsuser[modvarsuser!="NA"]
modvarsmachine
modvarslegend

if (ACEmoderation == TRUE) {
modvarsACElong <- data.frame(cbind(modvarsACEmachine,modACElong))
modvarslegend <- merge(modvarslegend, as.data.frame(modvarsACElong), by.x = "modvarsmachine", by.y = "modvarsACEmachine", all = TRUE)
}

if (Betamoderation == TRUE) {
modvarsBetalong <- data.frame(cbind(modvarsBetamachine,modBetalong))
modvarslegend <- merge(modvarslegend, as.data.frame(modvarsBetalong), by.x = "modvarsmachine", by.y = "modvarsBetamachine", all = TRUE)
}
modvarslegend

modmeanlabel <- NULL
for (i in 1:length(modvarsmachine)) {
  effect_t1 <- c(paste0("bMean",modvarsmachine[i],1:nv),rep(NA,nv))
  effect_t2 <- c(rep(NA,nv),paste0("bMean",modvarsmachine[i],1:nv))
  modmeanlabel <- rbind(modmeanlabel,effect_t1,effect_t2)
}
colnames(modmeanlabel) <- rep(acevars,2)

modmeanfree <- !is.na(modmeanlabel)
modmeanfree[modmeanfree = TRUE] <- FALSE

if (nv == 1) {
  modmeanfree <- modmeanfree[1:5,nv]
  modmeanfree <- matrix(modmeanfree, nrow = 5, ncol = nv, dimnames = list(modvarslegend$modvarsmachine,acevars))
}
if (nv > 1) {
    modmeanfree <- modmeanfree[1:5,1:nv]
    row.names(modmeanfree) <- modvarslegend$modvarsmachine
}

for (i in 1:(dim(modvarslegend)[1])) {
for (j in acevars) {
  if (!is.null(modACEuniv)) {
  if (grepl(j,modvarslegend[i,"avuniv"])) {
    modmeanfree[paste0("Mod",i),j] <- TRUE
    row.names(modmeanfree) <- modvarslegend$modvarsuser
  if (rownames(modmeanfree)[i] %in% acevars) {  
#  if (rownames(modmeanfree)[i] == j) {
    modmeanfree[i,] <- FALSE
  }
    row.names(modmeanfree) <- modvarslegend$modvarsmachine
  }
  }
  if (!is.null(modACEbiv)) {
  if (grepl(j,modvarslegend[i,"avbiv"])) {
    print(modvarslegend[i,"avbiv"])
    modmeanfree[paste0("Mod",i),j] <- TRUE
    row.names(modmeanfree) <- modvarslegend$modvarsuser
  if (rownames(modmeanfree)[i] %in% acevars) {  
#  if (rownames(modmeanfree)[i] == j) {
    modmeanfree[i,] <- FALSE
  }
    row.names(modmeanfree) <- modvarslegend$modvarsmachine
  }
    }
if (Betamoderation == TRUE) {
  if (grepl(j,modvarslegend[i,"avbeta"])) {
    modmeanfree[paste0("Mod",i),j] <- TRUE
        row.names(modmeanfree) <- modvarslegend$modvarsuser
  if (rownames(modmeanfree)[i] %in% acevars) {  
#  if (rownames(modmeanfree)[i] == j) {
    modmeanfree[i,] <- FALSE
  }
        row.names(modmeanfree) <- modvarslegend$modvarsmachine
  }
}
}
}
modmeanfree
modvarslegend

modmeanfree <- cbind(modmeanfree,modmeanfree)
modmeanfree <- modmeanfree[rep(1:nrow(modmeanfree), times = c(2,2,2,2,2)),]

modmeanfree[c(1,3,5,7,9),(length(acevars)+1):(length(acevars)*2)] <- FALSE
modmeanfree[c(2,4,6,8,10),1:(length(acevars))] <- FALSE

modmeanval <-  modmeanfree
modmeanval[!is.null(modmeanval)] <- 0

# matrix with effect sizes of moderation of the means stored in M
pathModMean <- mxMatrix(type = "Full", nrow = 5*2, ncol = ntv, byrow = FALSE, # 5 as maximum number of moderators (at the time)
                            free = (modmeanfree),
                            values = (modmeanval),
                            labels = (modmeanlabel),
                            name = "pMMod",
                        dimnames = list(NULL,NULL))
pathModMean
labeldefMModean <- NULL
modvarslegend[is.na(modvarslegend)] <- "NA"
# Matrix of definition variables for mean moderation
for (i in 1:length(modvarsmachine)) {
  print(i)
if (ACEmoderation == TRUE) {
if (modvarslegend$modACElong[i]== "FALSE") {  # wide-formatted moderators
l1 <- paste0("data.",modvarslegend$modvarsuser[i],sep,"1")
l2 <- paste0("data.",modvarslegend$modvarsuser[i],sep,"2")
}  else if (modvarslegend$modACElong[i]== "TRUE") {
l1 <- paste0("data.",modvarslegend$modvarsuser[i])
l2 <- paste0("data.",modvarslegend$modvarsuser[i])
} else if (modvarslegend$modACElong[i]== "NA") {
    l1 <- NA
  l2 <- NA
}
} 
if (Betamoderation == TRUE) {
if (modvarslegend$modBetalong[i] == "FALSE") {
l1 <- paste0("data.",modvarslegend$modvarsuser[i],sep,"1")
l2 <- paste0("data.",modvarslegend$modvarsuser[i],sep,"2")
} else if (modvarslegend$modBetalong[i]== "TRUE") {
l1 <- paste0("data.",modvarslegend$modvarsuser[i])
l2 <- paste0("data.",modvarslegend$modvarsuser[i])
} else if (modvarslegend$modBetalong[i]== "NA") {
  l1 <- NA
  l2 <- NA
}
}
  labeldefMModean <- rbind(labeldefMModean,l1,l2)
  row.names(labeldefMModean) <- NULL
  
}
labeldefMModean

defMModean      <- mxMatrix( type="Full", nrow=1, ncol=10, free=FALSE, labels=labeldefMModean, name="defMMod" )

# Matrix effects on means
effMModean <- mxAlgebra(expression = defMMod%*%pMMod, name = "effMMod")

# Vector of latent variables (all set to zero) just need to concatenate them to the manifests to get the dimensions right
latmodmeans <- mxMatrix(type = "Full", nrow = 1, ncol = c+l, name = "lmodmeans") # theoretisch müssten Kovariate in der Kovarianzmatrix möglich sein, daher das c

# Concatenate the manifest with the latent means vector
effMModeanFull <- mxAlgebra(expression = cbind(effMMod,lmodmeans), name = "effMModFull")

# Matrix of moderated means
modMean <- mxAlgebra(expression =M+t(effMModFull), name = "modM")
if (covariance == FALSE) {
modMean <- mxAlgebra(expression =M+t(effMModFull)+t(effMCovFull), name = "modM")
  }

# Matrix of expected means
mean <- mxAlgebra(expression = t(Filter%*%solve(I-A)%*%modM), name = "expMean")

#modvarslegend <- modvarslegend[,1:2]
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
binaryflag <- c(binaryflag[1:length(acevars)],rep(FALSE,length(acevars)))
nrowfilterbinary <- sum(binaryflag)
print("HERE")
print(binaryflag)
print(nrowfilterbinary)
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
print("lala")
print(filtermatvalues)
filtermatbin <- mxMatrix(type = "Full", values = filtermatvalues, name = "fmatbin")
print(filtermatbin)
binarycov <- mxAlgebra(expression = fmatbin %*%expCovMZ %*% t(fmatbin), name = "binCov")

if (moderation == TRUE) {
binarycov <- mxAlgebra(expression = fmatbin %*%expCovMZnMod %*% t(fmatbin), name = "binCov") 
}

one <- mxMatrix(type = "Unit", nrow = nrowfilterbinary, ncol = 1, name = "Unit")
print(one)
var1 <- mxConstraint(expression = diag2vec(binCov)==Unit , name = "VConstraint1")
binary <- c(filtermatbin,binarycov,one,var1)
if (moderation == TRUE) {
binary <- c(binary,modbinconstrpars)  
}
}
}
if (covariance == FALSE) {
  variables <- c(acevarswide)
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
pars      <- list(pathB, pathA, pathC, pathE, pathZ, matA, pathMan, pathACE,
                  covCovariates, covMan, covManCov, covManCovACE, covV,matSMan,
                  filterI, filterZ, matF, matI,
                  matM, mean, pathBottom)
if (!is.null(ordinal)) {
pars <- c(pars,c(thinG,inc,threG))
}

if (!is.null(covvars)) {
if (covariance == FALSE) {
covmeanpars <- c(effMeanCov,latentmeans,effMeanCovFull,modMeanCov,pathCov)  
pars <- c(pars,covmeanpars)
}
if (covariance == TRUE) {
pars <- c(pars,pathCov)
}
}

if (moderation == TRUE) {
if (ACEmoderation == FALSE) {
  modacepars <- NULL
  defacepars <- NULL
}
  if (Betamoderation == FALSE) {
  modbetapars <- NULL
  defbetapars <- NULL
}
  modmeanpars <- c(pathModMean, effMModean,latmodmeans,effMModeanFull,modMean)
  defmodpars <- c(defMModean,defacepars,defbetapars)
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
if (covariance == FALSE & !is.null(covariance)) {
  modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, defCovMean,name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun,defCovMean,name="DZ")

}
if (moderation == TRUE) {
if (!is.null(ordinal)) {
if (2 %in% ordinallength) {
  modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, binary,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars,name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars,name="DZ")
}
}
if (covariance == FALSE & !is.null(covariance)) {
modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun, defCovMean,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars,name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun,defCovMean,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars,name="DZ")

}
modelMZ   <- mxModel(pars, covCMZ, covMZ, matSMZ, dataMZ, expMZ, fitfun,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars, name="MZ")
modelDZ   <- mxModel(pars, covCDZ, covDZ, matSDZ, dataDZ, expDZ, fitfun,defacepars,defbetapars,defmodpars,modacepars,modbetapars,modmeanpars, name="DZ")
}
multi     <- mxFitFunctionMultigroup(c("MZ","DZ"))

# overall model object
modelACE  <- mxModel("ACE",  modelMZ, modelDZ, multi)
# run model

# uncomment to control the matrices involved in the moderation part
#if (moderation == TRUE) {
#sink("example.txt")
#print(variables)
#print(modvarslegend)
#print(modmeanpars)
#print(defmodpars)
#print(pathB)
#print(pathA)
#print(pathC)
#print(pathE)
#sink()
#stop("Stop!")
#}

if (is.null(optimizer)) {
mxOption(NULL , 'Default optimizer' , 'SLSQP')
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
}
sumACE <- summary(fitACE)
if (moderation == FALSE) {
  modvarslegend <- NULL
}
functionresult <- list("ModelFitted" = fitACE, "ModelSummary" = sumACE,"ParameterLegend" = NULL, "ModeratorLegend" = modvarslegend)
return(functionresult)
}