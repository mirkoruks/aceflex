# Implement moderation 
rm(list = ls())

# load packages
library(OpenMx)
library(dplyr)
library(OpenMx)
library(xtable)

data <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
# create binary
data <- data %>% 
  mutate(schoolbin_1 = ifelse(schoolhigh_1 %in% c(1,2), 1,
                       ifelse(schoolhigh_1 %in% c(3,4), 2, NA))) %>% 
  mutate(schoolbin_2 = ifelse(schoolhigh_2 %in% c(1,2), 1,
                       ifelse(schoolhigh_2 %in% c(3,4), 2, NA))) %>% rename(zy = zyg)
summary(data)
## INTEGRATE BETA  MODERATION
# moderation only of paths from ACE components on decomposed variables and from decomposed variables on decomposed variables

# new arguments:
  # modACEuniv = NULL ## specify a string vector which univariate ACE paths you want to interact with which moderator 
    # Example: modACEuniv <- c("posbez BY schoolhigh + sex","negbez BY schoolhigh + iseiempmean")
    # Explanation: Moderate univariate ACE paths of posbez by schoolhigh and sex; moderate univariate ACE paths of negbez by schoolhigh and iseiempmean
  # modACEbiv = NULL ## if moderation = TRUE -> specify here if you want to moderate bivariate ACE paths, and if yes, of which variable pairs 
    # input element = character vector of variable relationships
      # e.g.: Both, X and Y are ACE decomposed variables and the user wants to moderate the effect of the ACE components of X on Y -> 
      # the user needs to write: modACEbiv = "X -> Y"
  # modBeta = NULL ## if moderation = TRUE -> specify here if you want to moderate a phenotypic paths between acevars, and if yes, of which variables

# necessary stuff to keep it running
acevars <- c("posbez","negbez") # three variables to be decomposed into ACE components
variables <- c("posbez_1","posbez_2","negbez_1","negbez_2")
nv <- length(acevars) # Vars per twin
ntv <- nv*2
sep = "_"
modACEuniv <- c("negbez BY iseiempmean + age","posbez BY age") #,"negbez BY schoolhigh + iseiempmean")
modACEbiv <- c("posbez -> negbez BY iseiempmean + age")
modBeta <- NULL #c("posbez->negbez BY sex","negbez -> iseiempmean BY sex")
type = "aceb"
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


###############################################################################
###############################################################################
# From here on: Code to be implemented! 

# Check if moderation at all
if (is.null(modBeta) & is.null(modACEuniv) & is.null(modACEbiv)) {
  moderation <- FALSE
} else {
  moderation <- TRUE
}

###############################################################################
# Prepare moderatorvars vector and some other stuff
###############################################################################
if (moderation == TRUE) {
# Split function to get the moderators out of the the strings supplied to the function
splitit <- function(input) {
if (grepl("+",input) == TRUE) {
result <- unlist(strsplit(input, split = "\\->|BY|\\+"))
} else {
result <- unlist(strsplit(input, split = "\\->|BY"))
}
}

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
# Moderation of the paths
###############################################################################
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
}
}


###############################################################################
# Moderation of the mean
###############################################################################

## ONLY ENTER THE MEAN IF MODERATOR NOT PART OF ACEVARS!!
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

