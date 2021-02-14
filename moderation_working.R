# Implement moderation 
  library(OpenMx)
  
  ## INTEGRATE UNI AND BIVARIATE MODERATION OF ACE
  ## INTEGRATE BETA  MODERATION
# moderation only of paths from ACE components on decomposed variables and from decomposed variables on decomposed variables

# new arguments:
  # moderation = NULL ## TRUE if moderation needed
  # modACEuniv = NULL ## if moderation = TRUE -> specify here if you want to moderate univariate ACE paths, and if yes, of which variables
    # input element = character vector of variables: The univariate ACE paths of each variable specified in the vector will be moderated 
  # modACEbiv = NULL ## if moderation = TRUE -> specify here if you want to moderate bivariate ACE paths, and if yes, of which variable pairs 
    # input element = character vector of variable relationships
      # e.g.: Both, X and Y are ACE decomposed variables and the user wants to moderate the effect of the ACE components of X on Y -> 
      # the user needs to write: modACEbiv = "X -> Y"
  # modBeta = NULL ## if moderation = TRUE -> specify here if you want to moderate a phenotypic paths between acevars, and if yes, of which variables

# necessary stuff to keep it running
acevars <- c("X","Z","Y") # three variables to be decomposed into ACE components
nv <- length(acevars) # Vars per twin
modACEbiv <- c("X   -> Y","Z -> Y")
modBeta <- c("X->Y")
modACEbiv2 <- c("X   -> Y BY K","Z -> Y BY K + Mod2")
modACEbiv2
#grepl("+",modACEbiv2)
#strsplit(modACEbiv2, split = "\\->|BY|\\+")

splitit <- function(input) {
if (grepl("+",input) == TRUE) {
result <- unlist(strsplit(input, split = "\\->|BY|\\+"))
} else {
result <- unlist(strsplit(input, split = "\\->|BY"))
}
  
}
test <- lapply(modACEbiv2,splitit)

moderatedbivACE <- vector(mode = "list", length = length(modACEbiv2))
moderatorbivACE <- vector(mode = "list", length = length(modACEbiv2))
for (i in 1:length(modACEbiv2)) {
test[[i]] <- trimws(test[[i]])
moderatedbivACE[[i]] <- test[[i]][1:2] # save moderated vars
moderatorbivACE[[i]] <- test[[i]][3:length(test[[i]])] # save moderators
}
test
moderatedbivACE
moderatorbivACE
modvarsACEbiv <- unique(unlist(moderatorbivACE))
modvarsACEbiv

# create matrix of interaction effects
modmat <- matrix(0,nrow = 3,ncol = 3, dimnames = list(acevars,acevars))
modpathACEfree <- modmat!=0

freeList <- list()
count <- 1 
for (j in modvarsACEbiv) {
freevector <- modpathACEfree
for (i in test) {
if (j %in% i) {
index <- paste0("Mod",count)
freevector[as.vector(i)[2],as.vector(i)[1]] <- TRUE
freeList[[index]] <-  freevector
}
}
count <- count+1
}
View(freeList)

labelList <- list()
nvstring <- as.character(1:nv)
for (i in c("a","c","e")) {
pathModlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("b",i,x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathModlabel[upper.tri(pathAModlabel, diag = FALSE)] <- NA
labelList[[i]] <-pathModlabel
}
View(labelList)

count <- 1
pathModStore <- list()
for (i in names(freeList)) {
for (j in names(labelList)) {
  index <- paste0("pathMod",count,j)
  name <- paste0("pMod",count,j)
matrixstore <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = labelList[[j]],
                       free = freeList[[i]],
                       values = 0,
                       name = name)
pathModStore[[index]] <- matrixstore
}
count <- count+1
}
View(pathModStore)

#assign(paste0("path", i,j), pathModStore)

# create list with matrices of moderators of bivariate ACE paths
defMod <- list()
count = 1
for (i in modvarsACEbiv) {
index <- paste0("dMod",count)
label <- paste0("data.",i) 
defMod[[index]]      <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=label, name= index)
count <- count +1
}
View(defMod)
View(pathModStore)

def <- NULL
for (i in names(defMod)) {
assign(paste0("def", i), defMod[[i]])
}

path <- NULL
for (i in names(pathModStore)) {
assign(paste0("path", i), pathModStore[[i]])
}

modvars <- modvarsACEbiv
# generate moderated paths
  # for long formatted moderators -> save with wide name but long label -> e.g. aMod11 and aMod12 with label "data.Mod1", that is, they are identical!

if (length(modvars == 1)) {
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod11, name = "aMod1")
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod12, name = "aMod2")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod11, name = "cMod1")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod12, name = "cMod2")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod11, name = "eMod1")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod12, name = "eMod2")
} else if (length(modvars == 2)) {
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod11 + pMod2a*defMod21, name = "aMod1")
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod12 + pMod2a*defMod22, name = "aMod2")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod11 + pMod2c*defMod21, name = "cMod1")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod12 + pMod2c*defMod22, name = "cMod2")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod11 + pMod2e*defMod21, name = "eMod1")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod12 + pMod2e*defMod22, name = "eMod2") 
} else if (length(modvars == 3)) {
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod11 + pMod2a*defMod21 + pMod3a*defMod31, name = "aMod1")
  aModFull <- mxAlgebra(expression = a + pMod1a*defMod12 + pMod2a*defMod22 + pMod3a*defMod31, name = "aMod2")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod11 + pMod2c*defMod21 + pMod3c*defMod31, name = "cMod1")
  cModFull <- mxAlgebra(expression = c + pMod1c*defMod12 + pMod2c*defMod22 + pMod3c*defMod31, name = "cMod2")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod11 + pMod2e*defMod21 + pMod3e*defMod31, name = "eMod1")
  eModFull <- mxAlgebra(expression = e + pMod1e*defMod12 + pMod2e*defMod22 + pMod3e*defMod31, name = "eMod2") 
}































































if (!is.na(modvars[4])) {
 if (modvarswide[4]==FALSE) {
pathaMod4 <- mxAlgebra(expression = pathMod4a*defMod4, name = "aMod4")
pathcMod4 <- mxAlgebra(expression = pathMod4c*defMod4, name = "cMod4")
patheMod4 <- mxAlgebra(expression = pathMod4e*defMod4, name = "eMod4")
} else if (modvarswide[4]==TRUE) {
pathaMod41 <- mxAlgebra(expression =pathMod4a*defMod41, name = "aMod41")
pathcMod41 <- mxAlgebra(expression =pathMod4c*defMod41, name = "cMod41")
patheMod41 <- mxAlgebra(expression =pathMod4e*defMod41, name = "eMod41")  
pathaMod42 <- mxAlgebra(expression =pathMod4a*defMod42, name = "aMod42")
pathcMod42 <- mxAlgebra(expression =pathMod4c*defMod42, name = "cMod42")
patheMod42 <- mxAlgebra(expression =pathMod4e*defMod42, name = "eMod42") 
}  
}
if (!is.na(modvars[5])) {
 if (modvarswide[5]==FALSE) {
pathaMod5 <- mxAlgebra(expression = pathMod5a*defMod5, name = "aMod5")
pathcMod5 <- mxAlgebra(expression = pathMod5c*defMod5, name = "cMod5")
patheMod5 <- mxAlgebra(expression = pathMod5e*defMod5, name = "eMod5")
} else if (modvarswide[5]==TRUE) {
pathaMod51 <- mxAlgebra(expression =pathMod5a*defMod51, name = "aMod51")
pathcMod51 <- mxAlgebra(expression =pathMod5c*defMod51, name = "cMod51")
patheMod51 <- mxAlgebra(expression =pathMod5e*defMod51, name = "eMod51")  
pathaMod52 <- mxAlgebra(expression =pathMod5a*defMod52, name = "aMod52")
pathcMod52 <- mxAlgebra(expression =pathMod5c*defMod52, name = "cMod52")
patheMod52 <- mxAlgebra(expression =pathMod5e*defMod52, name = "eMod52") 
}  
}

else if (length(modvars == 1) & modvarswide[1]==TRUE) {

modpars <- c(aMod1,cMod1,eMod1,pathaMod,pathcMod,patheMod)
  
} else if (length(modvars == 2)) {
aMod1 <- mxAlgebra(expression = pathMod1a*defMod1, name = "aM1")
cMod1 <- mxAlgebra(expression = pathMod1c*defMod1, name = "cM1")
eMod1 <- mxAlgebra(expression = pathMod1e*defMod1, name = "eM1")
aMod2 <- mxAlgebra(expression = pathMod2a*defMod2, name = "aM2")
cMod2 <- mxAlgebra(expression = pathMod2c*defMod2, name = "cM2")
eMod2 <- mxAlgebra(expression = pathMod2e*defMod2, name = "eM2")
pathaMod <- mxAlgebra(expression = a +aM1+aM2, name = "aMod")
pathcMod <- mxAlgebra(expression = a +cM1+cM2, name = "cMod")
patheMod <- mxAlgebra(expression = a +eM1+eM2, name = "eMod")
} else if (length(modvars == 3)) {
aMod1 <- mxAlgebra(expression = pathMod1a*defMod1, name = "aM1")
cMod1 <- mxAlgebra(expression = pathMod1c*defMod1, name = "cM1")
eMod1 <- mxAlgebra(expression = pathMod1e*defMod1, name = "eM1")
aMod2 <- mxAlgebra(expression = pathMod2a*defMod2, name = "aM2")
cMod2 <- mxAlgebra(expression = pathMod2c*defMod2, name = "cM2")
eMod2 <- mxAlgebra(expression = pathMod2e*defMod2, name = "eM2")
aMod3 <- mxAlgebra(expression = pathMod3a*defMod3, name = "aM3")
cMod3 <- mxAlgebra(expression = pathMod3c*defMod3, name = "cM3")
eMod3 <- mxAlgebra(expression = pathMod3e*defMod3, name = "eM3")
pathaMod <- mxAlgebra(expression = a +aM1+aM2+aM3, name = "aMod")
pathcMod <- mxAlgebra(expression = a +cM1+cM2+cM3, name = "cMod")
patheMod <- mxAlgebra(expression = a +eM1+eM2+eM3, name = "eMod")
} else if (length(modvars == 4)) {
aMod1 <- mxAlgebra(expression = pathMod1a*defMod1, name = "aM1")
cMod1 <- mxAlgebra(expression = pathMod1c*defMod1, name = "cM1")
eMod1 <- mxAlgebra(expression = pathMod1e*defMod1, name = "eM1")
aMod2 <- mxAlgebra(expression = pathMod2a*defMod2, name = "aM2")
cMod2 <- mxAlgebra(expression = pathMod2c*defMod2, name = "cM2")
eMod2 <- mxAlgebra(expression = pathMod2e*defMod2, name = "eM2")
aMod3 <- mxAlgebra(expression = pathMod3a*defMod3, name = "aM3")
cMod3 <- mxAlgebra(expression = pathMod3c*defMod3, name = "cM3")
eMod3 <- mxAlgebra(expression = pathMod3e*defMod3, name = "eM3")
aMod4 <- mxAlgebra(expression = pathMod4a*defMod4, name = "aM4")
cMod4 <- mxAlgebra(expression = pathMod4c*defMod4, name = "cM4")
eMod4 <- mxAlgebra(expression = pathMod4e*defMod4, name = "eM4")
pathaMod <- mxAlgebra(expression = a +aM1+aM2+aM3+aM4, name = "aMod")
pathcMod <- mxAlgebra(expression = a +cM1+cM2+cM3+cM4, name = "cMod")
patheMod <- mxAlgebra(expression = a +eM1+eM2+eM3+eM4, name = "eMod")
} else if (length(modvars == 5)) {
aMod1 <- mxAlgebra(expression = pathMod1a*defMod1, name = "aM1")
cMod1 <- mxAlgebra(expression = pathMod1c*defMod1, name = "cM1")
eMod1 <- mxAlgebra(expression = pathMod1e*defMod1, name = "eM1")
aMod2 <- mxAlgebra(expression = pathMod2a*defMod2, name = "aM2")
cMod2 <- mxAlgebra(expression = pathMod2c*defMod2, name = "cM2")
eMod2 <- mxAlgebra(expression = pathMod2e*defMod2, name = "eM2")
aMod3 <- mxAlgebra(expression = pathMod3a*defMod3, name = "aM3")
cMod3 <- mxAlgebra(expression = pathMod3c*defMod3, name = "cM3")
eMod3 <- mxAlgebra(expression = pathMod3e*defMod3, name = "eM3")
aMod4 <- mxAlgebra(expression = pathMod4a*defMod4, name = "aM4")
cMod4 <- mxAlgebra(expression = pathMod4c*defMod4, name = "cM4")
eMod4 <- mxAlgebra(expression = pathMod4e*defMod4, name = "eM4")
aMod5 <- mxAlgebra(expression = pathMod5a*defMod5, name = "aM5")
cMod5 <- mxAlgebra(expression = pathMod5c*defMod5, name = "cM5")
eMod5 <- mxAlgebra(expression = pathMod5e*defMod5, name = "eM5")
pathaMod <- mxAlgebra(expression = a +aM1+aM2+aM3+aM4+aM5, name = "aMod")
pathcMod <- mxAlgebra(expression = a +cM1+cM2+cM3+cM4+cM5, name = "cMod")
patheMod <- mxAlgebra(expression = a +eM1+eM2+eM3+eM4+eM5, name = "eMod")
}




for (i in names(defMod)) {
  for (j in names(pathModStore)) {
    if (grepl(i,j)) {
      namevec <- paste0("M",j)
      moderatedACE[[namevec]] <- mxAlgebra(expression = pathModStore[[j]]@name * defMod[[i]]@name, name = namevec)
    }
  }
}
View(moderatedACE)
##################################################################
##################################################################

# extract variables for bivariate moderation
getmodvars <- function(vector) {
varsnoarrow <- unlist(strsplit(vector, "->")) # remove arrow
onlyvars <- trimws(varsnoarrow) # remove white space (if there is some)
xvar <- onlyvars[1] # save X var
yvar <- onlyvars[2] # save Y var
result <- c(xvar,yvar)
}

if (!is.null(modACEbiv)) {
modACEbivvars <- lapply(modACEbiv, getmodvars)
}
if (!is.null(modBeta)) { # don't forget the additional condition: & type = "aceb" !!
modBetavars <- lapply(modBeta, getmodvars)
}

# the function (I still need to figure out how to implement the self referencing element here as in the loop)
#setfree <- function(condition, target) {
#condition <- as.vector(condition)
#target[condition[2],condition[1]] <- TRUE
#target
#}
# lapply(modACEbivvars, setfree, target = freemodpathB)

# modACEbiv

modmat <- matrix(0,nrow = 3,ncol = 3, dimnames = list(acevars,acevars))
modpathACEfree <- modmat!=0
for (i in modACEbivvars) {
modpathACEfree[as.vector(i)[2],as.vector(i)[1]] <- TRUE
modpathACEfree
}
modpathACEfree




pathAMod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = pathAModlabel,
                       free = modpathACEfree,
                       values = 0,
                       name = "aMod")
pathCMod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = pathAModlabel,
                       free = modpathACEfree,
                       values = 0,
                       name = "bMod")
pathEMod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = pathAModlabel,
                       free = modpathACEfree,
                       values = 0,
                       name = "bMod")

# define matrix of definition variables

# modBeta
modmat <- matrix(0,nrow = 3,ncol = 3, dimnames = list(acevars,acevars))
freemodpathBeta <- modmat!=0
for (i in modBetavars) {
freemodpathBeta[as.vector(i)[2],as.vector(i)[1]] <- TRUE
freemodpathBeta
}
freemodpathBeta



