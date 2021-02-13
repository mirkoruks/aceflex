# Implement moderation 
  library(OpenMx)
# moderation only of paths from ACE components on decomposed variables and from decomposed variables on decomposed variables

# new arguments:
  # moderation = FALSE ## TRUE if moderation needed
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

nvstring <- as.character(1:nv)
pathAModlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("ba",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathAModlabel[upper.tri(pathAlabel, diag = FALSE)] <- NA
pathAModlabel
nvstring <- as.character(1:nv)
pathCModlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("bc",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathCModlabel[upper.tri(pathAlabel, diag = FALSE)] <- NA
pathCModlabel
#if (type = "chol") {
nvstring <- as.character(1:nv)
pathAModlabel <- matrix(apply(expand.grid(nvstring, nvstring), 1, function(x) paste("be",x[2], x[1], sep="")), nrow = nv, ncol = nv, byrow = TRUE)
pathAModlabel[upper.tri(pathAlabel, diag = FALSE)] <- NA
pathAModlabel
#}

pathAMod <- mxMatrix(type = "Full", nrow = nv, ncol = nv, byrow = TRUE, 
                       labels = pathAModlabel,
                       free = modpathACEfree,
                       values = 0,
                       name = "bMod")
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



