#################################################################################################
#################################################################################################
#################################################################################################
## CHECK IF VARIABLES EXIST AND/OR ARE CONSTANT IF WITHIN VARIANCE
#################################################################################################
#################################################################################################
#################################################################################################

sep <- "_t"
# Build elements to construct expected covariance matrix (RAM Notation)
acevars <- c("posbez")

acevars1 <-    paste0(acevars,sep,"1") # ACE vars twin 1
acevars2 <-    paste0(acevars,sep,"2") # ACE vars twin 2
acevarswide <- c(acevars1, acevars2) # ACE vars variable vector (input for SEM)

covvars <- c("age","negbez")
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
existence_check_acevars


existence_check_acevars2 <- existence_check_acevars[grepl(c(suf1), existence_check_acevars)]


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
covvars <- c("age",)
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

#################################################################################################
#################################################################################################
#################################################################################################







# Check if acevars and covvars have within-twin-pair-variance -> give warning if correlation > .9
acevarscor <- mapply(cor,data[,varswide1],data[,varswide2], use = "pairwise.complete.obs")
#acevarscor[3] <- .98
s = attr(acevarscor, "names")
s1 = unlist(sapply(strsplit(s, split=sep, fixed=TRUE), function(x) (x[1])))
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


pathB <- mxMatrix(type = "Lower", nrow = nv, ncol = nv, byrow = TRUE,
                  free = freepathB,
                  values = valuespathB,
                  labels = pathBlabel,
                  name = "b")

pathZ <- mxMatrix(type = "Zero", nrow = nv, ncol = nv, name = "pZ")


cor(data_orig[,acevarswide], use = "pairwise.complete.obs")
acevarscor <- mapply(cor,data_orig[,acevars1],data_orig[,acevars2], use = "pairwise.complete.obs")
acevarscor[3] <- .98
s = attr(acevarscor, "names")
s1 = unlist(sapply(strsplit(s, split='_', fixed=TRUE), function(x) (x[1])))
attr(acevarscor, "names") <- s1
str(acevarscor)
lowcorrelation <- acevarscor < 0.9
alllowcorrelation <- all(lowcorrelation)
if (alllowcorrelation == FALSE) {
    print("Oups! For the following variables the within-pair correlation is >= 0.9. There might be estimation problems due to (multi-)collinearity")
    highcorrelation <- acevarscor[acevarscor >=.9]
    highcorrelationnames <- attributes(highcorrelation)
    print(highcorrelation)
    proceedcollinearity <- readline("Do you want to proceed? If so, type in 'yes', if not, type in 'no'") 
    if (proceedcollinearity == "no") {
        break
    }
}



loop_length <- seq_along(acevarscor)
for (i in loop_length) {
    if (acevarscor[i]  > 0.9) {
   print(paste0("Ups! The within-pair correlation for ", acevars[i], " is > 0.9"))
} else {
       print(paste0("Ok! The within-pair correlation for ", acevars[i], " is not > 0.9"))
}
=======
# ablage 



cor(data_orig[,acevarswide], use = "pairwise.complete.obs")
acevarscor <- mapply(cor,data_orig[,acevars1],data_orig[,acevars2], use = "pairwise.complete.obs")
acevarscor[3] <- .98
s = attr(acevarscor, "names")
s1 = unlist(sapply(strsplit(s, split='_', fixed=TRUE), function(x) (x[1])))
attr(acevarscor, "names") <- s1
str(acevarscor)
lowcorrelation <- acevarscor < 0.9
alllowcorrelation <- all(lowcorrelation)
if (alllowcorrelation == FALSE) {
    print("Oups! For the following variables the within-pair correlation is >= 0.9. There might be estimation problems due to (multi-)collinearity")
    highcorrelation <- acevarscor[acevarscor >=.9]
    highcorrelationnames <- attributes(highcorrelation)
    print(highcorrelation)
    proceedcollinearity <- readline("Do you want to proceed? If so, type in 'yes', if not, type in 'no'") 
    if (proceedcollinearity == "no") {
        break
    }
}



loop_length <- seq_along(acevarscor)
for (i in loop_length) {
    if (acevarscor[i]  > 0.9) {
   print(paste0("Ups! The within-pair correlation for ", acevars[i], " is > 0.9"))
} else {
       print(paste0("Ok! The within-pair correlation for ", acevars[i], " is not > 0.9"))
}
>>>>>>> 6b2c1e6ff10920cb8194250c757ae35b4510bad8
}