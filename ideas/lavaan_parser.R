library(lavaan)

# example model
mzmodel <- '

# phenotypic regressions twin 1
w3x1 ~ w2x1
w2x1 ~ w1x1
w3y1 ~ w2y1
w2y1 ~ w1y1
# phenotypic regressions twin 2
w3x2 ~ w2x2
w2x2 ~ w1x2
w3y2 ~ w2y2
w2y2 ~ w1y2

# ACE decomposition Within twin

# A
Aw1x1 =~ aw1xw1x*w1x1 + aw2xw1x*w2x1 + aw2yw1x*w2y1 
Aw1y1 =~ aw1yw1y*w1y1 + aw2xw1y*w2y1 + aw2xw1y*w2x1

Aw2x1 =~ aw2xw2x*w2x1 + aw3xw2x*w3x1 + aw3yw2x*w3y1 
Aw2y1 =~ aw2yw2y*w2y1 + aw3xw2y*w3y1 + aw3xw2y*w3x1

Aw3x1 =~ aw3xw3x*w3x1 
Aw3y1 =~ aw3yw3y*w3y1 

Aw1x2 =~ aw1xw1x*w1x2 + aw2xw1x*w2x2 + aw2yw1x*w2y2 
Aw1y2 =~ aw1yw1y*w1y2 + aw2xw1y*w2y2 + aw2xw1y*w2x2

Aw2x2 =~ aw2xw2x*w2x2 + aw3xw2x*w3x2 + aw3yw2x*w3y2 
Aw2y2 =~ aw2yw2y*w2y2 + aw3xw2y*w3y2 + aw3xw2y*w3x2

Aw3x2 =~ aw3xw3x*w3x2 
Aw3y2 =~ aw3yw3y*w3y2 

# C
Cw1x1 =~ cw1xw1x*w1x1 + cw2xw1x*w2x1 + cw2yw1x*w2y1 
Cw1y1 =~ cw1yw1y*w1y1 + cw2xw1y*w2y1 + cw2xw1y*w2x1

Cw2x1 =~ cw2xw2x*w2x1 + cw3xw2x*w3x1 + cw3yw2x*w3y1 
Cw2y1 =~ cw2yw2y*w2y1 + cw3xw2y*w3y1 + cw3xw2y*w3x1

Cw3x1 =~ cw3xw3x*w3x1 
Cw3y1 =~ cw3yw3y*w3y1 

Cw1x2 =~ cw1xw1x*w1x2 + cw2xw1x*w2x2 + cw2yw1x*w2y2 
Cw1y2 =~ cw1yw1y*w1y2 + cw2xw1y*w2y2 + cw2xw1y*w2x2

Cw2x2 =~ cw2xw2x*w2x2 + cw3xw2x*w3x2 + cw3yw2x*w3y2 
Cw2y2 =~ cw2yw2y*w2y2 + cw3xw2y*w3y2 + cw3xw2y*w3x2

Cw3x2 =~ cw3xw3x*w3x2 
Cw3y2 =~ cw3yw3y*w3y2 

# E
Ew1x1 =~ ew1xw1x*w1x1 
Ew1y1 =~ ew1yw1y*w1y1 

Ew2x1 =~ ew2xw2x*w2x1 
Ew2y1 =~ ew2yw2y*w2y1 

Ew3x1 =~ ew3xw3x*w3x1 
Ew3y1 =~ ew3yw3y*w3y1 

Ew1x2 =~ ew1xw1x*w1x2 
Ew1y2 =~ ew1yw1y*w1y2 

Ew2x2 =~ ew2xw2x*w2x2 
Ew2y2 =~ ew2yw2y*w2y2 

Ew3x2 =~ ew3xw3x*w3x2 
Ew3y2 =~ ew3yw3y*w3y2 

# ACE Covariances Within Twin Between Trait
# A
Aw1x1 ~~ raw1*Aw1y1
Aw2x1 ~~ raw2*Aw2y1
Aw3x1 ~~ raw3*Aw3y1

Aw1x2 ~~ raw1*Aw1y2
Aw2x2 ~~ raw2*Aw2y2
Aw3x2 ~~ raw3*Aw3y2

# C
Cw1x1 ~~ rcw1*Cw1y1
Cw2x1 ~~ rcw2*Cw2y1
Cw3x1 ~~ rcw3*Cw3y1

Cw1x2 ~~ rcw1*Cw1y2
Cw2x2 ~~ rcw2*Cw2y2
Cw3x2 ~~ rcw3*Cw3y2

# E
Ew1x1 ~~ rew1*Ew1y1
Ew2x1 ~~ rew2*Ew2y1
Ew3x1 ~~ rew3*Ew3y1

Ew1x2 ~~ rew1*Ew1y2
Ew2x2 ~~ rew2*Ew2y2
Ew3x2 ~~ rew3*Ew3y2

# ACE Covariances Between Twin
Aw1x1 ~~ 1*Aw1x2
Aw2x1 ~~ 1*Aw2x2
Aw3x1 ~~ 1*Aw3x2
Aw1x1 ~~ raw1*Aw1y2 
Aw2x1 ~~ raw1*Aw2y2
Aw3x1 ~~ raw1*Aw3y2

Aw1y1 ~~ 1*Aw1y2
Aw2y1 ~~ 1*Aw2y2
Aw3y1 ~~ 1*Aw3y2
Aw1y1 ~~ raw1*Aw1x2 
Aw2y1 ~~ raw1*Aw2x2
Aw3y1 ~~ raw1*Aw3x2

Cw1x1 ~~ 1*Cw1x2
Cw2x1 ~~ 1*Cw2x2
Cw3x1 ~~ 1*Cw3x2
Cw1x1 ~~ rcw1*Cw1y2 
Cw2x1 ~~ rcw1*Cw2y2
Cw3x1 ~~ rcw1*Cw3y2

Cw1y1 ~~ 1*Cw1y2
Cw2y1 ~~ 1*Cw2y2
Cw3y1 ~~ 1*Cw3y2
Cw1y1 ~~ raw1*Cw1x2 
Cw2y1 ~~ raw1*Cw2x2
Cw3y1 ~~ raw1*Cw3x2

'
# the dz model-object (only parameters that differ to those in the mz object)
dzmodel <-     "
# AC Covariances
AEMO_1~~0.5*AEMO_2
CEMO_1~~1*CEMO_2
ALOC_1~~0.5*ALOC_2
CLOC_1~~1*CLOC_2"

# parse input using the lavaan functions
mzparsed <- lavParseModelString(mzmodel)
mzpartable <- lavaanify(mzparsed, auto.var = TRUE, auto.cov.lv.x = FALSE)


dzparsed <- lavParseModelString(dzmodel)
dzpartable <- lavaanify(dzparsed)
dzpartable <- dzpartable[dzpartable$user == 1,]


ov <- lavaanNames(mzpartable, type = "ov")
lv <- lavaanNames(mzpartable, type = "lv")
dimvars <- c(ov,lv)


mzpartable$free <- ifelse(mzpartable$free != 0, 1,
                ifelse(mzpartable$free == 0, 0, NA))

# prepare parameter list for A matrix
Apar <- mzpartable[mzpartable$op %in% c("=~","~"),]
Apar$DV <- ifelse(Apar$op == "=~", Apar$rhs,
                  ifelse(Apar$op == "~", Apar$lhs,NA))
Apar$IV <- ifelse(Apar$op == "=~", Apar$lhs,
                  ifelse(Apar$op == "~", Apar$rhs,NA))

# prepare parameter list for S matrix
SparMZ <- mzpartable[mzpartable$op %in% c("~~"),]
SparMZ$rowvar <- SparMZ$lhs
SparMZ$colvar <- SparMZ$rhs
SparMZ$plabel <- paste0("cov",SparMZ$lhs,SparMZ$rhs)
SparDZ <- dzpartable[dzpartable$op %in% c("~~"),]
SparDZ$rowvar <- SparDZ$lhs
SparDZ$colvar <- SparDZ$rhs
SparDZ$plabel <- paste0("cov",SparDZ$lhs,SparDZ$rhs)

matrixtemplate <- matrix(NA, ncol = length(dimvars), nrow = length(dimvars), dimnames = list(dimvars, dimvars))

# construct A matrix
    # A labels matrix
label <- "NA" # information given by the user!
Alabel <- matrixtemplate

if (label == "all") {
labelmatrix <- paste(rep(paste0("b",dimvars), each = length(dimvars)), dimvars, sep = "")
Alabel <- matrix(labelmatrix, nrow = length(dimvars), ncol = length(dimvars), byrow = TRUE, dimnames = list(dimvars, dimvars)) # label matrix for matrix A. Explication: bDVIV with b = Effect symbol; DV = dep var; IV = indep var
} 
for (i in 1:dim(Apar)[1]) {
    if (Apar[i,"label"] != "") {
Alabel[Apar[i,"DV"],Apar[i,"IV"]] <- Apar[i,"label"]
    }
}
    # A free parameters
Afree <- matrixtemplate
Afree[is.na(Afree)] <- FALSE
Afree
for (i in 1:dim(Apar)[1]) {
    if (Apar[i,"free"] == 1) {
Afree[Apar[i,"DV"],Apar[i,"IV"]] <- TRUE
    } 
}
    # A starting values (at the moment default = 0.5)
Aval <- matrixtemplate
Aval[is.na(Aval)] <- 0
for (i in 1:dim(Apar)[1]) {
    if (Apar[i,"free"] == 1 & is.na(Apar[i,"ustart"])) { # if free and no starting value provided -> 0.5 (insert here an algorithmus to generate plausible starting values )
Aval[Apar[i,"DV"],Apar[i,"IV"]] <- 0.5 
}   else if (Apar[i,"free"] == 1 & !is.na(Apar[i,"ustart"])) { # if free and starting value provided -> take value provided
Aval[Apar[i,"DV"],Apar[i,"IV"]] <- Apar[i,"ustart"]
}   else if (Apar[i,"free"] == 0) { # if fixed -> take value provided
Aval[Apar[i,"DV"],Apar[i,"IV"]] <- Apar[i,"ustart"]
} 
}

# construct S matrix
    # S labels matrix
label <- "all" # information given by the user!
Slabel <- matrixtemplate

if (label == "all") {
labelmatrix <- paste(rep(paste0("cov",dimvars), each = length(dimvars)), dimvars, sep = "")
Slabel <- matrix(labelmatrix, nrow = length(dimvars), ncol = length(dimvars), byrow = TRUE, dimnames = list(dimvars, dimvars)) # label matrix for matrix A. Explication: bDVIV with b = Effect symbol; DV = dep var; IV = indep var
} 
for (i in 1:dim(SparMZ)[1]) {
    if (SparMZ[i,"label"] != "") {
Slabel[SparMZ[i,"rowvar"],SparMZ[i,"colvar"]] <- SparMZ[i,"label"]
Slabel[SparMZ[i,"colvar"],SparMZ[i,"rowvar"]] <- SparMZ[i,"label"]
    }
}
if (isSymmetric(Slabel) == TRUE) {
    print("We assume that the covariances don't differ within twin pairs! :-)")
} else {
    print("The maximum flexibility! We relax the assumption that the covariances don't differ within twin pairs! :-)")
}


    # S free parameters
Sfree <- matrixtemplate
Sfree[is.na(Sfree)] <- FALSE
for (i in 1:dim(SparMZ)[1]) {
    if (SparMZ[i,"free"] == 1) {
Sfree[SparMZ[i,"rowvar"],SparMZ[i,"colvar"]] <- TRUE
    } 
}

if (isSymmetric(Sfree) == TRUE) {
    print("Sfree is symmetric! :-)")
} else {
    stop("Ups, Sfree is not symmetric! :-(")
}

    # S starting values (at the moment default = 0.5)
        # mz
SvalMZ <- matrixtemplate
SvalMZ[is.na(SvalMZ)] <- 0
SvalDZ <- SvalMZ
for (i in 1:dim(SparMZ)[1]) {
    if (SparMZ[i,"free"] == 1 & is.na(SparMZ[i,"ustart"])) { # if free and no starting value provided -> 0.5 (insert here an algorithmus to generate plausible starting values )
SvalMZ[SparMZ[i,"rowvar"],SparMZ[i,"colvar"]] <- 0.5 
SvalMZ[SparMZ[i,"colvar"],SparMZ[i,"rowvar"]] <- 0.5 
}   else if (SparMZ[i,"free"] == 1 & !is.na(SparMZ[i,"ustart"])) { # if free and starting value provided -> take value provided
SvalMZ[SparMZ[i,"rowvar"],SparMZ[i,"colvar"]] <- SparMZ[i,"ustart"]
SvalMZ[SparMZ[i,"colvar"],SparMZ[i,"rowvar"]] <- SparMZ[i,"ustart"]
}   else if (SparMZ[i,"free"] == 0) { # if fixed -> take value provided
SvalMZ[SparMZ[i,"rowvar"],SparMZ[i,"colvar"]] <- SparMZ[i,"ustart"]
SvalMZ[SparMZ[i,"colvar"],SparMZ[i,"rowvar"]] <- SparMZ[i,"ustart"]
} 
}

for (i in 1:dim(SparDZ)[1]) {
    if (SparMZ[i,"free"] == 1 & is.na(SparMZ[i,"ustart"])) { # if free and no starting value provided -> 0.5 (insert here an algorithmus to generate plausible starting values )
SvalDZ[SparDZ[i,"rowvar"],SparDZ[i,"colvar"]] <- 0.5 
SvalDZ[SparDZ[i,"colvar"],SparDZ[i,"rowvar"]] <- 0.5 
}   else if (SparMZ[i,"free"] == 1 & !is.na(SparMZ[i,"ustart"])) { # if free and starting value provided -> take value provided
SvalDZ[SparDZ[i,"rowvar"],SparDZ[i,"colvar"]] <- SparDZ[i,"ustart"]
SvalDZ[SparDZ[i,"colvar"],SparDZ[i,"rowvar"]] <- SparDZ[i,"ustart"]
}   else if (SparMZ[i,"free"] == 0) { # if fixed -> take value provided
SvalDZ[SparDZ[i,"rowvar"],SparDZ[i,"colvar"]] <- SparDZ[i,"ustart"]
SvalDZ[SparDZ[i,"colvar"],SparDZ[i,"rowvar"]] <- SparDZ[i,"ustart"]
} 
}

if (isSymmetric(SvalMZ) == TRUE) {
    print("SvalMZ is symmetric! :-)")
} else {
    stop("Ups, SvalMZ is not symmetric! :-(")
}

if (isSymmetric(SvalDZ) == TRUE) {
    print("SvalDZ is symmetric! :-)")
} else {
    stop("Ups, SvalDZ is not symmetric! :-(")
}

# Matrix A and Matrix S as OpenMx-matrices 
matA <- mxMatrix(type = "Full", values = Aval, labels = Alabel, free = Afree, name = "A")
matSmz <- mxMatrix(type = "Full", values = SvalMZ, labels = Slabel, free = Sfree, name = "Smz")
matSdz <- mxMatrix(type = "Full", values = Svaldz, labels = Slabel, free = Sfree, name = "Sdz")

# start values
    # Paths: OLS Matrix (assuming for latent variables a variance of 1)
    # Variances: Empirical Variances for manifests, 1 for latents
    # Error variances: http://www.philender.com/courses/multivariate/notes/mr3.html



# missing:  
    # matrix M
    # threshold Matrix if there are ordinal endogenous variables
    # Moderation of A matrix
    # Covariates in the M matrix























### alte Sachen ab hier!
########################
# 1 delete white space
model1 <- gsub(" ", "", model)
model1
# 2 delete commented lines
model2 <- gsub("(\n#.*?\n)|(#.*?\n)","\n",model1) # funktioniert nicht bei zwei aufeinanderfolgendenen Zeilen mit Kommentaren
model2
if (grepl("#",model2)) {
    stop("Es gibt immer noch Kommentare")
    print("Alles ok")
}

model3 <- gsub("\n+", " ", model2)
model3 <- gsub("^ | $", "", model3)
model4 <- unlist(strsplit(model3, split = " "))
model4
latents <- grep("=~", model4, value = TRUE)
latents
# [:alnum:]


# A matrix: this string contains all the information needed to construct the A matrix 
# S matrix: this string contains all the within-twin (co-)variances, the cross-twin (co-)variances need to be added according to model type (cholesky or variance component)
# F matrix: this string contains all the information needed to construct the F matrix
# M matrix: this string contains all the information needed to construct the M matrix
# what I would need to add: information whether cholesky or variance component approach should be used 

# how to add definition variables
    # in the mean as covariates (e.g. age and sex as a covariate)
defmean <- '
# Mean structure
    # free means for manifest variables
loc0200 + loc0202 ~ mean*1
cor1101 + cor1102 + cor1103 ~ mean*1
mean = age + sex
'
# read "mean = age + sex" as: 
    # the parameter "mean" is a function of "age" and "sex"
        # that includes the following parameters:
            # meani = intercept (-> the value of mean when age=sex=0)
            # bmeanage = effect of age 
            # bmeansex = effect of sex
    # in the covariance matrix to moderate paths
acemod <- '
# ACE decomposition
A =~ a*var
C =~ c*var
E =~ e*var

# moderation of ACE paths
a = age + sex
c = age + sex
e = age + sex

# ACE variances
A~~1*A
C~~1*C
E~~1*E

# Mean structure
    # free means for manifest variables
loc0200 + loc0202 ~ mean*1
cor1101 + cor1102 + cor1103 ~ mean*1
mean = mu + bage*age + bsex*sex
'
# read "a = age" as: 
    # the parameter "a" is a function of "age" and "sex"
        # that includes the following parameters:
            # ai = intercept (-> the value of a when age=sex=0)
            # baage = effect of age 
            # basex = effect of sex
