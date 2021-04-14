library(lavaan)

# example model
mzmodel <- '
# measurement model
LOC_1 =~ b1*loc0200_1 + b2*loc0202_1
EMO_1 =~ b3*cor1101_1 + b4*cor1102_1 + b5*cor1103_1
LOC_2 =~ b1*loc0200_2 + b2*loc0202_2
EMO_2 =~ b3*cor1101_2 + b4*cor1102_2 + b5*cor1103_2

# error variances
loc0200_1 ~~ e1*loc0200_1
loc0202_1 ~~ e2*loc0202_1
cor1101_1 ~~ e3*cor1101_1
cor1102_1 ~~ e4*cor1102_1
cor1103_1 ~~ e5*cor1103_1

loc0200_2 ~~ e1*loc0200_2
loc0202_2 ~~ e2*loc0202_2
cor1101_2 ~~ e3*cor1101_2
cor1102_2 ~~ e4*cor1102_2
cor1103_2 ~~ e5*cor1103_2

# Error variances latent traits
LOC_1~~0*LOC_1
EMO_1~~0*EMO_1
LOC_2~~0*LOC_2
EMO_2~~0*EMO_2

# phenotypic structural model
EMO_1 ~ breg*LOC_1
EMO_2 ~ breg*LOC_2

# ACE decomposition of latent factors

ALOC_1 =~ a11*LOC_1 + a21*EMO_1
AEMO_1 =~ a22*EMO_1
CLOC_1 =~ c11*LOC_1 + c21*EMO_1
CEMO_1 =~ c22*EMO_1 
ELOC_1 =~ e11*LOC_1 
EEMO_1 =~ e22*EMO_1  
ALOC_2 =~ a11*LOC_2 + a21*EMO_2
AEMO_2 =~ a22*EMO_2
CLOC_2 =~ c11*LOC_2 + c21*EMO_2
CEMO_2 =~ c22*EMO_2 
ELOC_2 =~ e11*LOC_2 
EEMO_2 =~ e22*EMO_2  

    # ACE Variances 
ALOC_1~~1*ALOC_1
CLOC_1~~1*CLOC_1
ELOC_1~~1*ELOC_1
ALOC_2~~1*ALOC_2
CLOC_2~~1*CLOC_2
ELOC_2~~1*ELOC_2
    # ACE Variances 
AEMO_1~~1*AEMO_1
CEMO_1~~1*CEMO_1
EEMO_1~~1*EEMO_1
AEMO_2~~1*AEMO_2
CEMO_2~~1*CEMO_2
EEMO_2~~1*EEMO_2

    # AC Covariances
AEMO_1~~1*AEMO_2
CEMO_1~~1*CEMO_2
ALOC_1~~1*ALOC_2
CLOC_1~~1*CLOC_2

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
