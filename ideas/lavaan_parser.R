# models
rm(list = ls())
library(OpenMx)
setwd("C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Corona/Daten")
df_wide <- read.csv("data_wide.csv", header = TRUE)
df_long <- read.csv("data_long.csv", header = TRUE)

# model:
    # int. locus of control -> emotional stress
    # int. locus of control measured by: loc0200, loc0202	
    # emotional stress measured by: cor1101, cor1102, cor1103						

# what would be the lavaan code for the model for one twin
# ACE-beta model with latent traits
model <- '
# measurement model
LOC =~ NA*loc0200 + loc0202
EMO =~ NA*cor1101 + cor1102 + cor1103

# Variances latent traits
LOC~~1*LOC
EMO~~1*EMO

# phenotypic structural model
EMO ~ LOC

# ACE decomposition of latent factors
    # Locus of Control
ALOC =~ a11*LOC + a21*EMO
CLOC =~ c11*LOC + c21*EMO
ELOC =~ e11*LOC 

    # ACE Variances 
ALOC~~1*ALOC
CLOC~~1*CLOC
ELOC~~1*ELOC

    # Emotional stress
AEMO =~ a11*EMO
CEMO =~ c11*EMO 
EEMO =~ e11*EMO  

    # ACE Variances 
AEMO~~1*AEMO
CEMO~~1*CEMO
EEMO~~1*EEMO

# Mean structure
    # free means for manifest variables
loc0200 + loc0202 ~ 1
cor1101 + cor1102 + cor1103 ~ 1

    # fixed means (to 0) for latent variables
LOC + EMO + ALOC + CLOC + ELOC + AEMO + CEMO + EEMO ~ 0*1
'

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
