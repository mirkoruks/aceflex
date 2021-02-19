# example script using twinflex

# clear workspace
rm(list = ls())

# load packages for twinflex
library(OpenMx)
library(umx)
library(dplyr)
library(psych)
# load example twin data "Australian twin sample biometric data"
data("twinData")

# load twinflex
source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")

# data preparation
twinData <- twinData %>% mutate(zyg_r = ifelse(zygosity %in% c("MZFF","MZMM"), 1, # 1 = MZ!
                                    ifelse(zygosity %in% c("DZFF","DZMM","DZOS"), 2, NA))) # 2 = DZ!

# age as a moderator for univariate GxE
twinData <- twinData %>% select(zyg_r, wt1, wt2, ht1, ht2, age1, age2, age)
# data summary 
summary(twinData)

# rescale height variable (height in cm and not m)
twinData <- twinData %>% 
    mutate(ht1 = ht1*100) %>% 
    mutate(ht2 = ht2*100) 
GxE <- twinData %>% 
    mutate(wt1 = (wt1-mean(wt1, na.rm = TRUE)) / sd(wt1, na.rm = TRUE)) %>% 
    mutate(wt2 = (wt2-mean(wt2, na.rm = TRUE)) / sd(wt2, na.rm = TRUE)) %>% 
    mutate(ht1 = (ht1-mean(ht1, na.rm = TRUE)) / sd(ht1, na.rm = TRUE)) %>% 
    mutate(ht2 = (ht2-mean(ht2, na.rm = TRUE)) / sd(ht2, na.rm = TRUE)) 

# prepare data for umx
    # uni- and bivariate
mz    <- twinData %>% filter(zyg_r == 1) %>% select(ht1,ht2,wt1,wt2,age1,age2)
dz    <- twinData %>% filter(zyg_r == 2) %>% select(ht1,ht2,wt1,wt2,age1,age2)

    # univariate GxE
mzGxE1 <- twinData %>% filter(zyg_r == 1 & !is.na(age1) & !is.na(age2)) %>% select(wt1,wt2,age1,age2)
dzGxE1 <- twinData %>% filter(zyg_r == 2 & !is.na(age1) & !is.na(age2)) %>% select(wt1,wt2,age1,age2)

    # bivariate GxE
mzGxE2 <- GxE %>% filter(zyg_r == 1 & !is.na(wt1) & !is.na(wt2)) %>% select(wt1,wt2,ht1,ht2)
dzGxE2 <- GxE %>% filter(zyg_r == 2 & !is.na(wt1) & !is.na(wt2)) %>% select(wt1,wt2,ht1,ht2)

# prepare data for twinflex (actually you only have to make sure that there are no NAs on the definition variables)
    # univariate GxE
twinDataGxE1 <- twinData %>% select(-c(age1,age2)) %>% filter(!is.na(age))
    # bivariate GxE
twinDataGxE2 <- GxE %>% filter(!is.na(wt1) & !is.na(wt2)) 

##############################
# Comparison of model results 
##############################

# univariate continuous (Weight)
    # twinflex
ace_uni_c_twinflex <- twinflex(acevars = "wt",data = twinData,sep = "",zyg = "zyg_r")
    # umx
ace_uni_c_umx <- umxACE(selDVs = "wt",sep = "",mzData = mz, dzData = dz)

    # compare twinflex and umx
ace_uni_c_twinflex[["ModelSummary"]] # Fit (-2lnL units) = 55506.49
summary(ace_uni_c_umx)               # Fit (-2lnL units) = 55506.49


# bivariate (Cholesky) continous (Weight - wt - and Height - ht)
    # twinflex
ace_biv_c_twinflex <- twinflex(acevars = c("wt","ht"),data = twinData,sep = "",zyg = "zyg_r", type = "chol")
    # umx
ace_biv_c_umx <- umxACE(selDVs = c("wt","ht"),sep = "",mzData = mz, dzData = dz)

# compare twinflex and umx
ace_biv_c_twinflex[["ModelSummary"]] # Fit (-2lnL units) = 103212.6
summary(ace_biv_c_umx)               # Fit (-2lnL units) = 103212.6

# univariate GxE 
    # twinflex
ace_gxe_univ_twinflex <- twinflex(acevars = "wt", data = twinDataGxE1, sep = "", zyg = "zyg_r", modACEuniv = "wt BY age")
    # umx
ace_gxe_univ_umx <- umxGxE(selDVs = "wt", selDefs = "age", sep = "", mzData = mzGxE1, dzData = dzGxE1,lboundACE = 1e-04)
ace_gxe_univ_umx <- umxModify(ace_gxe_univ_umx, update = "quad11") # umxGxE estimates a quadratic effect of the moderator on the means, by default; twinflex doesn't -> so dropping this parameter from the umx-models is necessary to compare the models
ace_gxe_univ_twinflex[["ModelSummary"]] # Fit (-2lnL units) = 55338.39
summary(ace_gxe_univ_umx)               # Fit (-2lnL units) = 55338.39

# bivariate GxE
    
    # twinflex
ace_gxe_biv_twinflex <- twinflex(acevars = c("wt","ht"), data = twinDataGxE2, sep = "", zyg = "zyg_r", modACEuniv = "ht BY wt", modACEbiv = "wt -> ht BY wt")
    # umx
twin_gxe_biv_umx <- umxGxEbiv(selDVs = "ht", selDefs = "wt", sep = "", mzData = mzGxE2, dzData = dzGxE2)

# compare twinflex and umx
ace_gxe_biv_twinflex[["ModelSummary"]] # Fit (-2lnL units) = 32435.14
#      |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
# AIC:        3529.14               32469.14                 32469.31
# BIC:      -86055.82               32574.51                 32520.49
summary(twin_gxe_biv_umx)              # Fit (-2lnL units) = 33300.54
#      |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
# AIC:       4394.539               33334.54                 33334.71
# BIC:     -85190.421               33439.91                 33385.89
###
# --> twinflex results in a slightly better fitting model than umx for the bivariate GxE
###





