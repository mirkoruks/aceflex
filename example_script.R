# example script using twinflex

# clear workspace
rm(list = ls())

# load packages for twinflex
library(OpenMx)
library(umx)
library(dplyr)

# load example twin data "Australian twin sample biometric data"
data("twinData")
# summary of example twin data
summary(twinData)

# load twinflex
source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")

# data preparation
twinData <- twinData %>% mutate(zyg_r = ifelse(zygosity %in% c("MZFF","MZMM"), 1, # 1 = MZ!
                                    ifelse(zygosity %in% c("DZFF","DZMM","DZOS"), 2, NA))) # 2 = DZ!

# prepare data frames for umx
mz    <- subset(twinData, zyg_r==1,c("ht1","ht2","wt1","wt2"))
dz    <- subset(twinData, zyg_r==2,c("ht1","ht2","wt1","wt2"))



##############################
# Comparison of model results 
##############################

# univariate continuous (Weight)
    # twinflex
ace_uni_c_twinflex <- twinflex(acevars = "wt",data = twinData,sep = "",zyg = "zyg_r")
    # umx
ace_uni_c_umx <- umxACE(selDVs = "wt",sep = "",mzData = mz, dzData = dz)
    # results
summary(ace_uni_c_twinflex)
summary(ace_uni_c_umx) 


# bivariate (Cholesky) continous (Weight - wt - and Height - ht)
    # twinflex
ace_biv_c_twinflex <- twinflex(acevars = c("wt","ht"),data = twinData,sep = "",zyg = "zyg_r", type = "chol")
    # umx
ace_biv_c_umx <- umxACE(selDVs = c("wt","ht"),sep = "",mzData = mz, dzData = dz)

summary(ace_biv_c_twinflex)
summary(ace_biv_c_umx) 








