# example script using twinflex

rm(list = ls())

library(OpenMx)
library(umx)
library(dplyr)
data("twinData")

summary(twinData)

data <- twinData %>% mutate(zyg_r = )
source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")
twinflex_result <- twinflex(acevars = "posbez",data = twinData,sep = "_",zyg = "zyg", tryHard = FALSE)

table(twinData$zyg)
dz=twinData[twinData$zyg==2,]
mz=twinData[twinData$zyg==1,]




# univariate ACE: identical
umx_result <- umxACE(selDVs = "posbez",sep = "_",mzData = mz, dzData = dz)
summary(umx_result) 
twinflex_result


# bivariate ACE: identical
twinflex_result <- twinflex(acevars = c("posbez","negbez"),data = data,sep = "_",zyg = "zyg", tryHard = FALSE, type = "chol")
umx_result <- umxACE(selDVs = c("posbez","negbez"),sep = "_",mzData = mz, dzData = dz)
summary(umx_result) 
twinflex_result
