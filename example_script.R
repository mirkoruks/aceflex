# example script using twinflex

rm(list = ls())

library(OpenMx)
library(umx)
library(dplyr)
data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)

source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")
twinflex_result <- twinflex(acevars = "posbez",data = data_orig,sep = "_",zyg = "zyg", tryHard = FALSE)

table(data_orig$zyg)
dz=data_orig[data_orig$zyg==2,]
mz=data_orig[data_orig$zyg==1,]




# univariate ACE: identical
umx_result <- umxACE(selDVs = "posbez",sep = "_",mzData = mz, dzData = dz)
summary(umx_result) 
twinflex_result


# bivariate ACE: identical
twinflex_result <- twinflex(acevars = c("posbez","negbez"),data = data_orig,sep = "_",zyg = "zyg", tryHard = FALSE, type = "chol")
umx_result <- umxACE(selDVs = c("posbez","negbez"),sep = "_",mzData = mz, dzData = dz)
summary(umx_result) 
twinflex_result
