# test start values

rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)
library(umx)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)

# create binary variable
data_orig <- data_orig %>% mutate(schoolb_1 = ifelse(schoolhigh_1 <= 2, 1,
                                                     ifelse(schoolhigh_1 >= 3, 2, NA))) %>% 
    mutate(schoolb_2 = ifelse(schoolhigh_2 <= 2, 1,
                                                     ifelse(schoolhigh_2 >= 3, 2, NA))) %>% 
    mutate(kukabin_1 = ifelse(kultkapjahre_1 >= median(kultkapjahre_1, na.rm = TRUE),1,
                              ifelse(kultkapjahre_1 < median(kultkapjahre_1, na.rm = TRUE),2,NA))) %>% 
    mutate(kukabin_2 = ifelse(kultkapjahre_2 >= median(kultkapjahre_1, na.rm = TRUE),1,
                              ifelse(kultkapjahre_2 < median(kultkapjahre_1, na.rm = TRUE),2,NA)))

mz_data <- data_orig %>% filter(zyg ==1) 
dz_data <- data_orig %>% filter(zyg ==2)
ovars <- c("schoolhigh_1","schoolhigh_2","schoolb_1","schoolb_2","kukabin_1","kukabin_2")
mz_data[,ovars] <- umxFactor(mz_data[,ovars])
dz_data[,ovars] <- umxFactor(dz_data[,ovars])

# load twinflex function
source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")

    # twinflex
cov(data_orig[,c("posbez_1","negbez_1","posbez_2","negbez_2")], use = "pairwise")
var(data_orig[,c("posbez_1","negbez_1","posbez_2","negbez_2")], use = "pairwise")

source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")
datadef <- data_orig %>% filter(!(is.na(age) | is.na(iseiempmean_1) | is.na(iseiempmean_2)))

tf_uni_c_cc <- twinflex(acevars = c("posbez","negbez"),covvars = c("age","iseiempmean"),
                        data = datadef,sep = "_",
                        zyg = "zyg", covariance = FALSE) 
tf_uni_c_cc[["ModelSummary"]] 
