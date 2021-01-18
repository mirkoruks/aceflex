# example script using twinflex

rm(list = ls())

library(dplyr)
library(OpenMx)
library(xtable)

data_orig <- read.csv(file = "C:/Users/Besitzer/Documents/Arbeit/Twinlife/Artikel/Netzwerke/Git/netzwerke/Update/data_wide.csv",
                      header = TRUE)
summary(data_orig)

source("C:/Users/Besitzer/Documents/Arbeit/Twinlife/twinflex/Git/twinflex.R")
twinflex(acevars = c("kultkapjahre"),covvars = c("age","posbez"),data = data_orig,sep = "_",zyg = "zyg")

