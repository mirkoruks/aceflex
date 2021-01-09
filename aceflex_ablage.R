


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