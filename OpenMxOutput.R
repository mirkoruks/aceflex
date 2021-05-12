# Funktion für OpenMx-Output (aufbauend auf Kohler et al. (2011))

# das einzige Argument der Funktion ist das Objekt mit dem gefitteten OpenMx/umx-Modell

openMxOutput <- function(modelfit) {
sig.level <- function(coef,std.error){
  ifelse(abs(coef/std.error)>= qnorm(.0005,lower.tail=FALSE),"$^{***}$",
         ifelse(abs(coef/std.error) >= qnorm(.005,lower.tail=FALSE), "$^{**}$",
                ifelse(abs(coef/std.error) >= qnorm(.025,lower.tail=FALSE), "$^{*}$","")
        )
 )
}

nobs <- c("N",as.character(dim(data)[1]))
loli <- c("-2LL",as.character(round(summary(modelfit)[["Minus2LogLikelihood"]], digits = 3)))
aic <- c("AIC",as.character(round(summary(modelfit)[["informationCriteria"]][1,3], digits = 3)))
bic <- c("BIC",as.character(round(summary(modelfit)[["informationCriteria"]][2,3], digits = 3)))

np <- 6
ACE.printcoefs <-
  function(modelfit){
    namelist <- rbind(names(modelfit@output[["estimate"]][1:np]),
                      rep("",length(names(modelfit@output[["estimate"]][1:np]))))
    coeflist <- rbind(paste(round(modelfit@output[["estimate"]][1:np],3),
                            sig.level(modelfit@output[["estimate"]][1:np],
                                      modelfit@output[["standardErrors"]][1:np]),sep=""),
                      paste("(",round(modelfit@output[["standardErrors"]][1:np],3),")",sep="")
   )
    list(output = rbind(cbind(as.vector(namelist),as.vector(coeflist)),as.vector(loli),as.vector(aic),as.vector(bic)))
  }
return(ACE.printcoefs(fitBestModel))
}

#
openMxOutput(fitBestModel)
