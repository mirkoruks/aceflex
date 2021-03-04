# brain storming for dynamic start values of path coefficients (-> less problems when there are mayor differences in variances because then a starting value of e.g. 0.3 is not appropriate!)


library(OpenMx)
library(dplyr)
#library(haven)
library(foreign)

options(scipen=999)
options(digits=2)
mxOption(NULL,"Default optimizer","SLSQP")

# Daten laden
data <- read.csv(file = "C:/Users/Besitzer/Desktop/Neuer Ordner/blabla.csv",
                      header = TRUE)

reg <-lm(pop0100 ~ sex+ per0100+ per0200, data = data)
summary(reg)

data_1 <- data[complete.cases(data),]
Y <- data_1$pop0100
Y <- as.matrix(Y)
X <- data_1[,c("sex","per0100","per0200")]
X <- as.matrix(X)
dim(X)
X <- cbind(rep(1,dim(X)[1]),X)

# Formula to calculate betas for phenotypic beta effects between acevars and covariate effects
solve(t(X)%*%X)%*%t(X)%*%Y

# Calculate start values of effects of latent ACE factors on acevars
    # assumption: correlation = 0.3 for univariate effects for all ACE factors 
        # new argument: coruniv = c(0.3,0.3,0.3)
            # explication: assumed correlation of unique A factor with manifest = 0.3, of unique C factor = 0.3 of unique E factor = 0.3, the order is A, C, E
    # assumption: correlation = 0.2 for bivariate effects for all ACE factors
        # new argument: corbiv = c(0.2,0.2,0.2) 
            # explication: same as above

# formula of bivariate beta (regression): beta = cov(x,y)/var(x)    
    # we know that cor(x,y)=cov(x,y)/sd(x)*sd(y) and we assume that cor(x,y) = 0.3/0.2 while we constrain the variance of the latent factors always to 1 -> var(x)=sd(x)=1
    # so we get 0.3=cov(x,y)/sd(y) <=> cov(x,y)=0.3*sd(y) for the univariate effects
    # and 0.2=cov(x,y)/sd(y) <=> cov(x,y)=0.2*sd(y) for the bivariate effects
    # substituting the terms in the formula for the beta, leads us to:
        # beta = cov(x,y)/var(x) <=> beta = 0.3*sd(y) for the univariate effects
        # beta = cov(x,y)/var(x) <=> beta = 0.2*sd(y) for the bivariate effects

# additional argument: coracerandom = FALSE -> if TRUE the values provided in coruniv and corbiv are used as start values for a random draw of assumed correlations
    # coracerandomtrials = 10
    # coracerandomdist = "normal" # alternative values = "uniform" or "cautchy" or "log" -> values provided in coruniv and corbiv are used as mean of distribution as well, sd of distribution is set to 1
