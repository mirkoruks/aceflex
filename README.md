# aceflex
A flexible function for biometric twin modelling

Some remarks beforehand:
- As my markdown skills improve (or I find a proper cheat sheet), this readme becomes prettier and more readable.
- The explanation here remains short. Soon, you can find a more detailed documentation with examples [here](https://mirkoruks.github.io/).
- It's still a beta version. So if you find some bugs or have any comments, questions or suggestions, I am happy if you send me a mail to mirko.ruks@uni-bielefeld.de

At the moment you can use the function to estimate univariate or multivariate ACE-models. The multivariate ACE models follow the Cholesky parametrization. You can enter covariates to the model as well. At the moment, the covariates can be included as part of the covariance matrix of the model. There is also the option to include them into the mean vector as so called "definition variables". Currently, there is no option for a threshold model needed for ordinal/binary endogenous variables, but it will not take too much time to include it. 

So, this is the function with its arguments and the default values of the arguments (if existent):
twinflex(acevars, data, zyg, sep, covvars = NULL, optimizer = "SLSQP", tryHard = TRUE, tries = 10)

Some information about the arguments:
- acevars: a vector of strings, like e.g. c("extraversion", "education"), with the variables that you want to use for the ACE decomposition. It is important that you enter the variables without the twin specific suffix (see the argument "sep").
  - Two notes on the multivariate case (when there are more than one element in the acevars vector): 
      - The function follows the Cholesky parametrization, which implies directed paths from the ACE components of one phenotype on the other phenotype(s). The first element of the acevars vector is the first phenotype to be decomposed and would be the phenotype "on the left". The next element in the acevars vector is the next element in the causal chain, and so on. 
      - The default model is a ACE-beta-model (see Kohler et al. 2011), where there are no cross paths from the E components and insted there is a phenotypic effect (a "beta"). You can easily switch to the default Cholesky model by constraining the beta to zero and set the E-cross path free
- data: It's simply the data frame with the variables you want to decompose, etc.
- zyg: A string with the name of the variable with the zygosity information. Note, that I assume that 1 means MZ and 2 means DZ, so make sure that you recode the zygosity variable accordingly.
- sep: In the acevars argument you have to enter the variables without the twin-specific suffix. So instead of c("iqt1","iqt2), it's c("iq") and the separator ("t" in this case) can be entered as a string in the sep argument. If your variables are named like this: "iq_1" and "iq_2", the separator is "_"
- covvars: Here you can enter the covariates. Like for the acevars, you must not enter the variables with the twin specific suffix. You can enter variables with or without within-twin-pair variance. The function will check if there is within-pair variance automatically. The default is "covvars = NULL", so if you don't enter covariates, the function assumes that there are none. 

What's left are some OpenMx-specific arguments:
- optimizer: There are three options: "SLSQP", "NPSOL" and "CSOLNP". You can check the OpenMx documentation or the OpenMx forum to figure out the differences. "SLSQP" is the default.
- tryHard: If TRUE the function uses mxTryHard() to estimate the model, which means that it makes multiple attempts to fit the model with different starting values. It is the default, because sometimes bad starting values can ruin your model (although I tried to use wisely chosen starting values) 
- tries: Here you can specify how many attempts you want for mxTryHard


You can find the Kohler et al. paper here: https://doi.org/10.1080/19485565.2011.580619
