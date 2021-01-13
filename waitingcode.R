# fertiger code -> einfügen, wenn es passt

# path of constant covariates 
## HERE IF CONDITION: IF COVVARS == CONSTANT
pathCovlabelFunction <- function(string) {
paste0("b",rep(paste0(string,1:nv),2),rep(c(1,2),each=nv))
}
pathCovlabel <- as.vector(sapply(covvars,pathCovlabelFunction))
pathCovlabel

pathCovconstant <- mxMatrix(type = "Full", nrow = ntv, ncol = c, byrow = FALSE,
                            free = TRUE,
                            values = .3,
                            labels = pathCovlabel,
                            name = "pCov")