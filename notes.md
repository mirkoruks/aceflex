# to be added:
## the cogs and wheels of the moderation
- make sure to calculate the variance constraint for binary variables only with the unmoderated paths so that it also applies to the case with moderated paths and fixes the variance of the trait to 1 given a moderator value of 0
- if the moderator has within-pair variance and is not included in the covariance matrix as a ACE variable, make sure that twinflex calculates a extended univariate moderation model as suggested by van der Sluis et al. 2011

## output
- a list of the labels of *all* possible parameter, so in a Cholesky model include the beta paths (which are fixed to 0) as well, so the user can set them free while fixing any of the relevant cross-paths. Here It would be nice to include a RAM-ordered parameter list with a an A element (paths) and S matrix ((co-)variances) and a M element (means). These elements could be string matrices containing the labels of the parameters with the variable names as the labels of the dimensions. So the user can refer to the variance of X1 with a something like that: `parlist$S["X1","X1"]`
