#' Empirical CDF
#' 
#' Compute the empirical CDF of the sample \code{x} and evaluate it at datapoints \code{xi}.

empiricalcdf <- function(x,xi) {

	stopifnot(is.numeric(x),is.numeric(xi))

N <- length(x)
x <- sort(x)
F <- seq(1:N)/N

iu <- which(!duplicated(x,fromLast=TRUE))
F <- F[iu]
N <- length(F)
x <- x[iu]

Fi <- rep(1,times=length(xi))
for(j in N:1)
	Fi[which(xi <= x[j])] <- F[j]
Fi[which(xi < x[1])] <- 0

return(Fi)
}