#' Create an expanded sample \code{X_new} starting from a sample \code{X}
#'
#' This function create an expanded sample \code{X_new} starting from a sample \code{X} and using latin hypercube and the maximin criterion.
#'
#' @param X matrix (\code{N},\code{M}) of initial samples.
#' @param distr_fun string (eg: \code{"unif"}) if all inputs have the same pdf, or vector of \code{M} strings (eg: \code{c("unif", "norm")}) containing the probability distribution function of each input.
#' @param distr_par row vector if all input pdfs have the same parameters, list of \code{M} vectors otherwise containing the parameters of the probability distribution function.
#' @param N_new positive integer, new dimension of the sample (must be > \code{N})
#' @param nrep scalar, number of replicate to select the maximin hypercube (default value: 10)

#' @return \code{X_new} matrix (\code{N_new}, \code{M}) of expanded sample.  

#' @seealso \code{\link{AAT_sampling}}

#' @export
#' @examples
#' # Example 1: 2 inputs, both from Unif[0,3]
#' N <- 1000 
#' M <- 2 
#' distr_fun <- "unif" 
#' distr_par <- c(0, 3)
#' samp_strat <- "lhs"
#' X <- AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
#' # Plot results:
#' plot(X, xlab = expression(x[1]), ylab = expression(x[2]))
#' #
#' # Adding up new samples
#' N2 <- 500 # increase of base sample size
#' # (that means: N2 * (M+2) new samples that will need to be evaluated)
#' Xext <- AAT_sampling_extend(X, distr_fun, distr_par, 2 * (N + N2)) # extended sample 
#' # (it includes the already evaluated samples X and the new ones)
#' Xnew <- Xext[-(1:(2 * N)),] # extract the new input samples that need to be evaluated
#' # Plot results:
#' par(mfrow = c(1, 2))
#' plot(X, xlab = expression(x[1]), ylab = expression(x[2]), main = "X")
#' plot(Xnew, xlab = expression(x[1]), ylab = expression(x[2]), main = "X new")

AAT_sampling_extend <- function(X, distr_fun, distr_par, N_new, nrep = 10){

##############
# Check inputs
##############

stopifnot(is.matrix(X))

N <- nrow(X)
M <- ncol(X)

		stopifnot(length(distr_fun) %in% c(1, M),
			distr_fun %in% c('chisq', 'exp', 'geom', 'pois', 't',
				'beta', 'binom', 'f', 'gamma', 'lnorm', 'nbinom', 'norm', 'unif', 'weibull',
				'hyper', 'unifd'),
			is.scalar(N_new), N_new > N, N_new == floor(N_new),
			is.scalar(nrep), nrep > 0, nrep == floor(nrep)
		)
		
	if(!is.list(distr_par) & length(distr_fun) == 1){
		distr_par <- lapply(vector('list', M), function(h) distr_par)
	}
	
	stopifnot(is.list(distr_par), length(distr_par) == M)
	
	    if(length(distr_fun == 1)){
    	distr_fun <- rep(distr_fun, M)
    }

############################
# Map back the original sample into
# the uniform unit square
############################
	
	U <- sapply(1:M, function(i) do.call(sprintf("p%s", distr_fun[i]), c(list(X[, i]), as.list(distr_par[[i]])) ))
	
############################
# Add samples in the unit square
############################

U_new <- extend.lhcube(U = U, n = N_new, nrep = nrep)

row.names(U_new) <- paste('sample', 1:N_new, sep ="")

############################
# Map back into the specified distribution
# by inverting the CDF
############################
    	
    X_new <- sapply(1:M, function(i) do.call(sprintf("q%s", distr_fun[i]), c(list(U_new[, i]), as.list(distr_par[[i]])) ))
	
	return(X_new)
}