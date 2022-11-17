#' Sampling for the Fourier Amplitude Sensitivity Test (FAST)
#'
#' Implements sampling for the Fourier Amplitude Sensitivity Test (FAST; Cukier et al., 1978) and returns a matrix \code{X} of \code{N} input samples.
#'
#' @param  distr_fun string (eg: \code{'unif'}) if all inputs have the same pdf, or vector of \code{M} strings (eg: \code{c('unif','norm')}) containing the probability distribution function of each input.
#' @param distr_par row vector if all input pdfs have the same parameters, list of \code{M} vectors otherwise containing the parameters of the probability distribution function.
#' @param M scalar, number of inputs
#' @param Nharm scalar, interference factor, i.e.the number of higher harmonics to be considered (default is 4, taken from Saltelli et al. (1999; page 42))
#' @param omega vector (\code{M}), angular frequencies associated to inputs (default values computed by function \code{\link{generate_FAST_frequency}})
#' @param N odd scalar, number of samples (default is \code{2 * Nharm * max(omega) + 1} which is the minimum sampling size according to Cukier et al. (1978))
#'
#' @return \code{X} matrix of \code{N} input samples.

#' @references Cukier, R.I., Levine, H.B., and Shuler, K.E. (1978), Nonlinear Sensitivity Analyis of Multiparameter Model SYstems, Journal of Computational Physics, 16, 1-42.
#'
#' Saltelli, A., Tarantola, S. and Chan, K.P.S. (1999), A Quantitative Model-Independent Method ofr Global Sensitivty Analysis of Model Output, Technometrics, 41(1), 39-56.

#' @seealso \code{\link{FAST_sampling_unif}} \code{\link{generate_FAST_frequency}}

#' @examples 

#' fun_test  <- 'sobol_g_function'
#' # Define input distribution and ranges:
#' M <- 5 # may range from 5 to 11
#' distr_fun <- 'unif'
#' distrpar <-  c(0, 1)
#' a <- 9
#'
#' ## Step 1: Choose sampling size
#'
#' # suggested values of Nfast
#'
#' Nfast <- 2^13 + 1 
#' ## Step 2: Approximate first-order variance-based indices by FAST
#'
#' # FAST
#' Fsamp <- FAST_sampling(distr_fun, distrpar, M, N = Nfast)
#' X <- Fsamp$X
#' s <- Fsamp$s

#' @export


FAST_sampling <- function(distr_fun, distr_par, M, Nharm = 4, omega = generate_FAST_frequency(M), N = 2 * Nharm * max(omega) + 1){

##############
# Check inputs
##############

stopifnot(is.scalar(M), M >0, M == floor(M),
			length(distr_fun) %in% c(1, M),
			distr_fun %in% c('chisq', 'exp', 'geom', 'pois', 't',
				'beta', 'binom', 'f', 'gamma', 'lnorm', 'nbinom', 'norm', 'unif', 'weibull',
				'hyper'),
			is.scalar(N), N >0, N == floor(N), N %% 2 != 0,
	is.scalar(Nharm), Nharm >0, Nharm == floor(Nharm),
	is.numeric(omega), length(omega) == M, omega > 0, omega == floor(omega)
)
		
	if(!is.list(distr_par) & length(distr_fun) == 1){
		distr_par <- lapply(vector('list', M), function(h) distr_par)
	}
	
	stopifnot(is.list(distr_par), length(distr_par) == M)

###################################
# Recover and check optional inputs
###################################

if (N < 2 * Nharm * max(omega) + 1){ # and finally check that is is consistent with omega
	Nuser <- N
    N <- 2 * Nharm * max(omega) + 1
    cat(sprintf('Sample size specified by user (%d) is smaller than minimum sample size. Using the latter (%d) instead', Nuser, N))
    }

###################################
# Uniformly sample the unit square using FAST sampling
###################################

fsu <- FAST_sampling_unif(M, N = N, Nharm = Nharm, omega = omega)

X <- fsu$X
s <- fsu$s

###################################
# Map back into the specified distribution
# by inverting the CDF
###################################

if(length(distr_fun == 1)){
	distr_fun <- rep(distr_fun, M)
}


for(i in 1:M){
	
    name <- distr_fun[i]
    
    # create a list of elements to be passed in the inverse cdf function.
    # it contains X[, i] and the parameters of the distribution distr_par[[i]] 
    
    arglist <- list()
    arglist[[1]] <- X[, i]
    arglist <- c(arglist, distr_par[[i]])
    
    X[ ,i] <- do.call(sprintf("q%s", name), arglist) 
  }

robj <- list(X = X, s = s)

return(robj)

}