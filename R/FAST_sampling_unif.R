#' Sampling for the Fourier Amplitude Sensitivity Test (FAST)
#'
#' Implements sampling for the Fourier Amplitude Sensitivity Test (FAST; Cukier et al., 1978) and returns a matrix \code{X} of \code{N} input samples.
#' Inputs are assumed to be uniformly distributed in the unit hypercube [0,1]^\code{M}. 
#' Samples are taken along the search curve defined by transformations
#'      \deqn{x_i(s) = G_i( \sin( \omega_i*s ) ) \; \;   i=1,...,M}
#' where \eqn{s} is a scalar variable that varies in \eqn{(-\frac{1}{\pi}, \frac{1}{\pi})}
#' Notes: Here we use the curve proposed by Saltelli et al. (1999):
#'         \deqn{x_i(s) = \frac{1}{2} + \frac{1}{\pi}* \arcsin(\sin(\omega_i*s ) )}
#'
#' @param M scalar, number of inputs
#' @param N odd scalar, number of samples (default is \code{2 * Nharm * max(omega) + 1} which is the minimum sampling size according to Cukier et al. (1978))
#' @param Nharm scalar, interference factor, i.e.the number of higher harmonics to be considered (default is 4, taken from Saltelli et al. (1999; page 42))
#' @param omega vector (\code{M}), angular frequencies associated to inputs (default values computed by function \code{\link{generate_FAST_frequency}})
#'
#' @seealso \code{\link{FAST_sampling}} \code{\link{generate_FAST_frequency}}
#'
#' @references Cukier, R.I., Levine, H.B., and Shuler, K.E. (1978), Nonlinear
#' Sensitivity Analyis of Multiparameter Model SYstems, Journal of
#' Computational Physics, 16, 1-42.
#'
#' Saltelli, A., Tarantola, S. and Chan, K.P.S. (1999), A Quantitative
#' Model-Independent Method ofr Global Sensitivty Analysis of Model Output,
#' Technometrics, 41(1), 39-56.

#' @export

FAST_sampling_unif <- function(M, Nharm = 4, omega = generate_FAST_frequency(M), N = 2 * Nharm * max(omega) + 1){

##############
# Check inputs
##############

stopifnot(is.scalar(M), M >0, M == floor(M),is.scalar(N), N >0, N == floor(N), N %% 2 != 0,
	is.scalar(Nharm), Nharm >0, Nharm == floor(Nharm),
	is.numeric(omega), length(omega) == M, omega > 0, omega == floor(omega)
)

if (N < 2 * Nharm * max(omega) + 1){ # and finally check that is is consistent with omega
	Nuser <- N
    N <- 2 * Nharm * max(omega) + 1
    warning('Sample size specified by user (%d) is smaller than minimum sample size. Using the latter (%d) instead',Nuser,N)
    }

###################################
# Perform sampling over the search curve
###################################

s <- pi / 2 * (2 * seq(1:N) - N - 1) / N 

###################################
# Map back sampled points in the input space
###################################

X <-  1 / 2 + 1 / pi * asin(sin(s %o% omega))

robj <- list(X = X, s = s)

return(robj)

}