#' Main effect (first-order) sensitivity index for the Fourier Amplitude Sensitivity Test (FAST)
#'
#' Computes main effect (first-order) sensitivity index according to the Fourier Amplitude Sensitivity Test (FAST; Cukier et al., 1978).
#'
#' @param Y vector (\code{N}), set of model output samples 
#' @param M scalar, number of inputs
#' @param Nharm scalar, interference factor, i.e.the number of higher harmonics to be considered (default is 4)
#' @param omega vector (\code{M}), angular frequencies associated to inputs (default values computed by function \code{\link{generate_FAST_frequency}})
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{Si} (vector of length \code{M}) of main effect (first-order) sensitivity indices
#'   \item \code{V} (scalar) total output variance
#'   \item \code{A} (vector \code{N}) Fourier coefficients
#'   \item \code{B} (vector \code{N}) Fourier coefficients
#'   \item \code{Vi} (vector of length \code{M}) output variances from each input 
#' }
#'
#' @references Cukier, R.I., Levine, H.B., and Shuler, K.E. (1978), Nonlinear Sensitivity Analyis of Multiparameter Model SYstems, Journal of Computational Physics, 16, 1-42.
#'
#' Saltelli, A., Tarantola, S. and Chan, K.P.S. (1999), A Quantitative Model-Independent Method ofr Global Sensitivty Analysis of Model Output, Technometrics, 41(1), 39-56.
#' @seealso \code{\link{FAST_sampling}} \code{\link{generate_FAST_frequency}}

#' @examples 

#' fun_test  <- 'sobol_g_function'
#' # Define input distribution and ranges:
#' M <- 5 # may range from 5 to 11
#' distr_fun <- 'unif'
#' distrpar <-  c(0, 1)
#' a <- 9
#'
#' ## Step 1: Choose sampling size # 
#'
#' Nfast <- 2^5 + 1 # this value just for a quick example 
#' # it would be preferred: Nfast<- 2^13 + 1 
#' ## Step 2: Approximate first-order variance-based indices by FAST
#'
#' # FAST
#' Fsamp <- FAST_sampling(distr_fun, distrpar, M, N = Nfast)
#' X <- Fsamp$X
#' s <- Fsamp$s
#' # Run the model and compute selected model output at sampled parameter
#' # sets:
#' Y <- model_evaluation(fun_test, X, a = a)
#' # Estimate indices: 
#' Si_fast <- FAST_indices(Y, M)

#' @export

FAST_indices <- function(Y,M, Nharm = 4, omega = generate_FAST_frequency(M)){

##############
# Check inputs
##############


 stopifnot(is.scalar(M), M >=0, M == floor(M), 
 is.scalar(Nharm), Nharm >0, Nharm == floor(Nharm),
	is.numeric(omega), length(omega) == M, omega > 0, omega == floor(omega), ncol(Y) == 1)
 
 stopifnot(is.numeric(Y))

#########################################################
# Compute Fourier coefficients (vectors A and B) from Y
#########################################################

# The code below implements the equations given in Appendix C
# of Saltelli et al. (1999).

N <- length(Y)


Ymat <- matrix(Y[-1], 2, (N-1) / 2)
baseplus <- colSums(Ymat)
baseminus <- - apply(Ymat, 2, diff)

A <- numeric(N)
B <- numeric(N)

evenj <- seq(2, N, by = 2)
oddj <- seq(1, N, by = 2)

A[evenj] <- (Y[1] + colSums(baseplus * cos((1:((N-1) / 2) %o% evenj) * pi / N))) / N
B[oddj] <- colSums(baseminus * sin((1:((N - 1) / 2) %o% oddj) * pi / N)) / N

###################################
# Compute main effect from A and B
###################################

# The code below implements the equations given in Appendix B
# of Saltelli et al. (1999) (here we use 'V' and 'Vi' for the output
# variances while in that paper they are called 'D' and 'Di')

V <- 2 * sum(A^2 + B^2) # total output variance

idx <- omega %o% (1:Nharm)
Vi <- 2 * rowSums(matrix(A[idx]^2 + B[idx]^2, nrow = M))

Si <- Vi / V 

cat('\n \t main \n')
cat(sprintf(' X%d:\t %1.4f\t \n', 1:M, Si))
cat(sprintf('\n sum:\t %1.4f\t \n\n', sum(Si)))


robj <- list(Si = Si, V = V, A = A, B = B, Vi = Vi)

return(robj)

}    