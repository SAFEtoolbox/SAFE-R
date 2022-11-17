#' Resampling strategy needed to build the approximators of the first-order and total order sensitivity indices
#'
#' This function implements the resampling strategy needed to build the
#' approximators of the first-order (main effects) and total order
#' sensitivity indices (e.g. Saltelli et al. 2008; 2010).
#' This function is meant to be used in combination with \code{\link{vbsa_indices}}. 
#'
#' @param X matrix \code{(NX, M)} of \code{NX} input samples
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{XA} matrix \code{(N, M)} first \code{N = NX / 2} rows of \code{X}
#'   \item \code{XB} matrix \code{(N, M)} last \code{N = NX / 2} rows of \code{X}
#'   \item \code{XC} matrix \code{(N*M, M)} Block matrix of \code{M} recombinations of \code{XA} and \code{XB}, each block \code{XCi} is a \code{(N, M)} matrix whose columns are all taken from \code{XB} exception made for \code{i}-th, which is taken from \code{XA}.
#'}
#'
#' @references Saltelli et al. (2008), Global Sensitivity Analysis, The Primer, Wiley.
#'
#' Saltelli et al. (2010), Variance based sensitivity analysis of model 
#' output. Design and estimator for the total sensitivity index, Computer 
#' Physics Communications, 181, 259-270.

#' @seealso \code{\link{vbsa_indices}}

#' @export

#' @examples
#'
#' fun_test  <- "ishigami_homma_function"
#' M <- 3 
#' distrfun <- "unif"
#' distrpar  <- c(-pi, pi)
#' N <- 1000
#' sampstrat <- "lhs"
#' X <- AAT_sampling(sampstrat, M, distrfun, distrpar, 2 * N)
#' XABC <- vbsa_resampling(X)

vbsa_resampling <- function(X){

NX <- nrow(X)
M <- ncol(X)

if (NX %% 2 != 0){
	NX <- NX - 1
	cat(sprintf('\n WARNING: input matrix X has an odd number of rows, using the first %d rows only.\n', NX))
	}

N  <- NX / 2    
XA <- X[1:N,]
XB <- X[-(1:N),]
XC <- matrix(t(XB), nrow = N * M, ncol= M, byrow = TRUE)


for(i in 1:M){
    XC[((i-1)*N+1):((i-1)*N+N), i] <- XA[,i]
}

return(list(XA = XA, XB = XB, XC = XC))

}