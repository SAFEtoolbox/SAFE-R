#' Sobol' g-function
#'
#' Implements the Sobol' g-function, a standard benchmark function in
#' the Sensitivity Analysis literature (see for instance Sec. 3.6 in 
#' Saltelli et al. (2008)).
#'
#'
#' @param x vector (\code{M}) of inputs \eqn{x(1), x(2), ... x(M)}. \eqn{x(i)} ~ \eqn{Unif(0,1)} for all \eqn{i}
#' @param a vector (\code{M}) or scalar
#'
#' @return \code{y} scalar output. Two attributes \code{V} scalar output variance (exact value computed analytically), \code{Si_ex} vector (3), first-order sensitivity indices (exact value computed analytically)              
#'
#'
#' @references Saltelli et al. (2008) Global Sensitivity Analysis, The Primer, Wiley.

#' @export

#' @examples
#'
#' a <- 9 # options for the (fixed) parameters
#' M <- 5 # options for the number of inputs
#' y <- sobol_g_function(runif(M), a)
#' Si_ex <- attributes(y)$Si_ex


sobol_g_function <- function(x, a){
	
	M <- length(x) 
	
	if(length(a) == 1) a <- rep(a, M)
	
	g <- (abs( 4 * x - 2 ) + a) / (1 + a)
	y <- prod(g)
	
	# By model definition:
	# Vi = VAR(E(Y|Xi)) = 1 / (3 * (1 + ai)^2)
	# VARy = VAR(Y) = - 1 + prod( 1 + Vi )
	
	Vi <- 1 /(3 * (1 + a)^2) 
	V  <- - 1 + prod(1 + Vi)
	Si <- Vi / V
	
	attributes(y) <- list(V = V, Si_ex = Si)
	
	return(y)
	
}

