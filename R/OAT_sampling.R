#' One-At-the-Time sampling strategy
#'
#' This function builds a matrix \code{X} of input samples to be used for the Elementary Effects Test, using a One-At-the-Time sampling strategy as described in Campolongo et al. (2011).
#'
#' @param r  positive integer number , number of elementary effects 
#' @param M  positive integer number , number of inputs
#' @param  distr_fun probability distribution function of each input. list (eg: \code{"unif"}) if all inputs have the same pdf or a list of length \code{M} strings (eg: \code{list("unif", "norm")}) otherwise. See help of \code{\link{AAT_sampling}} to check supported PDF types.
#' @param distr_par parameters of the probability distribution function - row vector if all input pdfs have the same parameters - list of \code{M} vectors otherwise
#' @param samp_strat sampling strategy - string Options: \code{"rsu"}: random uniform, \code{"lhs"}: latin hypercube.                            
#' @param des_type  design type - string. Options: \code{"trajectory"}, \code{"radial"}

#' @return \code{X} matrix of sampling datapoints where \code{EE} must be computed. This is a matrix with \code{r * (M + 1)} rows and \code{M} columns. Each row is a point in the input space. Rows are sorted in \code{r} blocks, each including \code{M + 1} rows. Within each block, points (rows) differ in one component at the time. Thus, each block can be used to compute one Elementary Effect \code{(EE_i)} for each model input \code{(i = 1, ...,M)}.  

#' @references Campolongo F., Saltelli, A. and J. Cariboni (2011), From screening to quantitative sensitivity analysis. A unified approach, Computer Physics Communications, 182(4), 978-988.

#' @seealso \code{\link{AAT_sampling}}

#' @export

#' @examples
#' # Example 1: 2 inputs, both from Unif[0,3]
#' r <-  10
#' M <-  2
#' distr_fun <-  "unif"
#' distr_par <-  c(0, 3)
#' samp_strat <-  "lhs"
#' des_type <-  "trajectory"
#' X <-  OAT_sampling(r, M, distr_fun, distr_par, samp_strat, des_type)
#' # Plot results:
#' plot(X[,1], X[,2], col = rep(rainbow(r), each = M + 1), pch = 19,
#' xlab = expression(x[1]), ylab = expression(x[2]))
#' for(k in 0:(r-1)){
#'	 segments(X[c(1,3) + (M+1) * k, 1], X[c(1,3) + (M+1) * k, 2],
#' X[2 + (M+1) * k, 1], X[2 + (M+1) * k, 2], lty = 2, col = "gray")
#' }
#' # Example 2: 2 inputs, one from Unif[0,3], one from Unif[1,5]
#' distr_fun <-  "unif" 
#' distr_par <- list(c(0, 3), c(1, 5))
#' X <-  OAT_sampling(r, M, distr_fun, distr_par, samp_strat, des_type)
#' plot(X[,1], X[,2], col = rep(rainbow(r), each = M + 1), pch = 19,
#' xlab = expression(x[1]), ylab = expression(x[2]))
#' for(k in 0:(r-1)){
#'	 segments(X[c(1,3) + (M+1) * k, 1], X[c(1,3) + (M+1) * k, 2],
#' X[2 + (M+1) * k, 1], X[2 + (M+1) * k, 2], lty = 2, col = "gray")
#' }


OAT_sampling <- function(r, M, distr_fun, distr_par, samp_strat, des_type){

##############
# Check inputs
##############

		stopifnot(is.scalar(r), r >0, r == floor(r),
			is.scalar(M), M >0, M == floor(M),
			length(distr_fun) %in% c(1, M),
			distr_fun %in% c('chisq', 'exp', 'geom', 'pois', 't',
				'beta', 'binom', 'f', 'gamma', 'lnorm', 'nbinom', 'norm', 'unif', 'weibull',
				'hyper', 'unifd'),
			samp_strat %in% c('rsu', 'lhs'),
			des_type %in% c('trajectory', 'radial')
		)
		
	if(!is.list(distr_par) & length(distr_fun) == 1){
		distr_par <- lapply(vector('list', M), function(h) distr_par)
	}
	
	stopifnot(is.list(distr_par), length(distr_par) == M)
	
##################
# Perform sampling
##################

	AB <- AAT_sampling(samp_strat, M, distr_fun, distr_par, 2 * r)

	X <- matrix(nrow = r * (M + 1), ncol = M) # sampling points
	k <- 1 
	
	if(des_type == 'radial'){
		for (i in 1:r) {
			# Sample datapoints:
			ab  <- AB[c(2 * i - 1, 2 * i), ] 
			a <- ab[1,] 
			b <- ab[2,]
			
			for (j in 1:M){
				
				if(a[j] == b[j]){
					if(distr_fun[[j]] == 'unifd'){
						while(a[j] == b[j]){
							tmp <- AAT_sampling(samp_strat, M, distr_fun, distr_par, 2)
							b[j] <- tmp[1, j]
						}
					warning(sprintf('b(i = %d, j = %d) was randomly changed not to overlap with a(i = %d, j = %d) \n', i, j, i, j))
					} else {
						warning(sprintf('b(i = %d, j = %d) and a(i = %d, j = %d) are the same! \n', i, j, a[j], b[j]))
					}
					
				}
				}
			
			X[k,] <- a 
			k <- k + 1
			
			for (j in 1:M){
				x2  <- a
				x2[j] <- b[j]
				X[k,] <- x2 
				if (abs(X[k, j] - X[k - 1, j]) == 0) {
					warning(sprintf('X(%d,%d) and X(%d,%d) are equal\n', k, j, k-1, j))
				}
			k <- k+1
		}
        } 
	} else { # (des_type == 'trajectory')
		
		for (i in 1:r) {
			# Sample datapoints:
			ab  <- AB[c(2 * i - 1, 2 * i), ] 
			a <- ab[1,] 
			b <- ab[2,]
			X[k,] <- a 
			k <- k + 1
			x  <- a  
     		   	for (j in 1:M) { 
     		   		x2 <- x 
            		x2[j] <- b[j]  
            		X[k,] <- x2 
            		if (abs(X[k, j] - X[k - 1, j]) == 0) {
            			warning(sprintf('X(%d,%d) and X(%d,%d) are equal\n', k, j, k-1, j))
            		}
            	# move one step forward in the trajectory
            	x <- x2 # move one step forward in the trajectory
            	k <- k+1 
            }
     		}
     	}
	
	return(X)
}