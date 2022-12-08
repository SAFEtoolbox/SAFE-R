#' Morris Sampling
#'
#' This function builds a matrix X of input samples to be used for the Elementary Effects Test, using a One-At-the-Time sampling strategy as described in Campolongo et al. (2011).
#'
#' @param  r  positive integer number, number of elementary effects. 
#' @param  xmin  vector \code{(M)}. Lower bounds of input ranges.            
#' @param  xmax  vector \code{(M)}. Upper bounds of input ranges.     
#' @param  L  positive, even number. Number of levels in the sampling grid   

#' @return \code{X} matrix of sampling datapoints where \code{EE} must be computed. This is a matrix with \code{r * (M + 1)} rows and \code{M} columns. Each row is a point in the input space. Rows are sorted in \code{r} blocks, each including \code{M + 1} rows. Within each block, points (rows) differ in one component at the time. Thus, each block can be used to compute one Elementary Effect \code{(EE_i)} for each model input \code{(i = 1, ...,M)}.     

#' @references Morris, M.D. (1991), Factorial sampling plans for preliminary computational experiments, Technometrics, 33(2).  

#' @export

Morris_sampling<- function(r,xmin,xmax,L){   

##############
# Check inputs
##############

stopifnot(is.scalar(r), r > 0, r == floor(r),
			is.numeric(xmin), is.numeric(xmax), 
			length(xmin) == length(xmax),
			xmax > xmin,
			is.scalar(L), L > 0, L == floor(L))
			
	if( L %% 2 != 0){ 
		L <- ceiling(L / 2) * 2
		warning(sprintf('input L must be even!\n Using L=%g instead of user-defined value\n', L))
		}
		
	M <- length(xmax)	
	Dr <- matrix(xmax - xmin, nrow = M + 1, ncol = M, byrow = TRUE)
	xminMat <- matrix(xmin, nrow = M + 1, ncol = M, byrow = TRUE)

###################
# Perform sampling 
###################

X <- matrix(nrow = r * (M + 1), ncol = M) # sampling points

k <- 1
for (i in 1:r){
	# sample datapoints
	Bstar <- Morris_orientation_matrix(M, L)
	# resort to original ranges:
	Bstar <-  xminMat + Dr * Bstar
	
	X[((i - 1) * (M + 1) + 1 ):(i  * (M + 1)),] <- t(sapply(1:(M+1), function(j) Bstar[j,]))
	
    } 

return(X)

}