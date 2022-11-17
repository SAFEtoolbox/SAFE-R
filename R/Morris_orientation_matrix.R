#' Morris orientation matrix
#'
#' This function builds a Morris orientation matrix
#'
#' @param k integer positive scalar, number of inputs
#' @param p integer positive even scalar, number of levels 

#' @return \code{Bstar} matrix of \code{(k + 1)} datapoints in the \code{k}-dimensional input space (to be used for computing one Elementary Effect for each of the \code{k} inputs)   

#' @export
#' @examples
#' # Example in two-dimensional space (k = 2):
#' #
#' p <- 4
#' plot(1, xlim = c(0,1), ylim = c(0, 1), type ="n", xlab = "", ylab = "")
#' Bstar <- Morris_orientation_matrix(2,p)
#' lines(Bstar[,1], Bstar[,2], col = "red")
#' points(Bstar[,1],Bstar[,2], col = "red", pch = 20)
#' # if you want to generate more datapoints:
#' Bstar <- Morris_orientation_matrix(2,p)
#' lines(Bstar[,1], Bstar[,2], col = "red")
#' points(Bstar[,1],Bstar[,2], col = "red", pch = 20)
#' Bstar <- Morris_orientation_matrix(2,p)
#' lines(Bstar[,1], Bstar[,2], col = "red")
#' points(Bstar[,1],Bstar[,2], col = "red", pch = 20)
#' Bstar <- Morris_orientation_matrix(2,p)
#' lines(Bstar[,1], Bstar[,2], col = "red")
#' points(Bstar[,1],Bstar[,2], col = "red", pch = 20)

Morris_orientation_matrix <- function(k, p){

stopifnot(is.scalar(k), k >0, k == floor(k),
			is.scalar(p), p >0, p == floor(p))
			
	if( p %% 2 != 0){ 
		p <- ceiling(p / 2) * 2
		warning(sprintf('input p must be even!\n Using p=%g instead of user-defined value\n', p))
		}
		
	m <- k + 1
	Delta <- p / (2 * (p - 1))
	B <- lower.tri(matrix(1, m, k)) * 1
	
	# Create diagonal matrix with (-1, 1) elements
	
	D <- diag(sample(c(-1, 1), k, replace = TRUE))
	
	# Create base value vector
	In <- seq(0, p / 2 - 1) / (p - 1)
	
	tmp <- sample.int(length(In), size = k, replace = TRUE)
	
	x <- In[tmp]
	
	# Create random permutation matrix
	
	P <- diag(k)
	idx <- sample(k, k)
	P <- P[idx,]
	
	Jm1 <- matrix(1, m, 1)
	
	# Create a random orientation of B:
	
	Bstar <- (Jm1 %*% x + Delta / 2 * ((2 * B - 1) %*% D + 1)) %*% P

return(Bstar)
}