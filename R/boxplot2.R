#' Boxplot of vectors with two groups

#' Boxplots when the mean, lower values and upper values are specified for two groups a and b.

#' @param a vector (\code{M}) of mean values to be plotted
#' @param b vector (\code{M}) of mean values to be plotted
#' @param al vector (\code{M}) of lower values to be plotted. Default \code{al = NULL}
#' @param au vector (\code{M}) of upper values to be plotted. Default \code{au = NULL}
#' @param bl vector (\code{M}) of lower values to be plotted. Default \code{bl = NULL}
#' @param bu vector (\code{M}) of upper values to be plotted. Default \code{bu = NULL}
#' @param labels labels for the x-axis of the boxplot
#' @param leg names for the legend. Default \code{leg = NULL}
#' @param col colors for the plot. Default  \code{col = 2:3}
#' @param ylim limits of the y-axis. Default  \code{ylim = NULL}
#' @param ... others arguments to be passed in \code{\link{plot}}

#' @seealso \code{\link{boxplot1}} \code{\link{boxplot}} \code{\link{plot}}

#' @export

#' @examples

#' # Setup the model and define input ranges
#' myfun  <- "ishigami_homma_function"
#' M <- 3
#' DistrFun <- "unif"
#' DistrPar <-  c(-pi, pi)
#' # Sample parameter space using the resampling strategy proposed by 
#' # (Saltelli, 2008; for reference and more details, see help of functions
#' # vbsa_resampling and vbsa_indices) 
#' SampStrategy <- "lhs"
#' N <- 3000 # Base sample size.
#' # Comment: the base sample size N is not the actual number of input 
#' # samples that will be evaluated. In fact, because of the resampling
#' # strategy, the total number of model evaluations to compute the two
#' # variance-based indices is equal to N*(M+2) 
#' X <- AAT_sampling(SampStrategy, M, DistrFun, DistrPar, 2 * N)
#' XABC <- vbsa_resampling(X)
#' # Run the model and compute selected model output at sampled parameter
#' # sets:
#' YA <- model_execution(myfun, XABC$XA) # size (N,1)
#' YB <- model_execution(myfun, XABC$XB) # size (N,1)
#' YC <- model_execution(myfun, XABC$XC) # size (N*M,1)
#' # Compute main (first-order) and total effects:
#' ind <- vbsa_indices(YA, YB, YC)
#' Si <- ind[1,]
#' STi <- ind[2,] 
#' # Plot results:
#' boxplot2(Si, STi, leg = c("main effects", "total effects"))


boxplot2 <- function(a, b, al = NULL, au = NULL, bl = NULL, bu = NULL, labels = NULL, leg = NULL, col = 2:3, ylim = NULL, ...){
	
	stopifnot(length(a) == length(b))
	
	if(is.null(labels)){
		
		if(is.null(names(a))){ 
			names(a) <- paste("#", 1:length(a))
			}
			
		if(is.null(names(b))) names(b) <- names(a)	
			  
	} else {
		
		names(a) <- labels
		names(b) <- labels
	}
	
	if(is.null(al)){
	
	if(is.null(ylim)) ylim <- range(0, 1, a, b)
	
	dat <- data.frame(id2 = letters[1:length(a)], id = rep(c('a', 'b'), each = length(a)), ab = c(rbind(a, b)))
	
	boxplot(ab ~ id2 + id, data = dat, col = col, border = col, ylim = ylim, xaxt = 'n', ... )	

	axis(1, at = seq(1.5, 2 * length(a) + .5, by = 2), names(a))
	legend("topright", legend = leg, col = col, lty = 1, bg = 'white')
	
	} else {
		
		amat <- cbind(a, al, au)
		bmat <- cbind(b, bl, bu)
		
		if(is.null(ylim)) ylim <- range(0, 1, amat, bmat)
		
		rownames(bmat) <- NULL
		rownames(amat) <- NULL
		
		ab <- matrix(nrow = 2* nrow(amat), ncol = ncol(bmat))
		
		ab[seq(1, nrow(ab), by = 2), ] <- amat
		ab[seq(2, nrow(ab), by = 2), ] <- bmat
		
	dat <- data.frame(id2 = letters[1:length(a)], id = rep(c('a', 'b'), each = length(a)), ab = ab)
	boxplot(ab ~ id2 + id, data = dat, col = col, whisklty = 0, staplelty = 0, boxwex = .3, ylim = ylim, xaxt = 'n', ... )
	
	axis(1, at = seq(1.5, 2 * nrow(amat) + .5, by = 2), names(a))
	
	legend("topright", legend = leg, col =col, lty = 1, bg = 'white')
		
	}
	
	
}