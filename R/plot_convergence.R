#' Plot convergence

#' Draws the plot to analise the convergence

#' @param n vector \code{P} of number of model evaluations
#' @param est matrix \code{(P, M)} with the mean of the indices. 
#' @param estl matrix \code{(P, M)} with the lower bound of the estimates of the indices. Default \code{estl = NULL}
#' @param estu matrix \code{(P, M)} with the lower bound of the estimates  of the indices. Default \code{estu = NULL}
#' @param ex exact values to be drawn as a line. Default \code{ex = NULL}
#' @param labels labels for legend. Default \code{labels = NULL}
#' @param ylim limits of the y-axis. Default  \code{ylim = NULL}
#' @param col colors for the plot. Default  \code{col = 1:M}
#' @param lty line type. Default  \code{lty = 1:(M-1)}
#' @param ... others arguments to be passed in \code{\link{plot}}

#' @export

#'@examples
#' # Setup the model and define input ranges

#' fun_test  <- "ishigami_homma_function"
#' M <- 3
#' distr_fun <- "unif"
#' distrpar <-  c(-pi, pi)

#' # Compute the exact values of the output variance (V) and of 
#' # the first-order (Si_ex) and total-order (STi_ex) 
#' # variance-based sensitivity indices (this is possible in this
#' # very specific case because V, Si_ex and STi_ex can be computed
#' # analytically)

#' ihfun <- ishigami_homma_function(runif(M))

#' Si_ex <- attributes(ihfun)$Si_ex
#' STi_ex <- attributes(ihfun)$STi_ex

#' # Sample parameter space:

#' SampStrategy <- "lhs"
#' N <- 3000
#' X <- AAT_sampling(SampStrategy, M, distr_fun, distrpar, 2 * N)

#' # Apply resampling strategy for the efficient approximation of the indices:
#' XABC <- vbsa_resampling(X)

#' # Run the model and compute selected model output at sampled parameter
#' # sets:

#' YA <- model_evaluation(fun_test, XABC$XA) # size (N,1)
#' YB <- model_evaluation(fun_test, XABC$XB) # size (N,1)
#' YC <- model_evaluation(fun_test, XABC$XC) # size (N*M,1)

#' # Analyze convergence of sensitivity indices:
#' NN <- seq(N / 5, N, by = N/5)

#' conv <- vbsa_convergence(c(YA, YB, YC), M, NN)

#' Si <- conv$Si
#' STi <- conv$STi

#' par(mfrow = c(1, 2))

#' plot_convergence(NN * (M + 2), Si, ex = Si_ex, 
#' xlab = "model evals", ylab = "main effect")
#' plot_convergence(NN * (M + 2), STi, ex = STi_ex, 
#' xlab = "model evals", ylab = "total effect")

#' # With bootstrap
#' Nboot <- 100

#' conv100 <- vbsa_convergence(c(YA, YB, YC), M, NN, Nboot)

#' Si100 <- conv100$Si
#' STi100 <- conv100$STi

#' Sil100 <- conv100$Si_lb
#' Siu100 <- conv100$Si_ub

#' STil100 <- conv100$STi_lb
#' STiu100 <- conv100$STi_ub

#' par(mfrow = c(1, 2))

#' plot_convergence(NN * (M + 2), Si100, Sil100, Siu100, ex = Si_ex, 
#' xlab = "model evals", ylab = "main effect")
#' plot_convergence(NN * (M + 2), STi100, STil100, STiu100, ex = STi_ex, 
#' xlab = "model evals", ylab = "total effect")



plot_convergence <- function(n, est, estl = NULL, estu = NULL, ex = NULL, labels = NULL, ylim =  NULL, col = 1:M, lty = 1:(M-1),  ...){
	
	M <- ncol(est)
	
	if(is.null(labels)){
		
		if(is.null(colnames(est))){ 
			colnames(est) <- paste("#", 1:M)
			}  
	} else {
		
		colnames(est) <- labels
	}
	
	if(is.null(estl)){
		
	if(is.null(ylim)) ylim <- range(0, 1, est)	
	
	matplot(n, est, type = 'l', xaxt = 'n', lwd = 2, ylim = ylim, lty = lty, col = col, ...)
	
	axis(1, at = n, n)
	
	if(!is.null(ex)) abline(h = ex, col = adjustcolor(col, alpha.f = .2), lty = 2)
	
	legend('bottomleft', legend = colnames(est), lty = lty, col = col, bg = 'white', lwd = 2)
	
	} else {
		
		if(is.null(ylim)) ylim <- range(0, 1, est, estl, estu)	
		
		matplot(n, est, type = 'n', xaxt = 'n', lty = 2, ylim = ylim, ...)
		axis(1, at = n, n)
		
		for (i in 1:M){
			polygon(c(n, rev(n)), c(estu[,i], rev(estl[,i])), col = adjustcolor(col[i], alpha.f = .1), border = NA)
			}
			
		if(!is.null(ex)) abline(h = ex, col = adjustcolor(col, alpha.f = .2), lty = 2)	
		
		matlines(n, est, lwd = 2, lty = lty, col = col)
		
		legend('bottomleft', legend = colnames(est), lty = lty, col = col, bg = 'white', lwd = 2)
		
	}
	
	
}