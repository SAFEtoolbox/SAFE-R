#' EET plot

#' Plot the mean and the standard deviation of the elementary effects.

#' @param m vector (\code{M}) mean of the elementary effects
#' @param s vector (\code{M}) standard deviation of the elementary effects
#' @param ml vector (\code{M}) lower bound of \code{mi} (at level \code{alfa}) across \code{Nboot} estimations. Default \code{ml = NULL}
#' @param mu vector (\code{M}) upper bound of \code{mi} (at level \code{alfa}) across \code{Nboot} estimations. Default \code{mu = NULL}
#' @param sl vector (\code{M}) lower bound of \code{sigma} across \code{Nboot} estimations. Default \code{sl = NULL}
#' @param su vector (\code{M}) upper bound of \code{sigma} across \code{Nboot} estimations. Default \code{su = NULL}
#' @param labels labels for legend. Default \code{labels = NULL}
#' @param xlim limits of the x-axis. Default \code{xlim = NULL}
#' @param ylim limits of the y-axis. Default  \code{ylim = NULL}
#' @param xlab label for the x-axis. Default \code{xlab = 'Mean of EEs'}
#' @param ylab label for the x-axis. Default \code{ylab = 'Sd of EEs'}
#' @param col colors for the plot. Default  \code{col = 1:length(m)}
#' @param pch symbols to use when plotting points. Default  \code{pch = 15:19}
#' @param ... others arguments to be passed in \code{\link{plot}}

#' @seealso \code{\link{EET_indices}} \code{\link{EET_convergence}} \code{\link{plot}}

#' @export

#' @examples

#' # See the demo
#' # demo("workflow_eet_hymod")
#' # or
#' # demo("workflow_eet_hbv")

EET_plot <- function(m, s, ml = NULL, mu = NULL, sl = NULL, su = NULL, labels = NULL, xlim = NULL, ylim = NULL,  xlab = 'Mean of EEs', ylab = 'Sd of EEs', col = 1:length(m), pch = 15:19,  ...){
	
	stopifnot(length(m) == length(s))
	
	if(is.null(labels)){
		
		if(is.null(names(m))){ 
			names(m) <- paste("#", 1:length(m))
			}
			
		if(is.null(names(s))) names(s) <- names(m)	
			  
	} else {
		
		names(m) <- labels
		names(s) <- labels
	}
	
	if(is.null(ml)){
	
	plot(m, s, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
	legend('topleft', names(s), col = col, pch = pch)
	
	} else {
		
		plot(0, type ='n', xlim = range(ml, mu), ylim = range(sl, su),  xlab = xlab, ylab = ylab, ...)
		
		for (i in 1:length(m)){
			polygon(c(ml[i], mu[i], mu[i], ml[i]), rep(c(sl[i], su[i]), each = 2), col = adjustcolor(i, alpha.f = .1), border = NA)
			}
	
	points(m, s, col = col, pch = pch)
	
	legend('topleft', names(s), col = col, pch = pch)

		
	}
	
	
}