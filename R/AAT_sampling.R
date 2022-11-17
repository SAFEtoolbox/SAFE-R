#' Random sampling in the M-dimension space of the input factors of a model.
#'
#' This function performs random sampling in the \code{M}-dimension space of the input factors of a model. To this purpose, the function first perform uniform sampling in the unit hypercube in R^\code{M}; then, project the sampled points along each direction.
#'
#' @param samp_strat string, sampling strategy. Options: \code{"rsu"}: random uniform, \code{"lhs"}: latin hypercube.                            
#' @param M positive integer, number of inputs. 
#' @param  distr_fun string (eg: \code{"unif"}) if all inputs have the same pdf, or vector of \code{M} strings (eg: \code{c("unif", "norm")}) containing the probability distribution function of each input.
#' @param distr_par row vector if all input pdfs have the same parameters, list of \code{M} vectors otherwise containing the parameters of the probability distribution function.
#' @param N positive integer, number of samples.

#' @return \code{X} matrix (\code{N}, \code{M}) of samples. Each row is a point in the input space. In contrast to \code{\link{OAT_sampling}}, rows do not follow any specific order, and all components (columns) differ from point (row) to point.   

#' @seealso \code{\link{OAT_sampling}}

#' @export
#' @examples
#' # Example 1: 2 inputs, both from Unif[0,3]
#' N <- 1000 
#' M <- 2 
#' distr_fun <- "unif"
#' distr_par <- c(0, 3)
#' samp_strat <- "lhs"
#' X <- AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
#' # Plot results:
#' plot(X, xlab = expression(x[1]), ylab = expression(x[2]))
#' 
#' # Example 2: 2 inputs, one from Unif[0,3], one from Unif[1,5]
#' distr_fun <- "unif"
#' distr_par <- list(c(0, 3), c(1, 5))
#' X <- AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
#' # Plot results:
#' plot(X, xlab = expression(x[1]), ylab = expression(x[2]))
#'
#' # Example 3: 2 inputs, one from Unif[0,3], one from discrete Unifd[1,5]
#' distr_fun <- c("unif","unifd")
#' distr_par <- list(c(0, 3), c(1, 5))
#' X <- AAT_sampling(samp_strat, M, distr_fun, distr_par, N)
#' 
#' 
#'##################################
#'# Some more comments about sampling
#'##################################
#'
#' # If you want to see the difference between Latin-Hypercube, Sobol and
#' # Monte Carlo sampling strategy, you can use the following code:
#'
#' N <- 100
#' dist_type  <- "unif"
#' dist_param <- c(0, 1)
#' X1 <- AAT_sampling("rsu", 2, dist_type, dist_param, N)
#' X2 <- AAT_sampling("lhs", 2, dist_type, dist_param, N)
#' par(mfrow = c(1, 2))
#' plot(X1, xlab = expression(x[1]), ylab = expression(x[2]), main = ("Random Uniform"))
#' plot(X2, xlab = expression(x[1]), ylab = expression(x[2]), main = ("Latin Hypercube"))
#' # Example with a Normal distribution on the vertical axis:
#' dist_type  <- c("unif", "norm")
#' dist_param <- list(c(0, 1), c(5, 1))
#' X1 <- AAT_sampling("rsu", 2, dist_type, dist_param, N)
#' X2 <- AAT_sampling("lhs", 2, dist_type, dist_param, N)
#' par(mfrow = c(1, 2))
#' plot(X1, xlab = expression(x[1]), ylab = expression(x[2]), main = ("Random Uniform"))
#' plot(X2, xlab = expression(x[1]), ylab = expression(x[2]), main = ("Latin Hypercube"))

AAT_sampling <- function(samp_strat, M, distr_fun, distr_par, N){

##############
# Check inputs
##############

		stopifnot(samp_strat %in% c('rsu', 'lhs'),
			is.scalar(M), M >0, M == floor(M),
			length(distr_fun) %in% c(1, M),
			distr_fun %in% c('chisq', 'exp', 'geom', 'pois', 't',
				'beta', 'binom', 'f', 'gamma', 'lnorm', 'nbinom', 'norm', 'unif', 'weibull',
				'hyper', 'unifd'),
			is.scalar(N), N >0, N == floor(N)
		)
		
	if(!is.list(distr_par) & length(distr_fun) == 1){
		distr_par <- lapply(vector('list', M), function(h) distr_par)
	}
	
	stopifnot(is.list(distr_par), length(distr_par) == M)

	
##################
# Perform sampling
##################

# Uniform sampling in the unit square:

if(samp_strat == 'rsu'){      # Uniform sampling
	X <- matrix(runif(N * M), nrow = N, ncol = M)
    } else { #(samp_strat == 'lhs') # Latin Hypercube sampling
 	X   <- lhcube(N, M)
    }
    
    row.names(X) <- paste('sample', 1:N, sep ="")
    
    if(length(distr_fun == 1)){
    	distr_fun <- rep(distr_fun, M)
    }
    	
    X <- sapply(1:M, function(i) do.call(sprintf("q%s", distr_fun[i]), c(list(X[, i]), as.list(distr_par[[i]])) ))
	
	return(X)
}