#'  Variance-based first-order and total effects indices
#'
#' This function computes the variance-based first-order indices
#' (or main effects) and total effects indices (Homma and Saltelli, 1996)
#' using an increasing number of output samples.
#'
#' @param Y vector of \code{N * (M + 2)} model output samples  \code{Y = c(YA, YB, YC)}
#' @param M scalar number of inputs 
#' @param NN vector \code{(R)}, sample sizes at which indices will be estimated (must be a vector of integer positive values not exceeding \code{N})
#' @param Nboot scalar, number of resamples used for boostrapping (default: 0) 
#' @param alfa scalar, significance level for the confidence intervals estimated by bootstrapping (default: 0.05)
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{Si} estimates of the main effects at different sampling size
#'   \item \code{STi}estimates of the total effects at different sampling size. 
#' }
#' If \code{Nboot > 1} it also contains
#' \itemize{
#'   \item \code{Si_sd} standard deviation of main effects at different sampling size
#'   \item \code{STi_sd} standard deviation of total effects at different sampling size
#'   \item \code{Si_lb} lower bound of main effects at different sampling size
#'   \item \code{STi_lb} lower bound of total effects at different sampling size
#'   \item \code{Si_ub} upper bound of main effects at different sampling size
#'   \item \code{STi_ub} upper bound of total effects at different sampling size. 
#' }
#' All output arguments are matrices of size \code{(R, M)}

#' @seealso \code{\link{vbsa_indices}} \code{\link{vbsa_resampling}}

#' @export

#' @examples
#'
#' fun_test  <- "ishigami_homma_function"
#' M <- 3 
#' distrfun <- "unif"
#' distrpar <- c(-pi, pi)
#' N <- 1000
#' sampstrat <- "lhs"
#' X <- AAT_sampling(sampstrat, M, distrfun, distrpar, 2 * N)
#' XABC <- vbsa_resampling(X)
#' YA <- model_execution(fun_test, XABC$XA)
#' YB <- model_execution(fun_test, XABC$XB)
#' YC <- model_execution(fun_test, XABC$XC)
#' Y <- c(YA, YB, YC)
#' NN <- seq(N/10, N, by = N / 10)
#' SiSTi <- vbsa_convergence(Y, M, NN)

vbsa_convergence <- function(Y, M, NN, Nboot = 0, alfa = 0.05){

##############
# Check inputs
##############

if(is.matrix(Y)) Y <- c(Y)

N <- length(Y)

 stopifnot(is.scalar(M), M >=0, M == floor(M), 
 NN >0, NN == floor(NN), diff(NN) > 0, max(NN) <= N,
 is.numeric(Y), length(Y) %% (M+2) == 0
 )

N <- N / (M + 2)

###################################
# Check optional inputs
###################################

stopifnot(is.scalar(Nboot), Nboot >=0, Nboot == floor(Nboot),
is.numeric(alfa), alfa <= 1, alfa >= 0)

#################
# Compute indices
#################

YA <- Y[1:N]
YB <- Y[(N + 1):(N + N)]
YC <- matrix(Y[-(1:(N + N))], N, M)

# the function vbsa_indices_ran allows to compute the indices on a random
# subset of nn samples, nn varies from 1 to NN.
# This function is called in the next row (vbsa_ind).
vbsa_indices_ran <- function(nn, Nboot, alfa){
    indx <- sample.int(N, nn) # randomly select a subset of nn samples
    vbsa_indices(YA[indx], YB[indx], YC[indx,], Nboot = Nboot, alfa = alfa)
}

vbsa_ind <- sapply(NN, vbsa_indices_ran, Nboot = Nboot, alfa = alfa)

if(Nboot == 0){

Si <- t(vbsa_ind[seq(1, 2 * M, by = 2), ])
STi <- t(vbsa_ind[seq(2, 2 * M, by = 2), ])

colnames(Si) <- paste('X', 1:M, sep = "")
rownames(Si) <- NN 

colnames(STi) <- paste('X', 1:M, sep = "")
rownames(STi) <- NN 

robj <- list(Si = Si, STi = STi)

} else {
	
	Si <- t(vbsa_ind[seq(1, 8 * M, by = 8), ])
	Si_sd <- t(vbsa_ind[seq(2, 8 * M, by = 8), ]) 
	Si_lb <- t(vbsa_ind[seq(3, 8 * M, by = 8), ])
	Si_ub <- t(vbsa_ind[seq(4, 8 * M, by = 8), ])
	
	STi <- t(vbsa_ind[seq(5, 8 * M, by = 8), ])
	STi_sd <- t(vbsa_ind[seq(6, 8 * M, by = 8), ])
	STi_lb <- t(vbsa_ind[seq(7, 8 * M, by = 8), ])
	STi_ub <- t(vbsa_ind[seq(8, 8 * M, by = 8), ])
	
	colnames(Si) <- paste('X', 1:M, sep = "")
	rownames(Si) <- NN 
	
	colnames(STi) <- paste('X', 1:M, sep = "")
	rownames(STi) <- NN 

robj <- list(Si = Si, Si_sd = Si_sd, Si_lb = Si_lb, Si_ub = Si_ub, STi = STi, STi_sd = STi_sd, STi_lb = STi_lb, STi_ub = STi_ub)
	
}


return(robj)
}