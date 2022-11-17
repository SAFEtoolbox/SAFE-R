#' Mean and standard deviation of Elementary Effects
#'
#' This function computes mean and standard deviation of Elementary Effects (EEs) \cr
#' Compute the sensitivity indices according to the Elementary Effects Test
#' (Saltelli, 2008) or 'method of Morris' (Morris, 1991). \cr
#' These are: the mean (\code{mi}) of the Elementary Effects (EEs) associated to
#' input 'i', which measures the input influence; and the standard deviation
#' (\code{sigma}) of the EEs, which measures its level of interactions with other inputs.
#' For the mean EE, we use the version suggested by Campolongo et al. (2007), 
#' where absolute values of the EEs are used (this is to avoid that EEs with opposite sign 
#' would cancel each other out).
#'
#' @param r scalar, number of sampling point
#' @param xrange list of \code{M} vectors of length 2 containing the input ranges
#' @param X matrix (\code{r * (M + 1)}, \code{M}) matrix of sampling datapoints where EE must be computed  
#' @param Y vector (\code{r * (M + 1)}) vector of output values
#' @param design_type  string, design type (string) Options: \code{"trajectory"}, \code{"radial"}
#' @param Nboot scalar, number of resamples used for boostrapping (default:0)    - scalar
#' @param alfa scalar, significance level for the confidence intervals estimated by bootstrapping (default: 0.05)
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{mi} vector (\code{M}) mean of the elementary effects
#'   \item \code{sigma} vector (\code{M}) standard deviation of the elementary effects
#'   \item \code{EE} matrix (\code{r}, \code{M}) of elementary effects
#' }
#' If \code{Nboot > 1} it also contains
#' \itemize{
#'   \item \code{mi_sd} vector (\code{M}) standard deviation of \code{mi} across \code{Nboot} estimations
#'   \item \code{sigma_sd} vector (\code{M}) standard deviation of \code{sigma} across \code{Nboot} estimations
#'   \item \code{mi_lb} vector (\code{M}) lower bound of \code{mi} (at level \code{alfa}) across \code{Nboot} estimations
#'   \item \code{sigma_lb} vector (\code{M}) lower bound of \code{sigma} across \code{Nboot} estimations
#'   \item \code{mi_ub} vector (\code{M}) upper bound of \code{mi} (at level \code{alfa}) across \code{Nboot} estimations
#'   \item \code{sigma_ub} vector (\code{M}) upper bound of \code{sigma} across \code{Nboot} estimations. 
#' }
#'
#' @seealso \code{\link{EET_convergence}} \code{\link{EET_plot}}
#'
#' @export
#' 
#' @references Morris, M.D. (1991), Factorial sampling plans for preliminary computational experiments, Technometrics, 33(2). \cr
#' Saltelli, A., et al. (2008), Global Sensitivity Analysis, The Primer, Wiley. \cr
#' Campolongo, F., Cariboni, J., Saltelli, A. (2007), An effective screening design for sensitivity analysis of large models. Environ. Model. Softw. 22 (10), 1509-1518.
#'
#' @examples

#' # See the demo
#' # demo("workflow_eet_hymod")
#' # or
#' # demo("workflow_eet_hbv")


EET_indices <- function(r, xrange, X, Y, design_type, Nboot = 0, alfa = 0.05)
{
#########
# Check inputs
#########

stopifnot(is.scalar(r), r >=0, r == floor(r))

n <- nrow(X)
M <- ncol(X)

Dr <- sapply(xrange, diff)

stopifnot(n == r * (M+1), nrow(Y) == n, length(Y) == n, Dr > 0, length(Dr) == M, design_type %in% c("radial","trajectory"))

###################################
# Check optional inputs
###################################

stopifnot(is.scalar(Nboot), Nboot >=0, Nboot == floor(Nboot),
is.numeric(alfa), alfa <= 1, alfa >= 0)

#######################
# Compute Elementary Effects
#######################

EE <- matrix(NA,r,M) # matrix of elementary effects
k  <- 1
ki <- 1

for (i in 1:r){
    for (j in 1:M ){
        if (design_type == 'radial'){ 
        		# radial design: EE is the difference 
            # between output at one point in the i-th block and output at
            # the 1st point in the block
            EE[i,j] <- abs( Y[k + 1] - Y[ki] ) / abs( X[k+1, j] - X[ki, j] ) * Dr[j]
            
            } else { # if (design_type == 'trajectory') 
            	
            # trajectory design: EE is the difference 
            # between output at one point in the i-th block and output at
            # the previous point in the block (the "block" is indeed a
            # trajectory in the input space composed of points that
            # differ in one component at the time)
            idx <- which( abs( X[k+1, ] - X[k,] ) > 0 )  # if using 'morris' 
            # sampling, the points in the block may not
            # be in the proper order, i.e. each point in the block differs
            # from the previous/next one by one component but we don't know
            # which one; this is here computed and saved in 'idx'
            
            EE[i, idx] <- abs( Y[k + 1] - Y[k] ) / abs( X[k+1, idx] - X[k, idx] ) * Dr[idx]
        }
         
        k = k + 1
    }
    k = k + 1
    ki = k
}

#######################
# Compute Mean and Standard deviation
#######################


if (Nboot > 1){
    
    B <- matrix(sample.int(r, r * Nboot, replace = TRUE), r, Nboot)

	mi_all <- apply(B, 2, function(b) colMeans(EE[b,]))
	sigma_all <- apply(B, 2, function(b) apply(EE[b,], 2, sd))

    mi <- rowMeans(mi_all)
    mi_sd <- apply(mi_all, 1, sd)
    mi_lb <- apply(mi_all, 1, sort) 
    mi_lb <- mi_lb[max(1,round(Nboot * alfa / 2)),]
    mi_ub <- apply(mi_all, 1, sort)  
    mi_ub <- mi_ub[round(Nboot *(1 - alfa / 2)),]

    sigma <- rowMeans(sigma_all)
    sigma_sd <- apply(sigma_all, 1, sd)
    sigma_lb <- apply(sigma_all, 1, sort)
    sigma_lb <- sigma_lb[max(1,round(Nboot * alfa / 2)),]
    sigma_ub <- apply(sigma_all, 1, sort) 
    sigma_ub  <- sigma_ub[round(Nboot * (1 - alfa / 2)),]

     cat(sprintf('\n\t mean(EE) std(EE)\n'))
     cat(sprintf('X%d:\t %2.3f\t %2.3f\n', 1:M, mi, sigma))
     
     robj <- list(mi = mi,sigma = sigma, EE = EE, mi_sd = mi_sd, sigma_sd = sigma_sd, mi_lb = mi_lb, sigma_lb = sigma_lb, mi_ub = mi_ub, sigma_ub = sigma_ub, mi_all = mi_all, sigma_all = sigma_all)
     
     } else {

    mi <- colMeans(EE)
    sigma <- apply(EE,2, sd)

robj <- list(mi = mi,sigma = sigma, EE = EE)
}

     cat(sprintf('\n\t mean(EE) std(EE)\n'))
     cat(sprintf('X%d:\t %2.3f\t %2.3f\n', 1:M, mi, sigma))

return(robj)

}
