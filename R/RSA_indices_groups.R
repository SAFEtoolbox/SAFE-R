#' Regional Sensitivity Analysis with grouping
#'
#' Computation function for Regional Sensitivity Analysis with grouping (as first proposed by Wagener et al. 2001). \cr
#' The function can handle discrete outputs. \cr
#' The function splits the samples in a dataset X into \code{ngroup} sub-sets corresponding to \code{ngroup} of equal size 
#' based on the value of \code{Y} (i.e. 'equiprobable' groups). \cr
#' Then it assesses the distance (i.e. maximum vertical distance called 'mvd' and area between CDFs called 'spread') 
#' between pairs of CDFs of \code{X} in the different sub-sets. It aggregates the values using a statistic (median,
#' mean and maximum) e.g. for mvd the function computes:
#' 
#' \eqn{mvd_median = median( max_x( | Fi(x) - Fj(x) | ) )}
#' 
#' \eqn{mvd_mean = mean( max_x( | Fi(x) - Fj(x) | ) )}
#'
#' or
#' 
#' \eqn{mvd_max = max( max_x( | Fi(x) - Fj(x) | ) )}
#'
#' where \eqn{Fi()} is the CDF of \code{X} in i-th dataset and \eqn{Fj()} is the CDF in the
#' j-th group.
#'
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, P)} set of output samples
#' @param ngroup integer, number of groups considered (default: 10)   
#' @param Nboot scalar, number of resamples used for boostrapping. Default \code{Nboot = 0}, i.e. no bootstrapping.
#' 
#' @return List containing: 
#' \itemize{
#' \item \code{mvd_median} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, median of mvd between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{mvd_mean} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, mean of mvd between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{mvd_max} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, max of mvd between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{spread_median} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, median of spread between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{spread_mean} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, mean of spread between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{spread_max} vector \code{(M)} if Nboot <= 1 or matrix \code{(Nboot,M)} if Nboot > 1, max of spread between pairs of inputs' CDFs estimated for 
#' the different sub-sets and for each bootstrap resample when bootstrapping is used (i.e.Nboot > 1)
#' \item \code{idxb} vector \code{N} respective group of the samples
#' \item \code{Yk} vector \code{ngroup + 1} range of \code{Y} in each group
#'}
#' You can easily derive the n_groups datasets \code{Xi} as: \code{Xi = X[idx == i]}
#'
#' \cr
#' NOTES: \cr
#' - When \code{Y} is discrete and when the number of values taken by \code{Y} (ny) is lower than 
#' the prescribed number of groups (\code{ngroup}), a group is created for each value of \code{Y} 
#' (and therefore the number of groups is set to ny). \cr
#' The function ensures that values of Y that are repeated several times belong to the same group. 
#' This may lead to a final number of group lower than \code{ngroup} and to a different number of 
#' data points across the groups.
#'
#' @references Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, H. V., and Sorooshian, S. (2001): A framework for development and application of hydrological models, Hydrol. Earth Syst. Sci., 5, 13-26.
#' 
#' @seealso \code{\link{RSA_plot_groups}} \code{\link{RSA_indices_thres}}
#' 
#' @export
#' 
#' @examples
#' # See the demo
#' # demo("workflow_rsa_hymod")

RSA_indices_groups <- function(X, Y, ngroup = 10, Nboot = 0) {
  
  ##############
  # Check inputs
  ##############
  
  stopifnot(is.matrix(X), is.numeric(X),
            is.scalar(ngroup), ngroup > 1, ngroup == floor(ngroup),
            is.scalar(Nboot), Nboot >=0, Nboot == floor(Nboot))

  N <- nrow(X)
  M <- ncol(X)
  
  if(!is.matrix(Y)) Y <- matrix(Y, ncol = 1)
  
  stopifnot(is.numeric(Y), nrow(Y) == N)
  
  ########################
  # Compute indices
  ########################
  
  if (Nboot > 1){
    
    B <- matrix(sample.int(N, N * Nboot, replace = TRUE), N, Nboot)
    
    ###### sono arrivata qui, cambiare da qui in poi ####
    all_stat <- apply(B, 2, function(h)  compute_indices_rsa_gr(X[h,], Y[h], ngroup, Nboot))
    
    mvd_median <- t(sapply(all_stat,  function(x) colMedians(x$mvd,na.rm = TRUE)))
    mvd_mean <- t(sapply(all_stat,  function(x) colMeans(x$mvd,na.rm = TRUE)))
    mvd_max <- t(sapply(all_stat,  function(x) colMaxs(x$mvd,na.rm = TRUE)))
    spread_median <- t(sapply(all_stat,  function(x) colMedians(x$spread,na.rm = TRUE)))
    spread_mean <- t(sapply(all_stat,  function(x) colMeans(x$spread,na.rm = TRUE)))
    spread_max <- t(sapply(all_stat,  function(x) colMaxs(x$spread,na.rm = TRUE)))
    
    RSAstat <- compute_indices_rsa_gr(X, Y, ngroup, Nboot = 0)
    idx <- RSAstat$idx
    Yk <- RSAstat$Yk
    
    robj <- list(mvd_median = mvd_median, mvd_mean = mvd_mean, mvd_max = mvd_max,
                 spread_median = spread_median, spread_mean = spread_mean, spread_max = spread_max, 
                 idx = idx, Yk = Yk)
    
  } else {
    
    RSAstat <- compute_indices_rsa_gr(X, Y, ngroup, Nboot = 0)
    mvd <- RSAstat$mvd
    spread <- RSAstat$spread
    idx <- RSAstat$idx
    Yk <- RSAstat$Yk
    
    # Calculate statistics across pairs of CDFs:
    mvd_median <- colMedians(mvd,na.rm = TRUE)
    mvd_mean <- colMeans(mvd,na.rm = TRUE)
    mvd_max <- colMaxs(mvd,na.rm = TRUE)
    spread_median <- colMedians(spread,na.rm = TRUE)
    spread_mean <- colMeans(spread,na.rm = TRUE)
    spread_max <- colMaxs(spread,na.rm = TRUE)

    robj <- list(mvd_median = mvd_median, mvd_mean = mvd_mean, mvd_max = mvd_max,
                 spread_median = spread_median, spread_mean = spread_mean, spread_max = spread_max, 
                 idx = idx, Yk = Yk)
    
  }
  
  return(robj)
}