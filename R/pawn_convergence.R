#' This function computes the PAWN sensitivity indices  
#' using sub-samples of the original sample \code{Y} of increasing size.
#'
#' The function splits the output sample to create the conditional output.
#' The splitting strategy is an extension of the strategy for uniformy
#' distributed inputs described in Pianosi and Wagener (2018) to handle inputs sampled from 
#' any distribution.
#' 
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, 1)} set of output samples 
#' @param n number of conditional intervals (default: 10) \cr
#' - integer if all inputs have the same number of groups \cr
#' - list of \code{M} integers otherwise
#' @param NN vector \code{(R)} of subsample sizes at which indices will be estimated (\code{max(NN)} must not exceed \code{N})
#' @param Nboot scalar, number of bootstrap resamples to derive confidence intervals
#' @param dummy boolean, if \code{dummy} is TRUE, an artificial input is added to the set of inputs 
#' and the sensitivity indices are calculated for the dummy input.
#' See help of \link{pawn_indices} for reference and further explanation on the usage of the dummy input. \cr
#' ADVANCED USAGE: \cr
#' for Regional-Response Global Sensitivity Analysis: \cr
#' @param output_condition function (see the help for \link{pawn_indices} for information on this optional input)
#' @param par list (see the help for \link{pawn_indices} for information on this optional input)
#'  
#' @return List containing:
#' \itemize{
#'   \item \code{KS_median} list of shape \code{(Nboot, M)} median KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_mean} list of shape \code{(Nboot, M)} mean KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_max} list of shape \code{(Nboot, M)} max KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_dummy} list of shape \code{(Nboot, 1)} KS of dummy input (one value for each bootstrap resample)
#'   }
#'   @seealso \code{\link{pawn_split_sample}} \code{\link{pawn_indices}}
#'   
#' @export
#'   
#' @references Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity analysis from a generic input-output sample, Env. Mod. & Soft., 108, 197-207. \cr
#' Pianosi, F. and Wagener, T. (2015), A simple and efficient method for global sensitivity analysis based on cumulative distribution functions, Env. Mod. & Soft., 67, 1-11. \cr
#' Khorashadi Zadeh et al. (2017), Comparison of variance-based and moment-independent global sensitivity analysis approaches by application to the SWAT model, Environmental Modelling & Software,91, 210-222.

pawn_convergence <- function(X, Y, n=10, NN, Nboot = 1, dummy = FALSE, output_condition = allrange, par = list()){
  
  stopifnot(is.matrix(X), is.numeric(X), is.numeric(Y), is.scalar(Nboot), Nboot >=1, Nboot == floor(Nboot), 
            is.logical(dummy))
  
  stopifnot(NN >= 0, NN - floor(NN) < .1^6, diff(NN) > 0, max(NN) <= N)
  if(is.matrix(Y) == FALSE){Y <- as.matrix(Y)}
  
  N <- nrow(X)
  M <- ncol(X)
  
  ###########################
  # Compute indices
  ###########################
  
  R1 <- length(NN)
  R <- ifelse(max(NN) == N, R1 - 1, R1)
  
  # Initialise variables
  KS_median <- list(list()) #matrix(NA,R1,M)
  KS_mean <- list(list()) #matrix(NA,R1,M)
  KS_max <- list(list()) # matrix(NA,R1,M)
  KS_dummy <- NA # In case dummy = FALSE
  if(dummy == TRUE){
    KS_dummy <- list() #matrix(NA,R1,1)
  }
  
  pawn_ind <- list(list())
  
  for(jj in 1:R){ # loop over sample sizes
    
    idx_new <- shrink.lhcube(X, NN[jj], type = "indices")
    Xj <- X[idx_new, , drop = FALSE]
    # drop rows while trying to maximise the spread between the points
    Yj <- Y[idx_new, ]
    
    pawn_ind[[jj]] <- pawn_indices(Xj, Yj, n, Nboot, dummy, output_condition, par)
    KS_median[[jj]] <- pawn_ind[[jj]]$KS_median
    KS_mean[[jj]] <- pawn_ind[[jj]]$KS_mean
    KS_max[[jj]] <- pawn_ind[[jj]]$KS_max
    
    if(dummy == TRUE){
      KS_dummy[[jj]] <- pawn_ind[[jj]]$KS_dummy
    }
    
  }
  
  if(max(NN) == N){
    pawn_ind_N <- pawn_indices(X, Y, n, Nboot, dummy, output_condition, par)
    KS_median <- c(KS_median, list(pawn_ind_N$KS_median))
    KS_mean <- c(KS_mean, list(pawn_ind_N$KS_mean))
    KS_max <- c(KS_max, list(pawn_ind_N$KS_max))
    
    if(dummy == TRUE){
      pawn_dummy <- pawn_ind_N$KS_dummy
      KS_dummy <- c(KS_dummy, list(pawn_ind_N$KS_dummy))
    }
  }
  
  robj <- (list(KS_median=KS_median, KS_mean=KS_mean, KS_max=KS_max, KS_dummy=KS_dummy))
  return(robj)
  
} 





