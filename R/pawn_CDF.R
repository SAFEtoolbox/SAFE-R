#' Compute the CDFs for the PAWN method
#' 
#' This function computes the unconditional output Cumulative Distribution Funtions 
#' (i.e. when all inputs vary) and the conditional CDFs (when one input is fixed to a given 
#' conditioning interval, while the other inputs vary freely).
#' 
#' The function splits the output sample to create the conditional output by calling internally 
#' the function \link{pawn_split_sample}. 
#' The splitting strategy is an extension of the strategy for uniformy distributed inputs 
#' described in Pianosi and Wagener (2018) to handle inputs sampled from any distribution.
#' (see help of \link{pawn_split_sample} for further explanation).
#' 
#' The sensitivity indices for the PAWN method (KS statistic) measures the distance between 
#' these conditional and unconditional output CDFs 
#' (see help of \link{pawn_indices} for further details and reference).
#' 
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, 1)} set of output samples 
#' @param n number of conditional intervals (default: 10) \cr
#' - integer if all inputs have the same number of groups \cr
#' - list of \code{M} integers otherwise
#' @param dummy boolean, if dummy is True, an articial input is added to the 
#' set of inputs and the sensitivity indices are calculated for the dummy input. 
#' See the help of \link{pawn_indices} for refence and further explanation 
#' on the usage of the dummy input.
#' @param verbose boolean, provides additional information on how the samples are split by the function \link{pawn_split_sample}
#' 
#' @return List containing: 
#' \itemize{
#'   \item \code{YF} vector \code{(P, 1)}, values of Y at which the Cumulative  Distribution Functions (CDFs) FU and FC are given 
#'   \item \code{FU} list of \code{M} elements, values of the empirical unconditional output CDFs.
#'   \item \code{FC} list of lists, values of the empirical conditional output CDFs for each input and each conditioning interval.
#'   \code{FC[[i]]} is a list of \code{n_eff[i]} CDFs conditional to the i-th input.
#'   \code{FC[[i]][k]} is obtained by fixing the i-th input to its k-th conditioning interval (while the other inputs vary freely)
#'   \item \code{xc}  list of \code{M} elements, subsamples' centers (i.e. mean value of the \code{Xi} over each conditioning interval)
#'   \code{xc[[i]]} vector of length \code{n_eff[i]} contains the centers of the \code{n_eff[i]} conditioning intervals for the i-th input.
#'   \item \code{idx_bootsize_min} scalar, used to estimate the smaller sample size across inputs, 
#'   so that the sensitivity index for the dummy input estimates the 'worst' approximation error 
#'   of the sensitivity index across the inputs.
#'   \item \code{FC_dummy} vector \code{(P, 1)}, conditional CDF for the dummy input.
#' }
#'   
#' @export
#'   
#' @references Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity analysis from a generic input-output sample, Env. Mod. & Soft., 108, 197-207. \cr


pawn_CDF <- function(X, Y, n=10, dummy = FALSE, verbose = TRUE){
  
  # This functions splits a generic input-output dataset to create the conditional samples
  # for the approximation of PAWN sensitivity indices
  pawn_samp <- pawn_split_sample(X, Y, n, verbose) 
  XX <- pawn_samp$XX
  YY <- pawn_samp$YY
  xc <- pawn_samp$xc
  NC <- pawn_samp$NC
  Xk <- pawn_samp$Xk
  n_eff <- pawn_samp$n_eff
  
  N <- nrow(X)
  M <- ncol(X)
  
  # Set points at which the CDFs will be evaluated:
  YF <- as.vector(sort(unique(Y)))
  
  # Compute conditional CDFs
  # (bootstrapping is not used to assess conditional CDFs):
  FC <- list(list())

  for(ii in 1:M){ # loop over inputs
    FC[[ii]] <- rep(list(NA),n_eff[ii])
    for(kk in 1:n_eff[ii]){ # loop over conditioning intervals
      FC[[ii]][[kk]] <- ecdf(YY[[ii]][[kk]])(YF)
    }
  }
  
  # Initialize unconditional CDFs:
  FU <- list()
  
  # M unconditional CDFs are computed (one for each input), so that for
  # each input the conditional and unconditional CDFs are computed using the
  # same number of data points (when the number of conditioning intervals
  # n_eff[i] varies across the inputs, so does the shape of the conditional
  # outputs YY[i]).
  
  # Determine the sample size for the unconditional output bootsize:
  bootsize <- as.vector(unlist(lapply(NC,min)))
  # bootsize is equal to the sample size of the conditional outputs NC, or
  # its  minimum value across the conditioning intervals when the sample size
  # varies across conditioning intervals as may happen when values of an
  # input are repeated several times (more details on this in the Note in the
  # help of the function).
  
  # To reduce the computational time (the calculation of empirical CDF is
  # costly), the unconditional CDF is computed only once for all inputs that
  # have the same value of bootsize[ii].
  bootsize_unique <- sort(unique(bootsize))
  N_compute <- length(bootsize_unique) # number of unconditional CDFs that will
  # be computed for each bootstrap resample
  
  # Determine the sample size of the subsample for the dummy input.
  # The sensitivity index for the dummy input will be estimated at this minimum sample size
  # so to estimate the 'worst' approximation error of the sensitivity index
  # across the inputs:
  
  idx_bootsize_min <- NA # In case dummy = FALSE
  
  if(dummy == TRUE){
    bootsize_min <- min(bootsize) # we use the smaller sample size across
    # inputs, so that the sensitivity index for the dummy input estimates
    # the 'worst' approximation error of the sensitivity index across the
    # inputs:
    idx_bootsize_min <- match(bootsize_min,bootsize)
  }
  
  if(N_compute > 1 && verbose){
    warning('The number of data points to estimate the conditional and 
    unconditional output varies across the inputs. The CDFs 
    for the dummy input were computed using the minimum sample 
    size to provide an estimate of the "worst" approximation
    of the sensitivity indices across input.')
  }
  
  for(kk in 1:N_compute){
    # Bootstrap resapling (Extract an unconditional sample of size
    # bootsize_unique[kk] by drawing data points from the full sample Y
    # without replacement
    idx_bootstrap <- sample.int(N, bootsize_unique[kk], replace=FALSE)
    
    # Compute unconditional CDF:
    FUkk <- ecdf(Y[idx_bootstrap])(YF)
    # Associate the FUkk to all inputs that require an unconditional
    # output of size bootsize_unique[kk]:
    idx_input <- which(bootsize == bootsize_unique[kk])
    
    
    for(ii in 1:length(idx_input)){
      FU[[idx_input[ii]]] <- FUkk
    }
  
  }
  
  FC_dummy <- NA # In case dummy = FALSE
  
  if(dummy == TRUE){
    # Bootstrap again from unconditional sample (the size of the
    # resample is equal to bootsize_min):
    idx_dummy <- sample.int(N, bootsize_min, replace=FALSE)
    # Compute empirical CDFs for the dummy input:
    FC_dummy <- ecdf(Y[idx_dummy])(YF)
  }
  
  robj <- list(YF=YF,FU=FU,FC=FC,xc=xc,idx_bootsize_min=idx_bootsize_min,FC_dummy=FC_dummy)
  return(robj)
}  


