#' Compute the Kolmogorov-Smirnov (KS) 
#' 
#' This function computes the KS between the unconditional and conditional output 
#' for each input and each conditioning interval, for ONE sample/bootstrap resample. \cr
#' This function is called internally in \link{pawn_indices} and \link{pawn_plot_ks}. \cr
#' The inputs necessary to call this function are estimated by \link{pawn_CDF}.
#' 
#' @param YF vector \code{(P, 1)}, values of Y at which the Cumulative  Distribution Functions (CDFs) FU and FC are given
#' @param FU list of \code{M} elements, values of the empirical unconditional output CDFs.
#' @param FC list of lists, values of the empirical conditional output CDFs for each input and each conditioning interval.
#' \code{FC[[i]]} is a list of \code{n_eff[i]} CDFs conditional to the i-th input.
#' \code{FC[[i]][k]} is obtained by fixing the i-th input to its k-th conditioning interval (while the other inputs vary freely)
#' 
#' ADVANCED USAGE: \cr
#' for Regional-Response Global Sensitivity Analysis: \cr
#' @param output_condition function (see the help for \link{pawn_indices} for information on this optional input)
#' @param par list (see the help for \link{pawn_indices} for information on this optional input)
#' 
#' @return \code{KS} list of \code{M} elements. KS-statistic calculated between conditional and unconditional output for the M inputs 
#' and the \code{n_eff} conditioning intervals. \code{KS[i]} contains the KS values for the i-th input and the 
#' \code{n_eff[i]} conditioning intervals.
#' 
#' @export

pawn_ks <- function(YF, FU, FC, output_condition = allrange, par = list()){
  
  stopifnot(is.function(output_condition), is.scalar(par) || is.list(par))
  
  M <- length(FC) # number of inputs
  
  # Initialise variable
  KS <- vector(mode = "list", length = M)
  
  # Find subset of output values satisfying a given condition
  idx <- output_condition(YF,par)
  
  # Calculate KS statistics:
  for(ii in 1:M){ # loop over inputs
    n_effi <- length(FC[[ii]])
    KS[ii] <- lapply(M,function(x) rep(NA,n_effi))
    
    for(kk in 1:n_effi){  # loop over conditioning intervals
      # Compute KS:
      KS[[ii]][[kk]] <- max(abs(FU[[ii]][idx] - FC[[ii]][[kk]][idx]))
    }
  }
  
  return(KS)

  }
  
  
  