#' PAWN split sample
#' 
#' This function splits a generic input-output dataset in order to create unconditional and conditional samples for the approximation of PAWN sensitivity indices (Pianosi and Wagener, 2018) \cr
#' This function extends the splitting strategy described in (Pianosi and Wagener, 2018), which consists of 
#' splitting the input and output sample into a number of equally spaced conditioning intervals (intervals with the same size) 
#' based on the value of each input \code{Xi}. Here, the conditioning intervals are equiprobable, i.e. there is
#' (approximately) the same number of data points in each interval. While equally spaced intervals are adapted for 
#' inputs that have a uniform distribution, equiprobable intervals allow to handle inputs that have any distribution. \cr
#' This function is called internally in \link{pawn_indices} and \link{pawn_plot_CDF}
#' 
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, 1)} set of output samples 
#' @param n number of conditional intervals (default: 10) \cr
#' - integer if all inputs have the same number of groups \cr
#' - list of \code{M} integers otherwise
#' @param verbose boolean, provides additional information on how the samples are split 
#' 
#' @return List containing: 
#' \itemize{
#'   \item \code{YY} list of \code{M} elements. Output samples to assess the conditional CDFs. 
#'   \code{YY[[i]]} is a list of \code{n_eff[i]} subsamples that can be used to assess \code{n_eff[i]}
#'   conditional distributions with respect to the i-th input.
#'   \code{YY[[i]][k]} is obtained by fixing the i-th input to its k-th conditioning interval 
#'   (while the other inputs vary freely), and its length is \code{NC[[i]][k]}
#'   \item \code{xc} list of \code{M} elements. Subsamples' centers (i.e. mean value of \code{Xi} over each conditioning interval).
#'   The vector \code{xc[[i]]} of dimensions \code{(n_eff[i],1)} contains the subsamples' centers of the i-th input \code{Xi}.
#'   \item \code{NC} list of \code{M} elements. Number of data points in each conditioning interval and for each input.
#'   \code{NC[[i]]} has length \code{n_eff[i]} and contains the sample sizes for the i-th input.
#'   \item \code{n_eff} list of \code{M} elements. Number of conditioning intervals actually used 
#'   for each inputs \code{Xi} (see (*) for further explanation).
#'   \item \code{Xk} list of \code{M} elements. Subsamples' edges for each input \code{Xi} 
#'   (i.e. bounds of \code{Xi} over each conditioning interval)
#'   \code{Xk[i]} is a list of length \code{n_eff[i]+1} and contains the edges of the conditioning intervals 
#'   for the i-th input.
#'   \item \code{XX} list of \code{M} elements. Conditional input samples corresponding to the samples in \code{YY}.
#'   \code{XX[[i]]} is a list of \code{n_eff[i]} subsamples for the i-th input.
#'   \code{XX[[i]][k]} is obtained by fixing the i-th input to its k-th conditioning interval (while the other inputs vary freely), 
#'   and it has shape \code{(NC[[i]][k],M])}
#'   } 
#'   NOTES: \cr
#' (*)  - When \code{Xi} is discrete and when the number of values taken by \code{Xi} \code{(nxi)} is
#' lower than the prescribed number of conditioning intervals \code{(n[i])}, 
#' a conditioning interval is created for each value of \code{Xi} (and therefore the actual number 
#' of conditioning intervals is set to \code{n_eff[i] = nxi}). \cr
#' - The function ensures that values of \code{Xi} that are repeated several times
#' belong to the same group. This may lead to a number of conditioning
#' intervals \code{n_eff[i]} lower than the prescribed value \code{n[i]} and to a 
#' different number of data points between the groups.
#'   
#' @export
#' 
#' @references Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity analysis from a generic input-output sample, Env. Mod. & Soft.

pawn_split_sample <- function(X,Y, n=10, verbose = TRUE){
  
  ##############
  # Check inputs
  ##############
  
  stopifnot(is.matrix(X), is.numeric(X), is.numeric(Y))
  
  N <- nrow(X)
  M <- ncol(X)
  n2 <- nrow(Y)
  m <- ncol(Y)
  
  stopifnot(N == n2, m == 1)
  
  if(is.numeric(n)){
    if(length(n) == 1){
      n <- rep(n,M)
    } else if(length(n) != M){
      stop("n must have 1 or M components")
    }
  } else if (is.list(n) && length(n) != M){
    stop("If n is a list, it must have M components.")
  } else {
    stop("Wrong data type")
  }
  
  ######################
  #Create sub-samples:
  ######################
  
  n_eff <- numeric(M)
  NC <- vector(mode = "list", length = M)
  xc <- vector(mode = "list", length = M)
  Xk <- vector(mode = "list", length = M)
  XX <- list(list())
  YY <- list(list())
  
  for (ii in 1:M){
    pawn_samp <- split_sample(X[,ii],n[ii])
    idx <- pawn_samp$idx #"idx" contains the indices for each group for the i-th input
    Xk[[ii]] <- pawn_samp$Zk
    xc[[ii]] <- pawn_samp$Zc
    n_eff[ii] <- pawn_samp$n_eff
    
    
    XX[[ii]] <- rep(list(NA),n_eff[ii]) # conditioning samples for i-th input
    YY[[ii]] <- rep(list(NA),n_eff[ii])
    NC[[ii]] <- rep(NA,n_eff[ii])
    
    for (kk in 1:n_eff[ii]){ # loop over the conditioning intervals
      idxk <- idx == kk # indices of the k-th conditioning interval
      XX[[ii]][[kk]] <- X[idxk,]
      YY[[ii]][[kk]] <- Y[idxk]
      NC[[ii]][[kk]] <- sum(idxk)
    }
    
    if(n_eff[ii] < n[ii] && verbose){
      warning(paste0("For X", ii,", ", n_eff[ii], " groups were used instead of ",n[ii], 
                    ", so that values that are repeated several times belong to the same group"))
    }
    
  }
  
  robj <- (list(YY=YY, xc=xc, XX=XX, NC=NC, n_eff=n_eff, Xk=Xk))
  return(robj)
  
}



