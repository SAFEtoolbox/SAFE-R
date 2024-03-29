#' PAWN sensitivity indices
#' 
#' Compute the PAWN sensitivity indices. The method was first introduced
#' in Pianosi and Wagener (2015). Here indices are computed following the
#' approximation strategy proposed by Pianosi and Wagener (2018), which can be
#' applied to a generic input/output sample. \cr
#' The function splits the generic output sample to create the conditional
#' output by calling internally \link{pawn_split_sample}. The
#' splitting strategy is an extension of the strategy for uniformy distributed
#' inputs described in Pianosi and Wagener (2018) to handle inputs sampled
#' from any distribution(see help of \link{pawn_split_sample} for further explanation). \cr
#' Indices are then computed in two steps: \cr
#' 1. compute the Kolmogorov-Smirnov (KS) statistic between the empirical
#' unconditional CDF and the conditional CDFs for different conditioning intervals \cr
#' 2. take a statistic (median, mean and max) of the results.
#' 
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, 1)} set of output samples 
#' @param n number of conditional intervals to assess the conditional CDFs (default: 10)
#' - integer if all inputs have the same number of groups
#' - list of \code{M} integers otherwise
#' @param Nboot scalar, number of bootstrap resamples to derive confidence intervals
#' @param dummy boolean, if \code{dummy} is TRUE, an artificial input is added to the set of inputs 
#' and the sensitivity indices are calculated for the dummy input.
#' The sensitivity indices for the dummy input are estimates of the approximation
#' error of the sensitivity indices and they can be used for screening, i.e. to 
#' separate influential and non-influential inputs as described in Khorashadi Zadeh et al. (2017)
#' (default: False) \cr
#' (see (*) for further explanation). \cr \cr
#' ADVANCED USAGE: \cr
#' for Regional-Response Global Sensitivity Analysis:
#' @param output_condition function, condition on the output value to be  used to calculate KS. \cr
#'  Use the function: \cr
#'  - \link{allrange} to keep all output values \cr
#'  - \link{below} to consider only output values below a threshold value (\code{Y <= Ythreshold}) \cr
#'  - \link{above} to consider only output values above a threshold value (\code{Y >= Ythreshold}) \cr
#' @param par list, specify the input arguments of the \code{output_condition} function, i.e. the threshold
#'  value when \code{output_condition} is 'above' or 'below'.
#'   For more sophisticate conditions, the user can define its own function \code{output_condition} with the following structure:
#'   \code{idx <- output_condition(Y, param)}
#'   where: \code{Y} output samples \code{(N, 1)}
#'   \code{param} parameters to define the condition (list of any size)
#'   \code{idx} logical values, True if condition is satisfied, False otherwise \code{(N,1)}
#' 
#' @return List containing: 
#' \itemize{
#'   \item \code{KS_median} array \code{(Nboot, M)} median KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_mean} array \code{(Nboot, M)} mean KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_max} array \code{(Nboot, M)} max KS across the conditioning intervals 
#'   (one value for each input and each bootstrap resample)
#'   \item \code{KS_dummy} array \code{(Nboot, 1)} KS of dummy input (one value for each bootstrap resample)
#'   } 
#' NOTES: \cr
#' (*) For screening influential and non-influential inputs, we recommend the use of 
#' the maximum KS across the conditioning intervals (i.e. output argument \code{KS_max}), 
#' and to compare \code{KS_max} with the index of the dummy input as in 
#' Khorashadi Zadeh et al. (2017). \cr 
#' (**) For each input, the number of conditioning intervals which is actually used (\code{n_eff[i]}) 
#' may be lower than the prescribed number of conditioning intervals (\code{n[i]}) to ensure that 
#' input values that are repeated several time belong to the same group.
#' See the help of \link{pawn_split_sample} for further details.
#' @export
#' 
#' @references Pianosi, F. and Wagener, T. (2018), Distribution-based sensitivity analysis from a generic input-output sample, Env. Mod. & Soft., 108, 197-207. \cr
#' Pianosi, F. and Wagener, T. (2015), A simple and efficient method for global sensitivity analysis based on cumulative distribution functions, Env. Mod. & Soft., 67, 1-11. \cr
#' Khorashadi Zadeh et al. (2017), Comparison of variance-based and moment-independent global sensitivity analysis approaches by application to the SWAT model, Environmental Modelling & Software,91, 210-222.
#' 
#' @examples
#' library(gridExtra)
#' library(ggplot2)
#' # Create a generic input-output sample:
#' N <- 5000 # number of samples
#' M <- 3 # number of inputs
#' distrpar <-  c(-pi, pi)
#' X <- AAT_sampling('lhs', M, 'unif', distrpar, N)
#' Y <- model_execution("ishigami_homma_function", X)
#' x_labels = c('x(1)', 'x(2)', 'x(3)')
#' 
#' # Compute PAWN sensitivity indices:
#' n <- 10 # number of conditioning intervals
#' pawn_ind <- pawn_indices(X, Y, n)
#' KS_median <- pawn_ind$KS_median
#' KS_mean <- pawn_ind$KS_mean
#' KS_max <- pawn_ind$KS_max
#' 
#' # Plot results:
#' dev.new()
#' p1 <- boxplot1(as.vector(KS_median), prnam = x_labels) + ylab("KS (median)") 
#' p2 <- boxplot1(as.vector(KS_mean), prnam = x_labels) + ylab("KS (mean)") 
#' p3 <- boxplot1(as.vector(KS_max), prnam = x_labels) + ylab("KS (max)") 
#' grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
#' 
#' # Compute sensitivity indices for dummy input as well:
#' pawn_ind <- pawn_indices(X, Y, n, dummy = TRUE)
#' KS_median <- pawn_ind$KS_median
#' KS_dummy <- pawn_ind$KS_dummy
#' KS_median_d <- c(KS_median,KS_dummy) 
#' dev.new()
#' boxplot1_dummy(mu = as.vector(KS_median_d), prnam = x_labels) + ylab("KS (median)") 


pawn_indices<- function(X, Y, n=10, Nboot = 1, dummy = FALSE, output_condition = allrange, par = list()){
  
  N <- nrow(X)
  M <- ncol(X)
  
  stopifnot(is.matrix(X), is.numeric(X), is.numeric(Y), is.scalar(Nboot), Nboot >=1, Nboot == floor(Nboot), is.logical(dummy))
  stopifnot(is.function(output_condition),is.scalar(par) || is.list(par))
  
  #### Compute indices
  
  # Initialize sensitivity indices
  KS_median <- matrix(NA,Nboot,M)
  KS_mean <- matrix(NA,Nboot,M)
  KS_max <- matrix(NA,Nboot,M)
  KS_dummy <- NA # In case dummy = FALSE
  if(dummy == TRUE){
    KS_dummy <- matrix(NA,Nboot,1)
  }
  
  # Compute conditional CDFs
  # (bootstrapping is not used to assess conditional CDFs):
  
  for(bb in 1:Nboot){
    
    # Compute empirical unconditional and conditional CDFs
    pawnCDF <- pawn_CDF(X, Y, n, dummy, verbose = bb == 1)
    YF <- pawnCDF$YF # Set points at which the CDFs will be evaluated
    FU <- pawnCDF$FU # Unconditional CDFs
    FC <- pawnCDF$FC # Conditional CDFs
    idx_bootsize_min <- pawnCDF$idx_bootsize_min
    FC_dummy <- pawnCDF$FC_dummy
    
    # Compute KS statistic between conditional and unconditional CDFs:
    KS_all <- pawn_ks(YF, FU, FC, output_condition, par)
    
    # KS_all is a list (M elements) and contains the value of the KS for
    # for each input and each conditioning interval. KS[i] contains values
    # for the i-th input and the n_eff[i] conditioning intervals, and it
    # is a matrix of shape (n_eff[i],1).
    
    #  Take a statistic of KS across the conditioning intervals:
    KS_median[bb,] <- sapply(KS_all,median)
    KS_mean[bb,] <- sapply(KS_all,mean)
    KS_max[bb,] <- sapply(KS_all,max)
    
    if(dummy == TRUE){
      # Compute KS statistic for dummy parameter:
      KS_dummy[bb] <- pawn_ks(YF,FU[idx_bootsize_min],list(list(FC_dummy)),output_condition,par)
    }
    
  }
  
  KS_dummy <- unlist(KS_dummy)
  
  robj <- (list(KS_median=KS_median, KS_mean=KS_mean, KS_max=KS_max, KS_dummy=KS_dummy))
  return(robj)
  
  
}