#' Compute and plot the Kolmogorov-Smirnov (KS) statistic
#' 
#' This function computes and plots the KS statistic
#' between conditional and unconditional output CDFs for each input and each
#' conditioning interval.
#' 
#' The unconditional and conditional CDF are computed by the function \link{pawn_CDF} \cr
#' See the help of \link{pawn_CDF} for further explanation.
#' 
#' @param YF vector \code{(P, 1)}, values of Y at which the Cumulative  Distribution Functions (CDFs) FU and FC are given
#' @param FU list of \code{M} elements, values of the empirical unconditional output CDFs.
#' @param FC list of lists, values of the empirical conditional output CDFs for each input and each conditioning interval.
#' \code{FC[[i]]} is a list of \code{n_eff[i]} CDFs conditional to the i-th input.
#' \code{FC[[i]][k]} is obtained by fixing the i-th input to its k-th conditioning interval (while the other inputs vary freely)
#' @param xc list of \code{M} elements, subsamples' centers (i.e. mean value of the \code{Xi} over each conditioning interval)
#' \code{xc[[i]]} vector of length \code{n_eff[i]} contains the centers of the \code{n_eff[i]} conditioning intervals for the i-th input.
#' @param n_col integer, number of panels per row in the plot (default: \code{min(5, M))}
#' @param x_labels string, legend for the x-axis (input name), (default: \code{["X1","X2",...,"XM"]}) \cr
#' ADVANCED USAGE: \cr
#' for Regional-Response Global Sensitivity Analysis: \cr
#' @param output_condition function (see the help for \link{pawn_indices} for information on this optional input)
#' @param par list (see the help for \link{pawn_indices} for information on this optional input) \cr
#' 
#' @details
#' NOTE: \cr
#' \code{YF}, \code{FU}, \code{FC} and \code{xc} are computed using the function \link{pawn_plot_CDF}.
#' 
#' @return \code{KS_all} list of \code{M} elements. KS-statistic calculated between conditional and unconditional output for the M inputs 
#' and the \code{n_eff} conditioning intervals. \code{KS_all[i]} contains the KS values for the i-th input and the 
#' \code{n_eff[i]} conditioning intervals.
#' 
#' @export
#' 
#' @examples
#'
#' # Create a generic input-output sample:
#' N <- 1000 # number of samples
#' M <- 3 # number of inputs
#' distrpar <-  c(-pi, pi)
#' X <- AAT_sampling('lhs', M, 'unif', distrpar, N)
#' Y <- model_execution("ishigami_homma_function", X)
#' x_labels = c('x(1)', 'x(2)', 'x(3)')
#' 
#' # Calculate CDFs:
#' n <- 10 # number of conditioning intervals
#' pawn_cdf <- pawn_CDF(X, Y, n)
#' YF <- pawn_cdf$YF
#' FU <- pawn_cdf$FU
#' FC <- pawn_cdf$FC
#' xc <- pawn_cdf$xc
#' 
#' # Calculate and plot KS statistics:
#' KS_all <- pawn_plot_ks(YF, FU, FC, xc)

pawn_plot_ks <- function(YF, FU, FC, xc, n_col=5, x_labels=NULL, output_condition = allrange, par = list()){
  
  ###########################################################################
  # Check optional inputs for plotting
  ###########################################################################
  
  stopifnot(is.scalar(n_col), n_col > 0)
  
  
  ###########################################################################
  # Calculate KS-statistic
  ###########################################################################
  
  KS_all <- pawn_ks(YF, FU, FC, output_condition = allrange, par = list()) 
  
  M <- length(KS_all)
  
  if(is.null(x_labels)){
    for (ii in 1:M){
      x_labels[ii] <-paste0('X',ii) 
    }
  }
  
  
  ###########################################################################
  # Plot
  ###########################################################################
  
  plot_list = list()
  
  for(ii in 1:M){
    
    dat_ks <- data.frame(x = unlist(xc[ii]),
                         y = KS_all[[ii]],
                         group = as.factor(rep(1:length(KS_all[[ii]])))
    ) 
    
    
    
    pp <- ggplot(data = dat_ks) + geom_line(mapping = aes(x=x, y=y)) +
      geom_point(mapping = aes(x=x, y=y, colour = group), size = 2) +
      scale_colour_grey(guide=FALSE) + xlab(x_labels[ii]) + ylab("KS") + ylim(0,1) +
      theme_bw() 
    
    
    plot_list[[ii]] <- pp
  }
  
  n_col <- min(n_col, M)
  n_row <- ceiling(M/n_col)
  do.call("grid.arrange", c(plot_list, ncol=n_col))
  
  return(KS_all)
  
}


