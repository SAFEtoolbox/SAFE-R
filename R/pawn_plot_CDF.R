#' Plot the CDFs estimated by the PAWN method 
#' 
#' This function plots the unconditional output Cumulative Distribution Funtions 
#' (i.e. when all inputs vary) and the conditional CDFs (when one input is fixed 
#' to a given conditioning interval, while the other inputs vary freely).
#' 
#' The unconditional and conditional CDF are computed by the function \link{pawn_CDF} \cr
#' See the help of \link{pawn_CDF} for further explanation.
#'
#' @param X matrix \code{(N, M)} set of inputs samples
#' @param Y matrix \code{(N, 1)} set of output samples 
#' @param n number of conditional intervals (default: 10) \cr
#' - integer if all inputs have the same number of groups \cr
#' - list of \code{M} integers otherwise
#' @param n_col integer, number of panels per row in the plot (default: \code{min(5, M))}
#' @param y_label string, legend for the horizontal axix (default: 'output y')
#' @param labelinput vector, label the vertical axis (default: \code{["X1","X2",...,"XM"]})
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
#' \code{xc[[i]]} vector of length \code{n_eff[i]} contains the centers of the \code{n_eff[i]} conditioning intervals for the i-th input.
#' }
#'   
#' @export

pawn_plot_CDF <- function(X, Y, n=10, n_col=5, y_label='output y', labelinput='', verbose = TRUE){
  
  stopifnot(is.matrix(X), is.numeric(X), is.numeric(Y))
  
  ###########################################################################
  # Check optional inputs for plotting
  ###########################################################################
  
  stopifnot(is.scalar(n_col), n_col > 0)
  
  if(length(labelinput) == 1 && labelinput == ""){
    for (ii in 1:M){
      labelinput[ii] <-paste0('X',ii) 
    }
  }

  
  ###########################################################################
  # Split the input sample
  ###########################################################################
  
  pawn_samp <- pawn_split_sample(X, Y, n, verbose) 
  XX <- pawn_samp$XX
  YY <- pawn_samp$YY
  xc <- pawn_samp$xc
  NC <- pawn_samp$NC
  Xk <- pawn_samp$Xk
  n_eff <- pawn_samp$n_eff
  
  N <- nrow(X)
  M <- ncol(X)
  
  
  ###########################################################################
  # Compute CDFs
  ###########################################################################
  
  pawnCDF <- pawn_CDF(X, Y, n, dummy=FALSE, verbose = 1)
  YF <- pawnCDF$YF
  FU <- pawnCDF$FU
  FC <- pawnCDF$FC
  
  
  ###########################################################################
  # Plot
  ###########################################################################
  
  plot_list = list()
  
  for(ii in 1:M){
    
    dat_cond <- data.frame(x = rep(YF,n_eff[ii]),
                           y = unlist(FC[[ii]]),
                           group = rep(1:n_eff[ii],each=length(YF)),
                           clab = as.factor(round(rep(unlist(xc[ii]),each=length(YF)),1))
    ) 
    
    # check if the rounding kept the correct number of conditional CDFs
    zz = 1
    while(nlevels(dat_cond$clab) < nlevels(as.factor(dat_cond$group))){
      
      dat_cond$clab = as.factor(round(rep(unlist(xc[ii]),each=length(YF)),1+zz))
      zz = zz + 1
    }
    
    dat_uncond <- data.frame(x = rep(YF,M),
                             y = unlist(FU[[ii]]),
                             group = rep(1:M, each = length(YF))
    )
    
    pp <- ggplot(data = dat_cond) + geom_line(mapping = aes(x=x, y=y, colour = clab), size = 1) +
      scale_colour_grey(name = "") + xlab(y_label) + labs(title=labelinput[ii]) + ylab("") +
      geom_line(data = dat_uncond, inherit.aes = FALSE, aes(x = x, y = y), size = 1, color = "red") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
      ylim(0,1) + theme(aspect.ratio = 1) +
      theme_bw() 
    
    
    plot_list[[ii]] <- pp
  }
  
  n_col <- min(n_col, M)
  n_row <- ceiling(M/n_col)
  do.call("grid.arrange", c(plot_list, ncol=n_col))
  
  robj <- list(YF=YF,FU=FU,FC=FC,xc=xc)
  return(robj)
}


