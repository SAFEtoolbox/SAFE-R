#' Boxplot of vectors

#' Boxplots when the mean, lower values and upper values are specified.

#' @param mu vector (\code{M}) of mean or median values to be plotted
#' @param lb vector (\code{M}) of lower values to be plotted. 
#' @param au vector (\code{M}) of upper values to be plotted. 
#' @param prnam labels for the x-axis of the boxplot
#' @seealso \code{\link{boxplot2}} \code{\link{boxplot}} \code{\link{plot}}
#' 
#' @export
#' 
#' @examples
#' library(ggplot2)
#' # Setup the model and define input ranges
#' myfun  <- "ishigami_homma_function"
#' M <- 3
#' DistrFun <- "unif"
#' DistrPar <-  c(-pi, pi)
#' SampStrategy <- "lhs"
#' N <- 3000 # Base sample size.
#' X <- AAT_sampling(SampStrategy, M, DistrFun, DistrPar, N)
#' x_labels <- c('x(1)', 'x(2)', 'x(3)')
#' 
#' # Run the model 
#' Y <- model_execution(myfun, X) 
#' 
#' # Compute indices with the PAWN method
#' n <- 10
#' pawn_ind <- pawn_indices(X, Y, n, dummy = TRUE)
#' KS_max <- pawn_ind$KS_max 
#' KS_dummy <- pawn_ind$KS_dummy 
#' KS_max_d <- c(KS_max,KS_dummy) 
#' 
#' # Plot results:
#' boxplot1_dummy(mu = KS_max_d, prnam = x_labels) + ylab("KS")

boxplot1_dummy <- function(mu, lb = NULL, ub = NULL, prnam = NULL){
  
  mu1 <- mu[1:(length(mu)-1)]
  lb1 <- lb[1:(length(lb)-1)]
  ub1 <- ub[1:(length(ub)-1)]
  mu2 <- tail(mu,n=1) 
  lb2 <- tail(lb,n=1)
  ub2 <- tail(ub,n=1)
  
  dat <- data.frame(x = factor( prnam, levels = prnam ),
                    mu = mu1)
  
  .pl <- ggplot(data = dat, mapping = aes(x = x, y = mu1))
  
  if( !is.null(lb1) && !is.null(ub1) ){
    dat$lb = lb1 
    dat$ub = ub1
    .pl <- .pl + geom_errorbar(mapping = aes(ymin = lb1, ymax = ub1), width = 0.5) 
  }
  
  .pl <- .pl + geom_point(color = 'red', size = 3) + theme_bw() + xlab(NULL) + ylab("mvd") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    geom_hline(yintercept=mu2, size = 0.25, color = "grey49") +
    geom_hline(yintercept=lb2, size = 0.25, linetype = "dashed", color = "grey49") +
    geom_hline(yintercept=ub2, size = 0.25, linetype = "dashed", color = "grey49") +
    annotate("text", x = length(dat$x)/2, y = 0, label = "Threshold for non-influential input factors",
             size = 3, color = "grey49")
  
  return( .pl )
  
}

