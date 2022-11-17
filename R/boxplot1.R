#' Boxplot of vectors

#' Boxplots when the mean, lower values and upper values are specified.

#' @param mu vector (\code{M}) of mean or median values to be plotted
#' @param lb vector (\code{M}) of lower values to be plotted. 
#' @param au vector (\code{M}) of upper values to be plotted. 
#' @param prnam labels for the x-axis of the boxplot
#' @seealso \code{\link{boxplot2}} \code{\link{boxplot}} \code{\link{plot}}
#' @export
#' @examples
#' # Setup the model and define input ranges
#' myfun  <- "ishigami_homma_function"
#' M <- 3
#' DistrFun <- "unif"
#' DistrPar <-  c(-pi, pi)
#' # Sample parameter space using the resampling strategy proposed by 
#' # (Saltelli, 2008; for reference and more details, see help of functions
#' # vbsa_resampling and vbsa_indices) 
#' SampStrategy <- "lhs"
#' N <- 3000 # Base sample size.
#' # Comment: the base sample size N is not the actual number of input 
#' # samples that will be evaluated. In fact, because of the resampling
#' # strategy, the total number of model evaluations to compute the two
#' # variance-based indices is equal to N*(M+2) 
#' X <- AAT_sampling(SampStrategy, M, DistrFun, DistrPar, 2 * N)
#' XABC <- vbsa_resampling(X)
#' # Run the model and compute selected model output at sampled parameter
#' # sets:
#' YA <- model_evaluation(myfun, XABC$XA) # size (N,1)
#' YB <- model_evaluation(myfun, XABC$XB) # size (N,1)
#' YC <- model_evaluation(myfun, XABC$XC) # size (N*M,1)
#' # Compute main (first-order) and total effects:
#' ind <- vbsa_indices(YA, YB, YC)
#' Si <- ind[1,]
#' STi <- ind[2,] 
#' # Plot results:
#' # plot main and total separately
#' par(mfrow = c(1, 2))
#' boxplot1(Si) + ggtitle("Si")
#' boxplot1(STi) + ggtitle("STi")
#'

boxplot1 <- function(mu, lb = NULL, ub = NULL, prnam = NULL){
  
  dat <- data.frame(x = factor( prnam, levels = prnam ),
                    mu = mu)
  
  .pl <- ggplot(data = dat, mapping = aes(x = x, y = mu))
  
  if( !is.null(lb) && !is.null(ub) ){
    dat$lb = lb 
    dat$ub = ub
    .pl <- .pl + geom_errorbar(mapping = aes(ymin = lb, ymax = ub), width = 0.5) 
  }
  
  .pl <- .pl + geom_point(color = 'red', size = 3) + theme_bw() + xlab(NULL) + ylab("mvd") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))
  
  return( .pl )
  
}

