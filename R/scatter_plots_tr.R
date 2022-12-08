#' Scatter plots of \code{y} against \code{X}
#'
#' This function produces scatter plots of the model ouput \code{y} against model inputs \eqn{x(1), x(2), ..., x(M)}
#'
#' @param X matrix \code{(N, M)} of \code{N} inputs samples
#' @param Y vector \code{N} of associated ouput samples
#' @param prnam vector of characters containing the parameters names
#' @param idx indices of datapoints to be highlighted
#' @export
#' @examples
#' #############################
#' # Step 1 (setup the model)
#' ############################
#' fun_test  <- "ishigami_homma_function"
#' M <- 3
#' distr_fun <- "unif"
#' distrpar <-  c(-pi, pi)
#' # ############################
#' # Step 2 (sampling and model execution)
#' # ############################
#' N <- 3000
#' X <- AAT_sampling("lhs", M, distr_fun, distrpar, N)
#' Y <- model_execution(fun_test, X)         
#' # ############################
#' ## Step 3 (Scatter plots)
#' # ############################
#' # Use scatter plots of inputs againts output to visually assess 
#' # direct effects:
#' scatter_plots(X,Y)
#'
scatter_plots_tr <- function(X, Y, prnam = NULL, idx = NULL){
  
  if( is.null(prnam) ){
    prnam <- as.character( 1:ncol(X) )
  }
  
  if( !is.null(idx)){
    ns <- nrow( X )
    npr <- ncol( X )
    dat <- data.frame(x = as.vector(X), 
                      y = rep(Y, npr), 
                      parnam = factor(rep(prnam, each = ns), levels = prnam))
    
    Xb <- X[idx,]
    Yb <- Y[idx]
    nsb <- nrow( Xb )
    nprb <- ncol( Xb )
    
    dat_tr <- data.frame(x = as.vector(Xb), 
                         y = rep(Yb, nprb), 
                         parnam = factor(rep(prnam, each = nsb), levels = prnam))
    
    

    .pl <- ggplot(data = dat, 
                  mapping = aes(x = x, y = y, parnam = parnam)) + 
      facet_grid(. ~ parnam, scales = "free") + 
      geom_point(data=dat, fill="blue", colour="blue") + ylab("Output") + xlab("Inputs") + 
      theme_bw()
    
    .pl <- .pl + geom_point(data=dat_tr, fill="red", colour="red")

    return( .pl ) 
  }
}