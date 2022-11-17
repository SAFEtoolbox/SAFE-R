#' Scatter plots of \code{y} against \code{X}
#'
#' This function produces scatter plots of the model ouput \code{y} against model inputs \eqn{x(1), x(2), ..., x(M)}
#'
#' @param X matrix \code{(N, M)} of \code{N} inputs samples
#' @param Y vector \code{N} of associated ouput samples
#' @param prnam vector of characters containing the parameters names
#' @param ngr number of groups into which the parameters are divided
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
#' # Step 2 (sampling and model evaluation)
#' # ############################
#' N <- 3000
#' X <- AAT_sampling("lhs", M, distr_fun, distrpar, N)
#' Y <- model_evaluation(fun_test, X)         
#' # ############################
#' ## Step 3 (Scatter plots)
#' # ############################
#' # Use scatter plots of inputs againts output to visually assess 
#' # direct effects:
#' scatter_plots(X,Y)
#'
scatter_plots <- function(X, Y, prnam = NULL, ngr = 0){
  
  if( is.null(prnam) ){
    prnam <- as.character( 1:ncol(X) )
  }
  
  if( ngr <= 1 ){
    ns <- nrow( X )
    npr <- ncol( X )
    dat <- data.frame(x = as.vector(X), 
                      y = rep(Y, npr), 
                      parnam = factor(rep(prnam, each = ns), levels = prnam))
    
    return(ggplot(data = dat, 
                  mapping = aes(x = x, y = y, parnam = parnam)) + 
             facet_grid(. ~ parnam, scales = "free") + 
             geom_point() + ylab("Output") + xlab("Inputs") + 
             theme_bw()) 
  } else {
    
    N <- length( Y )
    N <- floor(N / ngr) * ngr
    Y <- Y[ 1:N ]
    X <- X[1:N, ]
    np <- length( prnam )
    
    # First part
    ord <- order( Y )
    
    split <- seq(0, N, by = floor(N / ngr))
    
    idx <- numeric(N)
    X_ <- vector("list", ngr)
    Y_ <- matrix(nrow = floor(N / ngr), ncol = ngr)
    
    for (i in 1:ngr){
      idx[ord[(split[i] + 1):split[i + 1]]] <- i
      X_[[i]] <- matrix(NA, nrow = floor(N / ngr), ncol = np)
      X_[[i]] <- X[idx == i,]
      Y_[,i] <- Y[idx == i]
    }
    
    Yk <- apply(Y_, 2, max)
    
    # Second part
    ng <- length( X_ )
    ns <- nrow( X_[[1]] )
    dat <- list()
    for (ii in 1:np){
      dat[[ii]] <- data.frame(x = do.call("c", 
                                          lapply(1:ngr, 
                                                 function(.kk) X_[[.kk]][ , ii])), 
                              y = as.vector(Y_))
      dat[[ii]]$group <- as.factor( round(rep(Yk, each = ns), 2) )
      dat[[ii]]$parnam <- prnam[ii]
    }
    dat <- do.call("rbind", dat)
    dat$parnam <- factor( dat$parnam, levels = prnam )
    
    .pl <- ggplot(data = dat, 
                  mapping = aes(x = x, y = y, col = group, parnam = parnam)) + 
      facet_grid(.~parnam, scales = "free") + 
      geom_point() + ylab("Output") + xlab("Inputs") + 
      theme_bw() + labs(color='\t Upper \n bound \n output')
    
    return( .pl ) 
    
  }
}