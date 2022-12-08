#' This function can be used as input argument ("output_condition") when
#' applying pawn_indices, pawn_convergence, pawn_plot_ks
#' 
#' Returns all indices of \code{y}
#' This is often the default function
#' 
#' @export

allrange <- function(y,par){
  
  stopifnot(is.numeric(y), is.scalar(par) || is.list(par))
  if(is.matrix(y) && ncol(y) >1) stop("y should be a vector!") 
  y <- as.vector(y)
  
  idx <- rep(TRUE,length(y))
  return(idx)
}
