#' This function can be used as input argument ("output_condition") when
#' applying pawn_indices, pawn_convergence, pawn_plot_ks
#' 
#' Returns the indices of \code{y} below the scalar \code{par} threshold
#' 
#' @export

below <- function(y,par){
  
  stopifnot(is.scalar(par) || is.list(par), is.numeric(y))
  if(is.matrix(y) && ncol(y) >1) stop("y should be a vector!") 
  y <- as.vector(y)
  
  idx <- y <= par
  return(idx)
}
