#' Mean and standard deviation of Elementary Effects
#'
#' Compute mean and standard deviation of Elementary Effects (EEs) using an increasing number of samples.
#'
#' @param EE matrix (\code{r}, \code{M}) of \code{r} elementary effects
#' @param rr vector (\code{R}), sample sizes at which indices will be estimated (must be a vector of integer positive values not exceeding \code{r})
#' @param Nboot scalar, number of resamples used for boostrapping (default: 0) 
#' @param alfa scalar, significance level for the confidence intervals estimated by bootstrapping (default: 0.05)
#'
#' @return List containing: 
#' \itemize{
#'   \item \code{m_r} mean of EEs at different sampling size
#'   \item \code{s_r} standard deviation of EEs at different sampling size
#' }
#' If \code{Nboot > 1} it also contains
#' \itemize{
#'   \item \code{m_lb_r} lower bound of EEs mean at different sampling size
#'   \item \code{m_ub_r} upper bound of EEs mean at different sampling size
#'   \item \code{s_lb_r} lower bound of EEs std.dev. at different sampling size
#'   \item \code{s_ub_r} upper bound of EEs std.dev. at different sampling size
#'   \item \code{m_sd_r} standard deviation of EEs mean at different sampling size
#'   \item \code{s_sd_r} standard deviation of EEs std.dev. at different sampling size. 
#' }
#' All output arguments are matrices of size (\code{R},\code{M})
#' @seealso \code{\link{EET_indices}} \code{\link{EET_plot}}

#' @export

#' @examples

#' # See the demo
#' # demo("workflow_eet_hymod")
#' # or
#' # demo("workflow_eet_hbv")



EET_convergence <- function(EE, rr, Nboot = 0, alfa = 0.05){


#########
# Check inputs
#########

r <- rr[length(rr)]

R <- length(rr)
M <- ncol(EE)

m_r <- matrix(nrow = R, ncol = M)
s_r  <- matrix(nrow = R, ncol = M)

stopifnot(is.matrix(EE), is.numeric(EE), floor(rr) == rr, rr > 0, diff(rr) >0, r <= nrow(EE),
is.scalar(Nboot), Nboot >=0, Nboot == floor(Nboot),
is.numeric(alfa), alfa <= 1, alfa >= 0)

############
# Compute indices
############

if (Nboot > 1){
	
    m_sd_r <- matrix(nrow = R, ncol = M)
    m_lb_r <- matrix(nrow = R, ncol = M)
    m_ub_r <- matrix(nrow = R, ncol = M)
    s_sd_r <- matrix(nrow = R, ncol = M)
    s_lb_r <- matrix(nrow = R, ncol = M)
    s_ub_r <- matrix(nrow = R, ncol = M)
    
    for (j in 1:R){
    	
        idx <- sample.int(r, rr[j])
        EEi <- EE[idx,]
        
        bootsize <- rr[j]
        
        B <- matrix(sample.int(bootsize, bootsize * Nboot, replace = TRUE), bootsize, Nboot)
        
        	m_j <- apply(B, 2, function(b) colMeans(EEi[b,]))
        	s_j <- apply(B, 2, function(b) apply(EEi[b,], 2, sd))

        	 m_r[j,] <- rowMeans(m_j)
        	 m_sd_r[j,] <- apply(m_j, 1, sd)
        	 m_sorted <- apply(m_j, 1, sort) 
        	 m_lb_r[j,] <- m_sorted[max(1,round(Nboot * alfa / 2)),]
        	 m_ub_r[j,] <- m_sorted[round(Nboot * (1 - alfa / 2)),]
        	 
        	 s_r[j,] <- rowMeans(s_j)
        	 s_sd_r[j,] <- apply(s_j, 1, sd)
        	 s_sorted <- apply(s_j, 1, sort) 
        	 s_lb_r[j,] <- s_sorted[max(1,round(Nboot * alfa / 2)),]
        	 s_ub_r[j,] <- s_sorted[round(Nboot *(1 - alfa / 2)),]
   
    }
    
     robj <- list(m_r = m_r, s_r = s_r, m_sd_r = m_sd_r, s_sd_r = s_sd_r, m_lb_r = m_lb_r, s_lb_r = s_lb_r, m_ub_r = m_ub_r, s_ub_r = s_ub_r)
    
} else {
    
    
    for (j in 1:R){
        
        idx <- sample.int(r, rr[j])
        EEi <- EE[idx,]
        
        m_r[j,] <- colMeans(EEi)
        s_r[j,] <- apply(EEi, 2, sd)
    
    }
    
  robj <- list(m_r = m_r, s_r = s_r)
  
}

return(robj)

}