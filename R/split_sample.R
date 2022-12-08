#' Split sample
#' 
#' Split a sample in n equiprobable groups based on the sample values (each groups contains approximately the same number of data points).
#' 
#' @param Z vector \code{(N, 1)}  sample of a model input  or output
#' @param ngroup number of groups considered (default: 10) 
#' 
#' @return List containing: 
#' \itemize{
#'   \item \code{idx} vector \code{N} respective group of the samples
#'   \item \code{Zk} vector \code{n_eff + 1} groups' edges (range of Z in each group)
#'   \item \code{Zc} vector \code{n_eff} groups' centers (mean value of Z in each group)
#'   \item \code{n_eff} scalar number of groups actually used to split the sample
#'}
#'
#' @export

split_sample <- function(Z, ngroup=10){
  
  ##############
  # Check inputs
  ##############
  
  stopifnot(is.numeric(Z))
  if(is.matrix(Z) && ncol(Z) >1) stop("Z should be a vector!") 
  Z <- as.vector(Z)
  
  ######################
  # Create sub-samples
  ######################
  
  N <- length(Z)
  
  n_eff <- ngroup
  
  Zu <- unique(Z) # distinct values of Z
  
  if (length(Zu) < ngroup ){ # if number of distinct values less than the specified number of groups
    n_eff <- length(Zu) # groups' centers are the different values of Xi
    Zc <- sort(Zu)
    Zk <- c(Zc, tail(Zc,n=1)) # groups' edges
    
  }  else {
    
    # Sort values of Z in ascending order:
    Z_sort <- sort(Z)
    # Define indices for splitting Z into ngroup equiprobable groups
    # (i.e. with the same number of values):
    
    split <- floor(seq(1,N,length.out = n_eff+1))
    split[2:n_eff] <- split[2:n_eff] + 1
    
    Zk <- as.vector(Z_sort[split])
    
    # Check that values that appear several times in Z belong to the same group:
    # To do this, we check that no values are repeated in Zk
    
    idx_keep <- rep(TRUE, n_eff+1) # index of value of Zk to be kept (TRUE) or to be dropped (FALSE)       
    
    for (kk in 1:length(Zk)){
      if (idx_keep[kk] == TRUE){ # if the value was not already discarded
        if (sum(Zk[(idx_keep)] == Zk[kk]) > 2){
          # the value Zk[k] appear more than twice, in this case, we 
          # remove one value from Zk (the k-th value)  
          idx_keep[kk] <- FALSE
        } else if (sum(Zk[idx_keep] == Zk[kk]) == 2){
          # the value Zk[k] appear exactly twice, in this case:
          if (kk < length(Zk)-3){
            # if the last and one before last values are the same, we
            # we do not change anything (this case will be specifically
            # treated when determining the respective groups of the 
            # sample)
            # Otherwise, we change the subsequent value of the edge 
            # Zk[kk+1] so that it is equal to the value in the sample
            # immediately higher than Zk[kk] and remove the next edge 
            # (Zk[kk+2]) to avoid having groups with very small sizes:
            idx_keep[kk+2] <- FALSE
            Zk[kk+1] <- Z_sort[which(Z_sort > Zk[kk])[1]]
          } else if(kk == length(Zk)-3){
            Zk[kk+1] <- Z_sort[which(Z_sort > Zk[kk])[1]]}
        }
      }
    }
    
    Zk <- Zk[idx_keep]
    n_eff <- length(Zk)-1
    
  }
  
  # Determine the respective groups of the sample:
  idx <- rep(-1,N)
  for (kk in 1:n_eff){
    if (kk < n_eff)
    {idx[which(Z >= Zk[kk] & Z < Zk[kk+1])] <- kk}
    else
    {idx[which(Z >= Zk[kk] & Z <= Zk[kk+1])] <- kk}
  }
  
  #Zc <- (Zk[-length(Zk)] + Zk[-1])/2
  Zc <- rep(0,n_eff)
  for (kk in 1:n_eff){
    Zc[kk] <- mean(Z[idx == kk]) 
    
  }
  
  # Check that all samples were assigned to a group
  if(any(idx == -1)) stop("Some samples were not assigned to any group") 
  
  robj <- (list(idx=idx, Zk=Zk, Zc=Zc, n_eff=n_eff))
  return(robj)
  
}