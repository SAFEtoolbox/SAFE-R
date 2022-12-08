compute_indices_rsa_gr <- function(X, Y, ngroup, Nboot=0){
  
  # This function computes the sensitivity indices for Regional Sensitivity Analysis 
  # with grouping for ONE sample/bootstrap resample.
  
  N <- nrow(X)
  M <- ncol(X)
  
  # Split the output sample in equiprobable groups:
  split_sam <- split_sample(Y, ngroup)
  idx <- split_sam$idx
  Yk <- split_sam$Zk
  ngroup_eff <- split_sam$n_eff
  
  if(ngroup_eff < ngroup){
    warning(paste0(ngroup_eff," groups were used instead of ", ngroup, " so that values that are repeated several times belong to the same group"))
  }
  
  # Initialize arrays of indices
  mvd <- matrix(NA,ngroup_eff*(ngroup_eff-1)/2,M)
  spread <- matrix(NA,ngroup_eff*(ngroup_eff-1)/2,M)
  
  for(ii in 1:M){
    # Approximate CDF of the i-th parameter for each group:
    L <- length(unique(X[,ii]))
    CDF_ <- matrix(NA,L,ngroup_eff)
    xx <- unique(sort(X[,ii]))
    
    for(jj in 1:ngroup_eff){
      CDF_[,jj] <- empiricalcdf(X[idx == jj, ii], xx)
    }
    
    # Compute the distance between the different CDFs
    count = 1
    for(jj in 1:(ngroup_eff-1)){
      for(kk in seq(jj+1,ngroup_eff,1)){
        mvd[count,ii] <- max(abs(CDF_[,jj] - CDF_[,kk]))
        spread[count,ii] <- trapz(xx, pmax(CDF_[,jj], CDF_[,kk])) - trapz(xx, pmin(CDF_[,jj], CDF_[,kk]))
        count <- count + 1
      }
    }
  }
  
  if(Nboot > 0){
    robj <- (list(mvd=mvd, spread=spread))
    }
  else {
    robj <- (list(mvd=mvd, spread=spread, idx = idx, Yk = Yk))
  }
  
  return(robj)
}