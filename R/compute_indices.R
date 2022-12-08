compute_indices <- function(X, Y, threshold, flag){
	stopifnot(flag %in% 1:3)
	
	N <- nrow(X)
	M <- ncol(X)
	P <- ncol(Y)
	
	idxb <- colSums(t(Y) < threshold) == P
	stat <- numeric(M)

## Define below and above subsamples:
	Xb  <- X[idxb,] 
	Xb <- matrix(Xb, nrow = sum(idxb))
	Nb  <- length(Xb) ## number of behavioural parameterizations
	Xnb <- X[!idxb,]
	Nnb <- length(Xnb) # number of non-behavioural parameterizations

if(Nb <= 0){ 
	stop('Cannot find any output value below the threshold! Try increasing the threshold value')
} else {

	if(Nnb <= 0)
    stop('Cannot find any output value above the threshold! Try reducing the threshold value')

}

# Perform RSA
   
   if(flag < 3){
   
    for (i in 1:M){
        # Approximate behavioural and non-behavioural CDFs:
      
      aCDF <- approximatecdfpair(Xb[,i], Xnb[,i])
      
        # Compute distance between CDFs:
        if (flag == 1){
            stat[i] <- max(abs(aCDF$CDF1 - aCDF$CDF2))
            }
        if (flag == 2){

        		stat[i] <- trapz(aCDF$x, pmax(aCDF$CDF1, aCDF$CDF2)) - trapz(aCDF$x, pmin(aCDF$CDF1, aCDF$CDF2))
           
        }
    }
    
    } else { # if (flag == 3)
  
        # Ranges of input factors that produce ''behavioural'' output:
        xmin <- apply(Xb, 2, min)
        xmax <- apply(Xnb, 2, max)
        
        # Compute the relative reduction wrt the original ranges of input factors
        # (as they appear in the sample ''X''):
        
        stat <- (xmax - xmin) / apply(X, 2, function(h) diff(range(h)))
    }

	robj <- list(stat = stat, idxb = idxb)

return(robj)
}