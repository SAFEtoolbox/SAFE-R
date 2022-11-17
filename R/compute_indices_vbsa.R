compute_indices_vbsa <- function(YA, YB, YC){
#
#    YA - vector (N,1) #
#    YB - vector (N,1) #
#    YC - matrix (N,M) #

if(!is.matrix(YC)){
	YC <- matrix(YC, N, M) 
	} 

if(ncol(YC) == 1) YC <- matrix(YC, nrow = nrow(YA))

N <- nrow(YC)
M <- ncol(YC)

	f0  <- mean(YA)
	VARy <- mean(YA^2) - f0^2
		
	sumYaYC <- c(t(YA) %*% YC)
	
	Si <- (1 / N * sumYaYC - f0^2) / VARy
	
	# This is Eq (4.21) in Saltelli et al. (2008)
    # and also Eq. (12) in Saltelli et al. (2010), where the method is
    # attributed to Sobol (1993).
	
	STi <- 1 - (1 / N * c(t(YB) %*% YC) - f0^2) / VARy
	# This is Eq (4.22) in Saltelli et al. (2008)
    
    robj <- rbind(Si = Si, STi = STi)
    colnames(robj) <- paste("X", 1:M, sep ='')
    return(robj)
    
}