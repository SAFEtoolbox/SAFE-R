punifd <- function(q, min = 0, max = 1){
	
  out <- q * 0
  
  for(ii in 1:length(q)){
	 out[ii] <- ifelse(q[ii] > max, 1,  floor(q[ii]) / (max - min + 1))		
  }
		
	return(out)
		
}