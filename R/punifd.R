punifd <- function(q, min = 0, max = 1){
	
	if(q < min){
		robj <- 0
		} else {
			robj <- ifelse(q > max, 1,  floor(q) / (max - min + 1))		
		}
		
		return(robj)
		
}