approximatecdfpair <- function(x1,x2){

#     x1 = sample of variable 1  - vector (N,1)
#     x2 = sample of variable 2  - vector (M,1)
#     it returns:  
#	x = datapoints where the empirical CDFs are evaluated (it includes
#          all unique values in x1 and x2)  - vector (P,1)
#  CDF1 = empirical CDF of x1 evaluated at 'x'  - vector (P,1) 
#  CDF2 = empirical CDF of x2 evaluated at 'x'  - vector (P,1)

     # datapoints where the empirical CDFs are evaluated (it includes all unique values in x1 and x2)
      
     x <- unique(sort(c(x1, x2)))
      
      # empirical CDF of x1 and x2
      
      CDF1 <- ecdf(x1)
      CDF2 <- ecdf(x2)
      
      # empirical CDF of x1 and x2 evaluated at 'x'
      CDF1 <- CDF1(x)
      CDF2 <- CDF2(x)
      
      return(list(x = x, CDF1 = CDF1, CDF2 = CDF2))

}
