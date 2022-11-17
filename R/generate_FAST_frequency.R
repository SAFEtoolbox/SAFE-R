#' Generate FAST frequency
#'
#' Generate frequency set free of interferences through (at least) 4th order (vector (\code{M}))
#'
#' For \code{M > 4}, frequencies are computed based on the recursive algorithm by Cukier et al. (1975) which is free of interferences through the 4th order. 
#' For \code{M <= 4}, we use values from the literature that guarantee higher order interferences free: if \code{M = 2} use values from Sec. 3.1 in Xu, C. and G. Gertner (2007) (free of interference through 10th order), if \code{M = 4} use values from Table III in Cukier et al. (1975) (free of interferences through 6th order).
#'
#' @param M integer scalar between 4 and 50, number of inputs
#'
#' @references Cukier et al. (1975) Study of the sensitivity of coupled reaction systems to uncertainties in rate coefficients. III. Analysis of the  approximations, J. Chem. Phys. 63, 1140
#'
#' Xu, C. and G. Gertner (2007), Extending a global sensitivity analysis technique to models with correlated parameters, Computational Statistics and Data Analysis, 51, 5579-5590.

#' @export

#' @examples

#' M <- 5
#' generate_FAST_frequency(M)

generate_FAST_frequency <- function(M){

stopifnot(is.scalar(M), M > 1, M < 51, M == floor(M))

if (M == 2){# Use values from Sec. 3.1 in Xu, C. and G. Gertner (2007)
	# (free of interference through 10th order)

    omega <- c(5, 23)
     
     } else {
    	
    if (M == 4) {
    		# Use values from Table III in Cukier et al. (1975)
    		# (free of interferences through 6th order)
    omega <- c(13, 31, 37, 41)
    
    } else {# Use recursive algorithm in the same paper
    Omega <- c(0, 0, 1, 5, 11, 1, 17, 23, 19, 25, 41, 31, 23, 87, 67, 73, 85, 143, 149, 99, 119,
237, 267, 283, 151, 385, 157, 215, 449, 163, 337, 253, 375, 441, 673, 773, 875, 873,
587, 849, 623, 637, 891, 943, 1171, 1225, 1335, 1725, 1663, 2019)    
        
    d <- c(4, 8, 6, 10, 20, 22, 32, 40, 38, 26, 56, 62, 46, 76, 96, 60, 86, 126, 134, 112, 
92, 128, 154, 196, 34, 416, 106, 208, 328, 198, 382, 88, 348, 186, 140, 170, 284, 
568, 302, 438, 410, 248, 448, 388, 596, 216, 100, 488, 166)
    
    # above values taken from Table VI  
    
    omega <- c(Omega[M], cumsum(d[(M - 1):1]) + Omega[M])
    
    # equation (5.1)
   
   }
   }

return(omega)

}
