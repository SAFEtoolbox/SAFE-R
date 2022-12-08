#' Ishigami Homma function
#'
#' This function does  Ishigami Homma
#'
#' @param x vector \code{(M)} of inputs

#' @return \code{Y} matrix \code{(N, P)} of associated model ouputs (\code{P} being the number of scalar model outputs associated to each sampled input combination), \code{comp_time} total computing time for model evaluation (sec)  

#' @export

#' @examples
#' # See the demo
#' # demo("workflow_visual_ishigami_homma")
#' # or 
#' # demo("workflow_vbsa_ishigami_homma")

ishigami_homma_noise_function <- function(x){
  
  a <- 2
  b <- 1
  y <- sin(x[1]) + a * sin(x[2])^2 + b * x[3]^4 * sin(x[1]) + rnorm(1,0,20)
  
  # # By model definition, we should get:
  # # VARy = VAR(Y) = 1/2 + a^2/8 + b*pi^4/5 + b^2*pi^8/18
  # # V1 = VAR(E(Y|X1)) = 1/2 + b*pi^4/5 + b^2*pi^8/50
  # # V2 = VAR(E(Y|X2)) = a^2/8
  # # V3 = VAR(E(Y|X3)) = 0
  # # V13 = VAR(E(Y|X1,X3)) = b^2*pi^4/18 - b^2*pi^8/50
  # # V12 = V23 = V123 = 0
  # # and thus:
  # # ST1 = S1 + S13
  # # ST2 = S2
  # # ST3 = S13
  
  V <- 1 / 2 + a^2 / 8 + b * pi^4 / 5 + b^2 * pi^8 / 18 + rnorm(1,0,20)
  
  Si_ex <- numeric(3)
  
  Si_ex[1] <- (1 / 2 + b * pi^4 / 5 + b^2 * pi^8 / 50) / V
  Si_ex[2] <- a^2 / 8 / V
  
  STi_ex <- numeric(3)
  STi_ex[1] <- Si_ex[1] + (b^2 * pi^8 / 18 - b^2 * pi^8 / 50) / V
  STi_ex[2] <- Si_ex[2]
  STi_ex[3] <- (b^2 * pi^8 / 18 - b^2 * pi^8 / 50) / V
  
  attributes(y) <- list(V = V, Si_ex = Si_ex, STi_ex = STi_ex)
  
  return(y)
  
}