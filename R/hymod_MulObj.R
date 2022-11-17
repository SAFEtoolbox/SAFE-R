#' Hymod RMSE and BIAS
#'
#' This function runs the rainfall-runoff Hymod model and returns 2 metrics of model performance: RMSE and BIAS

#' @param x vector \code{5} of model parameters \code{(Smax, beta, alfa, Rs, Rf)}
#' @param dat dataset containing \code{rain} vector \code{T} time series of rainfall, \code{evap} vector \code{T} time series of potential evaporation, \code{flow} vector \code{T} time series of observed flow.  

#'@return \code{Y}  vector \code{(2)} of objective functions \code{(RMSE, BIAS)}

#' @seealso \code{\link{hymod_sim}}

#' @export

hymod_MulObj <- function(x, dat){

Qsim <- hymod_sim(dat$rain, dat$evaporation, x)

warmup <- 30 #  warmup period to be discarded

Qs <- Qsim[-(1:warmup)]

Qo <- dat$flow[-(1:warmup)]

Y <- numeric(2)
Y[1] <- sqrt(mean((Qs - Qo)^2)) # RMSE
Y[2] <- abs(mean(Qs - Qo)) # BIAS

return(Y)

}