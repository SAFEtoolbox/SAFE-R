#' Hymod Nash-Sutcliffe Efficiency
#'
#' This function runs the rainfall-runoff Hymod model and returns the associated Nash-Sutcliffe Efficiency

#'
#' @param x vector \code{5} of model parameters \code{(Smax, beta, alfa, Rs, Rf)}
#' @param dat dataset containing \code{rain} vector \code{T} time series of rainfall, \code{evap} vector \code{T} time series of potential evaporation, \code{flow} vector \code{T} time series of observed flow.  

#'@return \code{y}  Nash-Sutcliffe Efficiency 

#' @seealso \code{\link{hymod_sim}}

#' @export

hymod_nse <- function(x, dat){

Qsim <- hymod_sim(dat$rain, dat$evaporation, x)

warmup <- 30 #  warmup period to be discarded

Qs <- Qsim[-(1:warmup)]

Qo <- dat$flow[-(1:warmup)]

y <- 1 - var(Qs - Qo) / var(Qo)

return(y)

}