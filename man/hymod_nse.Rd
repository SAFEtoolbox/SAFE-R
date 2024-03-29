% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hymod_nse.R
\name{hymod_nse}
\alias{hymod_nse}
\title{Hymod Nash-Sutcliffe Efficiency}
\usage{
hymod_nse(x, dat)
}
\arguments{
\item{x}{vector \code{5} of model parameters \code{(Smax, beta, alfa, Rs, Rf)}}

\item{dat}{dataset containing \code{rain} vector \code{T} time series of rainfall, \code{evap} vector \code{T} time series of potential evaporation, \code{flow} vector \code{T} time series of observed flow.}
}
\value{
\code{y}  Nash-Sutcliffe Efficiency
}
\description{
This function runs the rainfall-runoff Hymod model and returns the associated Nash-Sutcliffe Efficiency
}
\seealso{
\code{\link{hymod_sim}}
}
