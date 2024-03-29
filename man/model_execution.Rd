% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_execution.R
\name{model_execution}
\alias{model_execution}
\title{Model execution}
\usage{
model_execution(fun_test, X, ...)
}
\arguments{
\item{fun_test}{name of the function implementing the model (string)}

\item{X}{matrix \code{(N, M)} of \code{N} sampled input factors}

\item{...}{other parameters to be passed in \code{fun_test}}
}
\value{
\code{Y} matrix \code{(N, P)} of associated model ouputs (\code{P} being the number of scalar model outputs associated to each sampled input combination), \code{comp_time} total computing time for model execution (sec)
}
\description{
This function executes the model
}
