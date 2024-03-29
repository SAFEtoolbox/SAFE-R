% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatter_plots.R
\name{scatter_plots}
\alias{scatter_plots}
\title{Scatter plots of \code{y} against \code{X}}
\usage{
scatter_plots(X, Y, prnam = NULL, ngr = 0)
}
\arguments{
\item{X}{matrix \code{(N, M)} of \code{N} inputs samples}

\item{Y}{vector \code{N} of associated ouput samples}

\item{prnam}{vector of characters containing the parameters names}

\item{ngr}{number of groups into which the parameters are divided}
}
\description{
This function produces scatter plots of the model ouput \code{y} against model inputs \eqn{x(1), x(2), ..., x(M)}
}
\examples{
#############################
# Step 1 (setup the model)
############################
fun_test  <- "ishigami_homma_function"
M <- 3
distr_fun <- "unif"
distrpar <-  c(-pi, pi)
# ############################
# Step 2 (sampling and model execution)
# ############################
N <- 3000
X <- AAT_sampling("lhs", M, distr_fun, distrpar, N)
Y <- model_execution(fun_test, X)         
# ############################
## Step 3 (Scatter plots)
# ############################
# Use scatter plots of inputs againts output to visually assess 
# direct effects:
scatter_plots(X,Y)

}
