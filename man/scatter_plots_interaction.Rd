% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scatter_plots_interaction.R
\name{scatter_plots_interaction}
\alias{scatter_plots_interaction}
\title{scatter plots of x(i) against all other x(j)}
\usage{
scatter_plots_interaction(X, Y, col = NULL, ...)
}
\arguments{
\item{X}{matrix \code{(N, M)} of \code{N} inputs samples}

\item{Y}{vector \code{(N)} of associated ouput samples}

\item{col}{vector of length two of column interactions. Default \code{col = NULL}, all the possible combinations are plotted.}

\item{...}{parameters to be passed in the plot see \code{\link{plot}}}
}
\description{
This function produces scatter plots of input \eqn{x(i)}  against all other inputs \eqn{x(j) (j = 1,...,M;  j! = i)}, where the color of the marker is proportional to the value of the model output \code{y}
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
# Step 3 (Scatter plots)
# ############################
# Use coloured scatter plots of one input against another on to assess
# interactions:
# plot x(i1) against x(i3)
dev.new()
scatter_plots_interaction(X, Y, col = c(1, 3))
# Put all possible combinations into one figure:
# Customize titles:
dev.new()
colnames(X) <- c("x(1)", "x(2)", "x(3)")
scatter_plots_interaction(X, Y)
}
