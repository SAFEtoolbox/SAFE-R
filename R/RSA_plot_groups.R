#' Plotting function for Regional Sensitivity Analysis with grouping
#'
#' Plotting function for Regional Sensitivity Analysis with grouping. Plot \code{Ng} CDFs of the samples in \code{X} with different colours.

#'
#' @param X matrix \code{(N, M)} set of input samples
#' @param idx vector \code{(N)} index of group to which input samples belong
#' @param Yk vector \code{(Ng + 1)} range of \code{Y} in each group 
#' @param prnam vector \code{(M)} labels for the horizontal axis (default: \code{c("X1", "X2",...)})
#' @seealso \code{\link{RSA_indices_groups}} \code{\link{RSA_plot_thres}}
#' @export
#' @examples
#' # See the demo
#' # demo("workflow_rsa_hymod")

RSA_plot_groups <- function(X, idx, Yk, prnam = NULL) {
  
  stopifnot(is.matrix(X), is.numeric(X),
            is.numeric(Yk),
            length(idx) == nrow(X))
  
  Ng <- length(Yk) - 1
  
  N <- nrow(X)
  M <- ncol(X)
  
  if(!is.null(prnam)){
    stopifnot(length(prnam) == M)
  } else {
    prnam <- paste("X", seq(1, M), sep ="")
  }
  
  xx <- lapply(1:M, function(.ii) { unique(sort(X[, .ii])) })
  
  CDF <- lapply(1:M, function(.ii) 
    lapply(1:Ng, function(.kk) {
      .tmp <- ecdf(X[idx == .kk, .ii])(xx[[.ii]])
      c(rep(.tmp[1:(length(.tmp)-1)], each = 2), tail(.tmp, 1) )
    }))
  
  xx <- lapply(xx, function(.inp) c(.inp[1], rep(.inp[-1], each = 2)))
  
  nux <- sapply(xx, length)
  
  dat <- data.frame(x = unlist(lapply(xx, function(.inp) rep(.inp, Ng))), 
                    CDF = unlist(lapply(CDF, function(.inp) unlist(.inp))),
                    parnam = unlist(lapply(1:M, function(.ii) rep(prnam[.ii], each = Ng * nux[.ii]))),
                    group = unlist(lapply(1:M, function(.ii) rep(round(Yk[2:(Ng+1)],2), each = nux[.ii])))
  )
  dat$group <- factor(dat$group, levels = round(Yk[2:(Ng+1)],2))
  dat$parnam <- factor(dat$parnam, levels = prnam)
  
  .pl <- ggplot(data = dat, mapping = aes(x = x, y = CDF, color = group, parnam = parnam)) + 
    facet_grid(. ~ parnam, scales = "free_x") + geom_line() + theme_bw() +
    labs(color='\t Upper \n bound \n output')
  
  return( .pl )
}