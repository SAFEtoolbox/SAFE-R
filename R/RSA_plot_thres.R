#' Plotting function for Regional Sensitivity Analysis
#'
#' This function plots Regional Sensitivity Analysis
#'
#' @param X matrix \code{(N, M)} set of input samples
#' @param idxb vector \code{(N)} indices of samples statisfying the condition
#' @param n_col scalar number of panels per row in the plot (default: \code{min(5,M)})
#' @param prnam vector \code{(M)} labels for the horizontal axis (default: \code{c("X1", "X2",...)})
#' @param str_legend vector (2) text for legend (default: no legend)

#' @seealso \code{\link{RSA_indices_thres}} \code{\link{RSA_convergence_thres}}

#' @export

#' @examples

#' # See the demo
#' # demo("workflow_rsa_hymod")

RSA_plot_thres <- function(X, idxb, prnam = NULL, threshold, str_legend = NULL) {
  
  stopifnot(is.matrix(X), is.numeric(X), is.logical(idxb),
            length(idxb) == nrow(X))
  
  Ng <- 2
  
  N <- nrow(X)
  M <- ncol(X)
  
  if(!is.null(prnam)){
    stopifnot(length(prnam) == M)
  } else {
    prnam <- paste("X", seq(1, M), sep ="")
  }
  
  # Define below and above subsamples:
  idxl <- as.numeric(idxb)
  idx <- idxl + 1
  Yk <- c(min(Y),threshold,max(Y))
  
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
    facet_grid(. ~ parnam, scales = "free_x") + geom_path(size = 1) + geom_point(size = 1) +
    scale_color_manual(values=c("grey","red")) + theme_bw() + theme(legend.position = "none")
                                                          
                                                          
  return( .pl )
}