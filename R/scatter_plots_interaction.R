#' scatter plots of x(i) against all other x(j)
#'
#' This function produces scatter plots of input \eqn{x(i)}  against all other inputs \eqn{x(j) (j = 1,...,M;  j! = i)}, where the color of the marker is proportional to the value of the model output \code{y}
#'
#' @param X matrix \code{(N, M)} of \code{N} inputs samples
#' @param Y vector \code{(N)} of associated ouput samples
#' @param col vector of length two of column interactions. Default \code{col = NULL}, all the possible combinations are plotted.  
#' @param ... parameters to be passed in the plot see \code{\link{plot}}

#' @export

#' @examples

#' #############################
#' # Step 1 (setup the model)
#' ############################
#' fun_test  <- "ishigami_homma_function"
#' M <- 3
#' distr_fun <- "unif"
#' distrpar <-  c(-pi, pi)
#' # ############################
#' # Step 2 (sampling and model evaluation)
#' # ############################
#' N <- 3000
#' X <- AAT_sampling("lhs", M, distr_fun, distrpar, N)
#' Y <- model_evaluation(fun_test, X)         
#' # ############################
#' # Step 3 (Scatter plots)
#' # ############################
#' # Use coloured scatter plots of one input against another on to assess
#' # interactions:
#' # plot x(i1) against x(i3)
#' dev.new()
#' scatter_plots_interaction(X, Y, col = c(1, 3))
#' # Put all possible combinations into one figure:
#' # Customize titles:
#' dev.new()
#' colnames(X) <- c("x(1)", "x(2)", "x(3)")
#' scatter_plots_interaction(X, Y)

scatter_plots_interaction <- function(X, Y, col = NULL, ...){

	stopifnot(is.matrix(X), is.numeric(X),
		is.numeric(Y),
		nrow(X) == length(Y))
		
	M <- ncol(X)	
	
	if(is.null(colnames(X))) colnames(X) <- paste("#", seq(1, M), sep ="")
	
	if(!is.null(col)){ 
		stopifnot(length(col) == 2, col <= M)	
		
		par(mar = c(5.1,4.1,4.1,10))
		
		nb <- 7
	
	Ybreak <- cut(Y, breaks = nb)
	colpoints <- rev(rainbow(nb, alpha = .7))[as.numeric(Ybreak)]
		
		plot(X[, col], pch = "o", col = colpoints, main = paste("(", colnames(X)[col[1]], " vs ", colnames(X)[col[2]], ")", sep =""), xlab =colnames(X)[col[1]], ylab =colnames(X)[col[2]])	
		
		legend(max(X[,col]) + .5, quantile(X[,col], .75) , legend = levels(Ybreak), lwd = 2, col = rev(rainbow(nb, alpha = .7)), xpd = TRUE, title = "Y ranges")
		
		par(mar=c(5.1,4.1,4.1,2.1))
		
		} else {
	
	par(mfrow = c(M - 1, M - 1))

	nb <- 7
	
	Ybreak <- cut(Y, breaks = nb)
	colpoints <- rev(rainbow(nb, alpha = .7))[as.numeric(Ybreak)]

	
	for(i in 1:(M-1)){
		for (j in 1:M){
			if (j!=i){
				
				if(j < i){	

					if(j == 1 & i  == 2){ 
						plot(1, type="n", axes=F, xlab="", ylab="", main = "Legend -- Y ranges")
						legend('center', legend = levels(Ybreak), lwd = 2, col = rev(rainbow(nb, alpha = .7)))
						} else {
							plot(1, type="n", axes=F, xlab="", ylab="")
						}
					
					} else{
			plot(X[, i], X[, j], pch = "o", col = colpoints, main = paste("(", colnames(X)[i], " vs ", colnames(X)[j], ")", sep =""), xlab =colnames(X)[i], ylab =colnames(X)[j])	
				}
			}
			}
		}
	
	par(mfrow = c(1,1))
	}
}