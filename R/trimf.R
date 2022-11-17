trimf <- function(x, a, b, c){
	stopifnot(c > a, c > b, b > a)
	pmax(pmin((x - a)/ (b - a), (c - x) / (c - b)), 0)
}