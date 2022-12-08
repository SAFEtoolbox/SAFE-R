qunifd <- function(p, min = 0, max = 1){
	stopifnot(p >= 0, p <= 1)
	ceiling(p * (max - min + 1))
}