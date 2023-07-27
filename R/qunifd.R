qunifd <- function(p, min = 0, max = 1){
	stopifnot(all(p >= 0), all(p <= 1))
	ceiling(p * (max - min + 1))
}