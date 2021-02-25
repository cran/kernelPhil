#' Square kernel
#'
#' Square kernel function - weights all points within bandwidth as 1, all others as 0.
#'
#' @param d Numeric vector of distances.
#' @param bandwidth Numeric bandwidth.
#' @return A numeric vector of kernel weights.
#' @keywords internal
square.kernel<-function(d,bandwidth){as.numeric(abs(d)<=bandwidth);}
#' Gaussian kernel
#'
#' Gaussian kernel function - weights all points by distance on a Gaussian curve with standard deviation of bandwidth.
#'
#' @param d Numeric vector of distances.
#' @param bandwidth Numeric bandwidth.
#' @return A numeric vector of kernel weights.
#' @keywords internal
gaussian.kernel<-function(d,bandwidth){return(exp(0-((d^2)/(2*(bandwidth^2)))));}
#' Gaussian square kernel
#'
#' Gaussian square kernel function - weights all points within bandwidth by distance on a Gaussian curve with standard deviation of bandwidth, but all points beyond bandwidth at 0.
#'
#' @param d Numeric vector of distances.
#' @param bandwidth Numeric bandwidth.
#' @return A numeric vector of kernel weights.
#' @keywords internal
gaussian.square.kernel<-function(d,bandwidth){
	return(exp(0-((d^2)/(2*(bandwidth^2))))*as.numeric(abs(d)<=bandwidth));
}
#' Triangular kernel
#'
#' Triangular kernel function - weights all points within bandwidth on a straight line such that where d=0 weight=1, and where d=bandwidth weight=0; all points beyond bandwidth are weighted 0.
#'
#' @param d Numeric vector of distances.
#' @param bandwidth Numeric bandwidth.
#' @return A numeric vector of kernel weights.
#' @keywords internal
triangular.kernel<-function(d,bandwidth){
	return(((bandwidth-abs(d))/bandwidth)*as.numeric(abs(d)<=bandwidth));
}
