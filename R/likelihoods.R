# likelihood calculation for MH step
#'
#' @description calculates likelihood for thetaT
#'
#' @param n sample size
#' @param R correlation matrix
#' @param j value of quadratic inner component of normal likelihood
#'
#' @return traceplots
#'
#' @examples
#' example(fdelt(n = 100,R = diag(c(1,1,1)), j = 10)
}
)
fdelt = function(n,R,j){ 
	return(-(n/2)*log(det(R))+(-0.5*j) )
	}
