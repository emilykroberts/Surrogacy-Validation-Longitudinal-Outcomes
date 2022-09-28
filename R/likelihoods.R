
# likelihood calculation for MH step
fdelt = function(n,R,j){ return(-(n/2)*log(det(R))+(-0.5*j) )}
