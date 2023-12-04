#' @title Pinv
#' @param z The latent class memberships.
#' @param n The number of nodes.
#' @return \item{X}{The optimal ressult of SDP}.
#' @noRd
#' @keywords internal

Pinv = function(z, n){
  mu = z[1:n]
  nu = z[(n+1): length(z)]
  x1 = 1/2/n * (diag(n) - (n-2)/n/(2*n-2)*matrix(1, n, n)) %*% mu
  x2 = 1/n/(2-2*n)*matrix(1, n, 1) %*% nu
  x3 = - 1/n/(2*n-2)*matrix(1, 1, n) %*% mu
  x4 = nu/(n - 1)
  X = c(x1 + x2, x3 + x4)
  return(X)
}
