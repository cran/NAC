#' @title projAXb
#' @param X0 The ground truth clustering matrix.
#' @param k The number of community.
#' @param n The number of nodes.
#' @return X The solution matrix of the SDP.
#' @noRd
#' @keywords internal

projAXb = function(X0, k, n){
  # k is the trace of X
  b = c(2*rep(1, n), k)
  X = X0 - Acs(Pinv(Ac(X0, n) - b, n), n)
  return(X)
}
