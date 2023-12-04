#' @title Ac
#' @param X A square matrix.
#' @param n The number of nodes.
#' @return z The result
#' @noRd
#' @keywords internal

Ac = function(X, n){
  z = c(2 * X %*% matrix(1, n, 1), sum(diag(X)))
  return(z)
}
