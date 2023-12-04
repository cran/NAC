#' @title Acs
#' @param z The input matrix.
#' @param n The number of nodes.
#' @return Z The latent class memberships.
#' @noRd
#' @keywords internal

Acs = function(z, n){
  #mu = z[, 1][1:n]
  mu = z[1:n]
  U = matrix(rep(mu, n), byrow = TRUE, nrow = n)
  V = matrix(rep(mu, n), byrow = FALSE, nrow = n)
  #Z = U + V + z[, 1][n + 1]*diag(n)
  Z = U + V + z[n + 1]*diag(n)
  return(Z)
}
