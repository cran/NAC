#' @title normalizeSym
#' @param A An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries being \eqn{1}.
#' @param L An optional input.
#' @return \item{An}{The result matrix}.
#' @noRd
#' @keywords internal

normalizeSym = function(A, L = NULL){
  if (is.null(L)) L = 0
  n = dim(A)[1]
  d = colSums(A)
  d[which(d == 0)] = 1e10
  if (L == 0) {
    d = 1/sqrt(d)
    D = diag(d)
    An = D %*% A %*% D
    An = 1/2*(An + t(An))
  }
  else{
    d = 1 / d
    D=diag(d)
    An=D %*% A;
  }
  return(An)
}
