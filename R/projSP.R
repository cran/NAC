#' @title projSP
#' @param X0 The ground truth clustering matrix.
#' @return \item{X}{The result matrix}.
#' @noRd
#' @keywords internal

projSp = function(X0){
  n = dim(X0)[1]
  tryCatch({print("no error");temp = eigen(0.5*(X0+t(X0))); U = temp$vectors; D = diag(temp$values)}, error = function(x)
  {print("error");temp  = svd(0.5*(X0+t(X0))); U = temp$u; V = temp$v;S = diag(temp$d); D = diag(diag(t(U) %*% X0 %*% U))})
  idx = as.integer(which(diag(D) >= 0))
  X = as.matrix(U[, idx]) %*% as.matrix(D[idx,idx]) %*% as.matrix(t(U[, idx]));
  return(X)
}
