#' @title Regularized Spectral Clustering
#' @param A An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries being \eqn{1}.
#' @param K A positive integer which is no larger than \eqn{n}. This is the predefined number of communities.
#' @param method The method of spectral clustering. 'pos' refers to regularized spectral clustering,
#'     'lap' refers to spectral clustering using normalized graph laplacian, and 'adj' refers to using
#'     adjacency matrix.
#' @param prior An optional input. The maximum of iteration.
#' @return \item{class}{A label vector}.
#' @noRd
#' @keywords internal

rsc = function(A, K, method, prior = NULL){
  # Regularized Spectral Clustering
  # Input: A: Adjacency matrix
  #        K: number of clusters
  #        method: 'pos' - Regularized Spectral Clustering
  #                'lap' - Spectral Clustering (use normalized graph laplacian)
  #                'adj' - Use adjacency matrix
  # Output: class labels

  if (is.null(prior)) prior = 1/K * matrix(1, K, 1)
  nv = dim(A)[1]
  tau = mean(A)
  if (method == "pos"){
    A_tau = A + tau * matrix(1, nv, nv)
    L_tau = normalizeSym(A_tau)
  }
  else if(method == "lap") L_tau = normalizeSym(A)
  else L_tau = (A + t(A))/2
  U1 = eigen(L_tau)$vectors[, 1: K]
  U1 = t(apply(U1, 1, function(x) x/sqrt(sum(x^2))))
  maxsum = Inf
  nrestart = 100
  for (i in 1: nrestart){
    clustering = kmeans(U1, K)
    class0 = clustering$cluster
    sumD = clustering$totss
    if (maxsum > sum(sumD)){
      maxsum = sumD
      class=class0
    }
  }
  return(class)
}
