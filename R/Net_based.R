#' @title Network-based Regularized Spectral Clustering.
#' @description \emph{Network-based Regularized Spectral Clustering} is a spectral clustering with regularized
#'   Laplacian method, fully established in fully established in \emph{Impact of Regularization on Spectral Clustering}
#'   of Joseph & Yu (2016).
#' @param Adj An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries
#'   being \eqn{1}.
#' @param tau An optional tuning parameter to add \eqn{J} to the adjacency matrix \eqn{A}, where \eqn{J} is a
#'   constant matrix with all entries equal to \eqn{1/n}. The default value is the mean of nodes' degrees.
#' @param K A positive positive integer which is no larger than \eqn{n}. This is the predefined number
#'   of communities.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. The number of times the algorithm should be run with different initial
#'   centroids. The default value is 10.
#'
#' @return \item{estall}{A factor indicating nodes' labels. Items sharing the same label are in the same community.}
#'
#' @importFrom stats kmeans runif
#'
#' @references Joseph, A., & Yu, B. (2016). Impact of Regularization on Spectral Clustering.
#'   \emph{The Annals of Statistics}, 44(4), 1765-1791. \cr\doi{10.1214/16-AOS1447}\cr
#'
#' @examples
#'
#' # Simulate the Network
#' n = 10; K = 2;
#' theta = 0.4 + (0.45-0.05)*(seq(1:n)/n)^2; Theta = diag(theta);
#' P  = matrix(c(0.8, 0.2, 0.2, 0.8), byrow = TRUE, nrow = K)
#' set.seed(2022)
#' l = sample(1:K, n, replace=TRUE); # node labels
#' Pi = matrix(0, n, K) # label matrix
#' for (k in 1:K){
#'   Pi[l == k, k] = 1
#' }
#' Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;
#' Adj = matrix(runif(n*n, 0, 1), nrow = n);
#' Adj = Omega - Adj;
#' Adj = 1*(Adj >= 0)
#' diag(Adj) = 0
#' Adj[lower.tri(Adj)] = t(Adj)[lower.tri(Adj)]
#' Net_based(Adj, 2)
#' @export

Net_based <- function(Adj, K, tau = NULL, itermax = 100, startn = 10){

  if(!isSymmetric(Adj)) stop("Error! Adjacency matrix is not symmetric!")
  if(any(Adj != 0 & Adj != 1)) stop("Error! Adjacency matrix contains values other than 0 or 1!")
  if(any(diag(Adj) != 0)) stop("Error! Adjacency matrix diagonals are not all zeros.")
  if(dim(Adj)[1] == 1) stop("Error! There is only one node. No need for clustering!")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")

  if(is.null(tau)) tau = mean(colSums(Adj));

  n <- dim(Adj)[1]
  A_tau = Adj + tau * matrix(1, n, n)/n
  s = rowSums(A_tau)
  s =  s^(-1/2)
  S = diag(s)
  Z = S %*% A_tau %*% S
  g.eigen <-  eigen(Z)
  R = g.eigen$vectors
  R = R[, 1: K]
  R <- t(apply(R, 1, function(x) x/sqrt(sum(x^2))))

  # apply Kmeans to assign nodes into communities
  result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix

  estall = as.factor(result$cluster)
  return(estall)
}
