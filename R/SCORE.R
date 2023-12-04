#' @title Spectral Clustering On Ratios-of-Eigenvectors.
#' @description Using ratios-of-eigenvectors to detect underlying communities.
#' @details \code{SCORE} is fully established in \emph{Fast community detection by
#'   SCORE} of Jin (2015). \code{SCORE} uses the entrywise ratios between the
#'   first leading eigenvector and each of the other \eqn{K-1} leading eigenvectors for
#'   clustering. It is noteworthy that SCORE only works on connected graphs.
#'   In other words, it does not allow for isolated vertices.
#' @param G An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive
#'   entries being \eqn{1}, where isolated nodes are not allowed.
#' @param K A positive integer which is no larger than \eqn{n}. This is the predefined number of communities.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. The number of times the algorithm should be run with different initial
#'   centroids. The default value is 10.
#'
#' @return \item{estall}{A factor indicating nodes' labels. Items sharing the same label are in the same community.}
#'
#' @importFrom stats kmeans runif
#'
#' @references Jin, J. (2015). Fast community detection by score.
#'   \emph{The Annals of Statistics}, 43 (1), 57–89.\cr\doi{10.1214/14-AOS1265}\cr
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
#' library(igraph)
#' is.igraph(Adj) # [1] FALSE
#' ix = components(graph.adjacency(Adj))
#' componentLabel = ix$membership
#' giantLabel = which(componentLabel == which.max(ix$csize))
#' Giant = Adj[giantLabel, giantLabel]
#' SCORE(Giant, 2)
#'
#' @export

SCORE = function(G, K, itermax = 100, startn = 10){
  # Inputs:
  # 1) G: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1. No isolated nodes.
  # 2) K: a positive integer which is no larger than n. This is the predefined number of communities.

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  # Remark:
  # SCORE only works on connected graphs, i.e., no isolated node is allowed.

  # exclude all wrong possibilities:
  if(!isSymmetric(G)) stop("Error! G is not symmetric!")
  if (any(G != 0 & G != 1)) stop("Error! Adjacency matrix contains values other than 0 or 1!")
  if (any(diag(G) != 0)) stop("Error! Adjacency matrix diagonals are not all zeros.")
  if(K > dim(G)[1]) stop("Error! More communities than nodes!")
  if(dim(G)[1] == 1) stop("Error! There is only one node. No need for clustering!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")

  g.eigen = eigen(G)
  if(sum(g.eigen$vectors[, 1]==0) > 0) stop("Error! Zeroes in the first column")
  R = g.eigen$vectors[, -1]
  R = R[, 1: (K-1)]
  R = R / g.eigen$vectors[, 1]
  n = dim(G)[1]
  R[R > sqrt(log(n))] = sqrt(log(n))
  R[R < -1*sqrt(log(n))] = -1*sqrt(log(n))

  # apply Kmeans to assign nodes into communities
  result = kmeans(R, K, iter.max = itermax, nstart = startn)


  estall = as.factor(result$cluster)
  return(estall)
}
