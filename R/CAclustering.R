#' @title Covariate Assisted Spectral Clustering.
#' @description \emph{CAclustering} clusters graph nodes by applying spectral clustering with the assistance from
#'   node covariates.
#' @details \code{CAclustering} is an algorithm designed for community detection in networks with node covariates,
#'   as introduced in the paper \emph{Covariate-assisted spectral clustering} of Binkiewicz et al. (2017).
#'   \code{CAclustering} applies
#'   \code{k-means} on the first \eqn{K} leading eigenvectors of \eqn{L_{\tau}+\alpha XX^{\prime}}, where
#'   \eqn{L_{\tau}} is the regularized graph Laplacian, \eqn{X} is the covariates matrix, and \eqn{\alpha} is
#'   a tuning parameter.
#' @param Adj An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries
#'   being \eqn{1}.
#' @param Covariate An \eqn{n \times p} covariate matrix. The rows correspond to nodes and the columns correspond
#'   to covariates.
#' @param K A positive integer which is no larger than \eqn{n}. This is the predefined number of communities.
#' @param alphan The number of candidate \eqn{\alpha}'s to try within the range \eqn{(\alpha_{min}, \alpha_{max})} given in
#'  Binkiewicz et al. (2017). An optimal \eqn{\alpha} is expected to achieve a balance between \eqn{L_{\tau}}
#'  and \eqn{X}.
#' @param itermax \code{k-means} parameter, indicating the maximum number of iterations allowed.
#'   The default value is 100.
#' @param startn \code{k-means} parameter. The number of times the algorithm should be run with different initial
#'   centroids. The default value is 10.
#'
#' @return \item{estall}{A factor indicating nodes' labels. Items sharing the same label are in the same community.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif
#'
#' @references Binkiewicz, N., Vogelstein, J. T., & Rohe, K. (2017). Covariate-assisted spectral clustering.
#'   \emph{Biometrika}, 104(2), 361-377. \cr\doi{10.1093/biomet/asx008}\cr
#' @examples
#'
#' # Simulate the Network
#' n = 10; K = 2; p =5; prob1 = 0.9;
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
#' Q = 0.1*matrix(sign(runif(p*K) - 0.5), nrow = p);
#' for(i in 1:K){
#'   Q[(i-1)*(p/K)+(1:(p/K)), i] = 0.3; #remark. has a change here
#' }
#' W = matrix(0, nrow = n, ncol = K);
#' for(jj in 1:n) {
#'   pp = rep(1/(K-1), K); pp[l[jj]] = 0;
#'   if(runif(1) <= prob1) {W[jj, 1:K] = Pi[jj, ];}
#'   else
#'   W[jj, sample(K, 1, prob = pp)] = 1;
#'   }
#' W = t(W)
#' D0 = Q %*% W
#' D = matrix(0, n, p)
#' for (i in 1:n){
#'   D[i,] = rnorm(p, mean = D0[,i], sd = 1);
#' }
#' CAclustering(Adj, D, 2)
#' @export

CAclustering = function(Adj, Covariate, K, alphan = 5, itermax = 100, startn = 10){

  # Inputs:
  # 1) Adj: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) Covariate: an n by p covariate matrix
  # 3) K: a positive integer which is no larger than n. This is the predefined number of communities.
  # 4) alphan: a positive number to tune the weight of covariate matrix

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  if(!isSymmetric(Adj)) stop("Error! Adjacency matrix is not symmetric!")
  if(any(Adj != 0 & Adj != 1)) stop("Error! Adjacency matrix contains values other than 0 or 1!")
  if(any(diag(Adj) != 0)) stop("Error! Adjacency matrix diagonals are not all zeros.")
  if(K > dim(Adj)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")
  if(dim(Adj)[1] != dim(Covariate)[1]) stop("Error! Incompatible!")
  if(!all(is.numeric(Covariate))) stop("Error! Covariate must contain only numeric values!")
  if (!missing(alphan) && (alphan < 1 || !is.integer(alphan))) {
    stop("Error! alphan should be a positive integer!")
  }

  s = rowSums(Adj)
  s = s + mean(s)
  s =  s^(-1/2)
  S = diag(s)
  Z = S %*% Adj %*% S
  net.eigen = eigen(Z%*%Z)
  ca = Covariate %*% t(Covariate);
  ca.eigen = eigen(ca);
  alphalower = (net.eigen$values[K] - net.eigen$values[K+1])/ca.eigen$values[1];
  alphaupper = net.eigen$values[1]/(ca.eigen$values[K] - ca.eigen$values[K+1]);
  d = rep(0, alphan);
  alpha = seq(alphalower, alphaupper, length.out = alphan);
  est = matrix(0, alphan, dim(Adj)[1])

  for(ii in 1:alphan){
    casc.eigen = eigen(Z%*%Z + alpha[ii]*ca);
    U = casc.eigen$vectors[,1:K];
    Unorm = apply(U, 1, Norm);
    indu = which(Unorm > 0);
    U = U[indu, ]/Unorm[indu]
    result = kmeans(U, K, iter.max = itermax, nstart = startn);
    d[ii] = result$tot.withinss;
    est[ii, indu] = as.factor(result$cluster)
  }
  estall = est[which.min(d), ]
  return(estall)
}
