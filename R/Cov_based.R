#' @title Covariates-based Spectral Clustering.
#' @description \emph{Covariates-based Spectral Clustering} is a spectral clustering
#'   method that focuses solely on the covariates structure, i.e., the \eqn{XX^{\prime}} where \eqn{X} is the
#'   covariates matrix, as introduced in Lee et al. (2010).
#' @param Covariate An \eqn{n \times p} covariate matrix. The rows correspond to nodes and the columns
#'   correspond to covariates.
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
#' @references Lee, A. B., Luca, D., Klei, L., Devlin, B., & Roeder, K. (2010).
#'   Discovering genetic ancestry using spectral graph theory.
#'   \emph{Genetic Epidemiology: The Official Publication of the International Genetic Epidemiology Society},
#'   34(1), 51-59. \cr\doi{10.1002/gepi.20434}\cr
#'
#' @examples
#'
#' # Simulate the Covariate Matrix
#' n = 10; p = 5; K = 2; prob1 = 0.9;
#' set.seed(2022)
#' l = sample(1:K, n, replace=TRUE); # node labels
#' Pi = matrix(0, n, K) # label matrix
#' for (k in 1:K){
#'   Pi[l == k, k] = 1
#' }
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
#' Cov_based(D, 2)
#' @export

Cov_based <- function(Covariate, K, itermax = 100, startn = 10){

  if(!all(is.numeric(Covariate))) stop("Error! Covariate must contain only numeric values!")
  if(dim(Covariate)[1] == 1) stop("Error! There is only one node. No need for clustering!")
  if(K > dim(Covariate)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K < 2) stop("Error: There must be at least 2 communities!")

  # matrix preparation
  New_X <- Covariate %*% t(Covariate)

  g.eigen = eigen(New_X)
  R = g.eigen$vectors
  R = R[, 1: K]
  # apply Kmeans to assign nodes into communities
  result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix

  estall = as.factor(result$cluster)
  return(estall)
}
