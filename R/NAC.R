#' @title Spectral Clustering on Network-Adjusted Covariates.
#' @description Using network-adjusted covariates to detect underlying communities.
#' @details \emph{Spectral Clustering Network-Adjusted Covariates (NAC)} is fully established in
#'   \emph{Network-Adjusted Covariates for Community Detection}
#'   of Hu & Wang (2023). This method is particularly effective in the analysis of multiscale networks with
#'   covariates, addressing the challenge of misspecification between networks and covariates. \code{NAC}
#'   relies on the construction of network-adjusted
#'   covariate vectors \eqn{y_i = \alpha_ix_i + \sum_{j: A_{ij}=1}x_j, i \in 1, \cdots, n}, where the first part has
#'   the nodal covariate information and the second part conveys network information.
#'   By constructing \eqn{Y = (y_1, \cdots, y_n)^{\prime} = AX + D_{\alpha}X} where \eqn{A} is the adjacency matrix,
#'   \eqn{X} is the covariate matrix, and \eqn{D_{\alpha}} is the diagonal matrix with diagonals
#'   as \eqn{\alpha_1, \cdots, \alpha_n}, \code{NAC} applies \code{K-means} on the first \eqn{K} normalized left
#'   singular vectors, treating each row as a data point.
#'   A notable feature of \code{NAC} is its tuning-free nature, where node-specific coefficient \eqn{\alpha_i} is
#'   computed given the \eqn{i}-th node's degree. \code{NAC} allows for
#'   user-specified \eqn{\alpha_i} as well. A generalization with uninformative covariates is considered by adjusting
#'   parameter \eqn{\beta}. As long as the covariates do provide information, the specification of \eqn{\beta} can be
#'   ignored.
#' @param Adj An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries being
#'   \eqn{1}.
#' @param Covariate An \eqn{n\times p} covariate matrix. The rows correspond to nodes and the columns correspond to
#'   covariates.
#' @param K A positive integer which is no larger than \eqn{n}. This is the predefined number of communities.
#' @param alpha An optional numeric vector to tune the weight of covariate matrix. The default value is
#'   \eqn{\dfrac{\bar{d}/2}{d_i/\text{log}n+1}}, where \eqn{d_i} is the degree of node \eqn{i} and \eqn{\bar{d}} is the
#'   average degree.
#' @param beta An optional parameter used when the covariate matrix \eqn{X} is uninformative. By default, \eqn{\beta}
#'   is set as 0 assuming \eqn{X} carries meaningful information. Otherwise, users can manually specify a positive value
#'   to weigh network information.
#' @param itermax \code{k-means} parameter, indicating the maximum number of iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. The number of times the algorithm should be run with different initial
#'   centroids. The default value is 10.
#'
#' @return \item{estall}{A factor indicating nodes' labels. Items sharing the same label are in the same community.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif median
#'
#' @references Hu, Y., & Wang, W. (2023). Network-Adjusted Covariates for Community Detection.
#'   \emph{arXiv preprint arXiv:2306.15616}. \cr\url{https://arxiv.org/abs/2306.15616}\cr
#'
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
#' NAC(Adj, D, 2)

#' @export

NAC = function(Adj, Covariate, K, alpha = NULL, beta = 0, itermax = 100, startn = 10){
  # Inputs:
  # 1) Adj: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) Covariate: an n by p covariate matrix
  # 3) K: a positive integer which is no larger than n. This is the predefined number of communities.
  # 4) alpha: a positive number to tune the weight of covariate matrix
  # 5) beta: a non-negative number of weigh network information when covariates are uninformative.

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
  if(!is.null(alpha)){
    if(!all(is.numeric(alpha))) stop("Error! alpha must contain only numeric values!")
    if(length(alpha) != dim(Adj)[1]) stop("Error! Incorrect length of alpha!")
  }
  if(!all(is.numeric(Covariate))) stop("Error! Covariate must contain only numeric values!")
  if(!is.numeric(beta)) stop("Error! beta much be numeric!")
  if(length(beta) != 1) stop("Error! beta should contain only one value!")
  #  if(alpha < 0) stop("Negative Alpha")

  #Regularity check
  estall = rep(NA, dim(Adj)[1]);
  netrow = rowSums(Adj);
  covrow = rowSums(abs(Covariate));
  ind_reg = which(netrow != 0 | covrow != 0)
  Adj = Adj[ind_reg, ind_reg];
  Covariate = Covariate[ind_reg, ];

  ##Algorithm
  n = dim(Adj)[1]; p = dim(Covariate)[2]
  d = rowSums(Adj);
  X = Adj %*% Covariate

  lambda = log(n)/(d + log(n));

  if (is.null(alpha)) {
    alpha = mean(d)/2
    D_alpha = alpha*diag(lambda)
  }
  else{
    D_alpha = diag(alpha)
  }

  Newmat = X + D_alpha%*%Covariate;
  zz = Newmat%*%t(Newmat);
  if(beta != 0){
    AA = Adj%*%Adj;
    zz = zz +  beta*n*AA;
  }

  c = eigen(zz)

  vec = c$vectors
  vecn = vec[,1:K]/apply(vec[,1:K], 1, Norm);
  result = kmeans(vecn, K, iter.max = itermax, nstart = startn);
  if (result$ifault==4) { result = kmeans(X, K,  iter.max = itermax, nstart = startn, algorithm="Lloyd"); }
  est = as.factor(result$cluster);

  estall[ind_reg] = est;

  return(estall)
}
