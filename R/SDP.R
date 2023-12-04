#' @title Semidefinite programming for Community Detection in Networks with Covariates.
#' @description Semidefinite programming (SDP) for optimizing the inner product between combined network and the
#'   solution matrix.
#' @details \code{SDP} is proposed in \emph{Covariate Regularized Community Detection in Sparse Graphs}
#'   of Yan & Sarkar (2021). This method relies on semidefinite programming relaxations for detecting
#'   the community structure in sparse networks with covariates.
#' @param Adj An \eqn{n \times n} symmetric adjacency matrix with diagonals being \eqn{0} and positive entries
#'   being \eqn{1}.
#' @param Covariate An \eqn{n \times p} covariate matrix. The rows correspond to nodes and the columns correspond to
#'   covariates.
#' @param lambda A tuning parameter to weigh the covariate matrix.
#' @param K A positive integer which is no larger than \eqn{n}. This is the predefined number of communities.
#' @param alpha The element-wise upper bound in the \code{SDP}.
#' @param rho The learning rate of \code{SDP}.
#' @param TT The maximum of iteration.
#' @param tol The tolerance for stopping criterion.
#' @param quiet An optional input, indicating whether to print result at each step.
#' @param report_interval An optional input. The frequency to print intermediate result.
#' @param r An optional input. The expected rank of the solution, leave \code{NULL} if no constraint is required.
#' @return \item{estall}{A factor indicating nodes' labels. Items sharing the same label are in the same community.}
#'
#' @importFrom pracma eig Norm
#' @importFrom stats kmeans runif rnorm
#'
#' @references Yan, B., & Sarkar, P. (2021). Covariate Regularized Community Detection in Sparse Graphs.
#'   \emph{Journal of the American Statistical Association}, 116(534), 734-745.
#'   \cr\doi{10.1080/01621459.2019.1706541}\cr
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
#' SDP(Adj, D, lambda = 0.2, K = 2, alpha = 0.5, rho = 2, TT = 100, tol = 5)

#' @export

SDP = function(Adj, Covariate, lambda, K, alpha, rho, TT, tol, quiet = NULL, report_interval = NULL, r = NULL){

  # Inputs:
  # 1) Adj: adjacency matrix
  # 2) Covariate: n x p covaraite matrix
  # 3) lambda: tuning parameter between graph and covariates
  # 4) K: number of clusters
  # 5) alpha: elementwise upper bound in the SDP
  # 6) rho: learning rate of ADMM
  # 7) TT:   max iteration
  # 8) tol: tolerance for stopping criterion
  # 9) quiet: whether to print result at each step
  # 10) report_interval: frequency to print intermediate result
  # 11) r: expected rank of the solution, leave blank if no constraint is required.

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
  if(!is.numeric(alpha)) stop("Error! alpha must be numeric!")
  if(!is.numeric(rho)) stop("Error! rho must be numeric!")
  if(TT %% 1 != 0) stop("Error! TT should be a relatively large positive integer!")
  if(TT < 1) stop("Error! TT should be a relatively large positive integer!")

  As = Adj + lambda* Covariate %*% t(Covariate)

  n = dim(As)[1]
  U <- V <- matrix(0, n, n)

  #Initialization - spectral with perturbation
  v = eigen(As)$vectors[, 1: K]
  e = diag(eigen(As)$values[1: K])
  X = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)
  Y = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)
  Z = v %*% t(v) + 0.1*matrix(rnorm(n*n), nrow = n)

  As_rescaled = (1/rho) * As;

  if(is.null(report_interval)) report_interval = 1
  if(is.null(quiet)) quiet = FALSE
  if(is.null(r)) r = Inf

  if (is.infinite(TT)) {
    delta = matrix(0, 1000, 1)
    infeas = matrix(0, 1000, 1)
  }
  else {
    delta = matrix(0, TT, 1)
    infeas = matrix(0, TT, 1)
  }

  dt = matrix(0, 1, 3);
  t = 1;
  CONVERGED = FALSE;

  while(CONVERGED == FALSE & t<= TT){
    Xold = X
    X = projAXb(0.5*(Z - U + Y - V + As_rescaled), K, n);
    Z = X + U
    Z[Z < 0] = 0
    Z[Z > alpha] = alpha

    Y = projSp(X + V);
    U = U + X - Z;
    V = V + X - Y;
    delta[t] = norm(X - Xold) / norm(Xold);
    infeas[t] = (sqrt(sum(diag(t(X - Y) * (X - Y)))) + sqrt(sum(diag(t(X - Z) * (X - Z))))) / sqrt(sum(diag(t(X)*X)));
    CONVERGED = max(delta[t], infeas[t]) < tol;
    t = t + 1;
  }

  T_term = t - 1
  estall = rsc(X, K, 'adj')
  return(estall)
}
