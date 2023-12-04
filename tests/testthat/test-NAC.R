context("NAC")

test_that("NAC stops when it should", {
  expect_error( runningmean(0, c(0,0)) )
})

library(pracma)
# Simulate the Network
n = 10; K = 2; p =5; prob1 = 0.9;
theta = 0.4 + (0.45-0.05)*(seq(1:n)/n)^2; Theta = diag(theta);
P  = matrix(c(0.8, 0.2, 0.2, 0.8), byrow = TRUE, nrow = K)
set.seed(2022)
l = sample(1:K, n, replace=TRUE); # node labels
Pi = matrix(0, n, K) # label matrix
for (k in 1:K){
 Pi[l == k, k] = 1
}
Omega = Theta %*% Pi %*% P %*% t(Pi) %*% Theta;
Adj = matrix(runif(n*n, 0, 1), nrow = n);
Adj = Omega - Adj;
Adj = 1*(Adj >= 0)
diag(Adj) = 0
Adj[lower.tri(Adj)] = t(Adj)[lower.tri(Adj)]
Q = 0.1*matrix(sign(runif(p*K) - 0.5), nrow = p);
for(i in 1:K){
 Q[(i-1)*(p/K)+(1:(p/K)), i] = 0.3; #remark. has a change here
}
W = matrix(0, nrow = n, ncol = K);
for(jj in 1:n) {
 pp = rep(1/(K-1), K); pp[l[jj]] = 0;
 if(runif(1) <= prob1) {W[jj, 1:K] = Pi[jj, ];}
 else
 W[jj, sample(K, 1, prob = pp)] = 1;
}
W = t(W)
D0 = Q %*% W
D = matrix(0, n, p)
for (i in 1:n){
  D[i,] = rnorm(p, mean = D0[,i], sd = 1);
}

test_that("This function returns a list of predicted membership of all nodes", {
  expect_length(NAC(Adj, D, K), n)
  expect_length(unique(NAC(Adj, D, K)), K)
})
