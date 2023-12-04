context("Net_based")

test_that("Net_based stops when it should", {
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

test_that("This function returns a list of predicted membership of all nodes", {
  expect_length(Net_based(Adj, K), n)
  expect_length(unique(Net_based(Adj, K)), K)
})
