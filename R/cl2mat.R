#' @title cl2mat
#' @param labels A vector, representing the labels of nodes.
#' @return mat The result matrix.
#' @noRd
#' @keywords internal

cl2mat = function(labels){
  n = length(labels);
  k = length(unique(labels));
  mat = matrix(0, n, k)
  for (j in 1: k){
    mat[, j] =  1*(labels == j);
  }
  return(mat)
}
