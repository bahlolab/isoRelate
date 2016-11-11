#' Internal Function
#'
#' \code{merge_lists()} is a function used to merge summary IBD results for multiple pairs when running the IBD analysis on
#' multiple cores
#'
#' @param A List with 2 objects for one pair. Object 1 is the IBD summaries for a pair and object 2 is the IBD posterior probabilities.
#' @param B List with 2 objects for another pair. Object 1 is the IBD summaries for a pair and object 2 is the IBD posterior probabilities.
#' @return A list with 2 objects containing merged lists from \code{A} and \code{B}    above.
merge_lists <- function(A,B){
  X1 <- mapply(rbind, A[1], B[1], SIMPLIFY=FALSE)[[1]]
  X2 <- A[[2]] + B[[2]]
  return(list(X1,X2))
}


#' @useDynLib isoRelate
#' @importFrom Rcpp sourceCpp
NULL
