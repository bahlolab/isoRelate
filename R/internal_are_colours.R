#' Internal function
#'
#' \code{areColors()} checks if colour names are valid
#' @param X vector of length 1 or higher containing numeric or character values for colours
areColors <- function(x) {
  vapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  }, FUN.VALUE=FALSE)
}
