# Interval functin
#
#' function to find IBD that overlap interval
#' @param region.1 IBD segment boundary
#' @param region.2 interest interval
getOverlap <- function(region.1, region.2) {
  a <- max(region.1[1], region.2[1])
  b <- min(region.1[2], region.2[2])
  return(c(a,b))
}
