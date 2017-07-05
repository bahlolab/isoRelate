#' Internal Function
#'
#' \code{iRfunction()} calculates the iR statistic used to assess the significance of excess IBD sharing at loci
#' in the genome. The final statistic, -log10 (P-values), is returned for each SNP.
#' @param locus.matrix A data frame containing binary IBD values. See \code{\link{getIBDMatrix}} for more details.
#' @param frequency A vector of population allele frequencies for each SNP.

iRfunction <- function(locus.matrix, frequency, plot_qq = FALSE) {
  # subtract the mean sharing from each pair
  y_ijk <- matrix(nrow=nrow(locus.matrix), ncol=ncol(locus.matrix))
  for (i in 1:ncol(locus.matrix)) {
    y_ijk[,i] <- locus.matrix[,i] - mean(locus.matrix[,i])
  }

  # for each SNP, calculate mean transformed value across all pairs and divide by the population allele frequency
  z_ijk <- matrix(nrow=nrow(y_ijk), ncol=ncol(y_ijk))
  for (i in 1:nrow(y_ijk)) {
    z_ijk[i,] <- (y_ijk[i,] - mean(y_ijk[i,]))/sqrt(frequency[i]*(1-frequency[i]))
  }

  # calculate a new summary statistic
  stat_k <- rowSums(z_ijk)/sqrt(ncol(z_ijk))

  # normalize new summary statistic
  stat_norm <- (stat_k - mean(stat_k))/sd(stat_k)
  stat_norm_sq <- stat_norm^2

  # take -log10(p-value)
  p_value <- 1 - pchisq(stat_norm_sq, df = 1)
  log10_p_value <- -log10(p_value)

  # log10_p_value == Inf (results from p-value = 0) are changed to max(log10_p_value)
  log10_p_value[is.infinite(log10_p_value)] <- max(log10_p_value[!is.infinite(log10_p_value)])

  # investigate quantiles - chisquared distirbution
  if(plot_qq == TRUE) {
    qx <- quantile(rchisq(length(stat_norm),df=1),seq(0,1,0.0001))
    qy <- quantile(stat_norm_sq,seq(0,1,0.0001))
    plot(qx,qy,main=country); abline(a=0,b=1,col="grey66")
  }
  return(cbind(stat_norm_sq, log10_p_value))
}
