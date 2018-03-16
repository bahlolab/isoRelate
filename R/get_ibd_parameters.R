#' to prevent notes
globalVariables("pair.i")
#' Parameter Estimation
#'
#' \code{getIBDparameters()} estimates the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD between all pairwise
#' combinations of isolates.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' Note the family IDs and isolate IDs in object 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param number.cores Positive integer. The number of cores used for parallel execution.
#' @return A data frame with the following eight columns:
#' \enumerate{
#' \item Family 1 ID
#' \item Isolate 1 ID
#' \item Family 2 ID
#' \item Isolate 2 ID
#' \item The number of meiosis separating the pair
#' \item Probability of sharing 0 alleles IBD
#' \item Probability of sharing 1 allele IBD
#' \item Probability of sharing 2 alleles IBD
#' }
#' where each row describes a unique pair of isolates. The data frame is headed \code{fid1, iid1, fid2, iid2, m, ibd0, ibd1} and \code{ibd2} respectively.
#' @importFrom foreach "%dopar%"
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile
#' @export
#' @seealso \code{\link{getGenotypes}} and \code{\link{getIBDsegments}}.
#' @examples
#' # following processing and filtering of genotype data,
#' # we estimate the proportion of genome shared IBD
#' my_parameters <- getIBDparameters(ped.genotypes = png_genotypes,
#'                                   number.cores = 1)
#'
#' head(my_parameters)
getIBDparameters <- function(ped.genotypes, number.cores = 1){

  # check input ped and genotypes
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) stop ("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) stop ("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  if (!is.data.frame(pedigree)) stop ("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")
  if (!is.data.frame(genotypes)) stop ("'ped.genotypes' has incorrect format - 'genotypes' is not a data.frame")

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6) stop ("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", "moi", "aff")))
    stop ("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) <= 6) stop ("'ped.genotypes' has incorrect format - minimum of 2 isolates required for analysis")
  if (nrow(genotypes) <= 1) stop ("'ped.genotypes' has incorrect format - too few SNPs for analysis")
  if (any(colnames(genotypes)[1:5] != c("chr", "snp_id", "pos_M","pos_bp", "freq")))
    stop("'ped.genotypes' has incorrect format - 'genotypes' must have columns labelled: chr, snp_id, pos_M, pos_bp and freq")

  # check numeric input parameters
  if (!is.vector(number.cores)) stop ("'number.cores' has incorrect format - must be a vector")
  if (!is.numeric(number.cores)) stop ("'number.cores' has incorrect format - must be numeric")
  if (length(number.cores) != 1) stop ("'number.cores' has incorrect format - must be a single numeric value")
  if (number.cores < 1) stop ("'number.cores' has incorrect format - must be >= 1")

  # create new isolate IDs from PED FIDs and IIDs
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  number.pairs  <- 1:nrow(isolate.pairs)

  # calculate quantiles from the number of pairs, for progress bar only
  pair.quantiles   <- unique(round(quantile(number.pairs,probs=seq(0,0.9999,0.01))))
  number.quantiles <- length(pair.quantiles)

  # create progress bar
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)

  # define number of cores
  doParallel::registerDoParallel(cores=number.cores)

  # for each subgroup of pairs belonging to the quantiles, get parameters (for progress bar)
  start <- 1
  ibd.estimates <- NULL
  for (quantile.group in 1:number.quantiles) {

    # assign pairs to a subgroup based on which quantile they're in (for progress bar)
    if (number.quantiles == nrow(isolate.pairs))
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      pair.group <- start:(start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles))
      pair.group <- start:nrow(isolate.pairs)

    # get IBD parameters for subgroups of pairs
    ibd.estimates.0 <- foreach::foreach(pair.i=pair.group, .combine='rbind') %dopar% {
      fid.1    <- as.character(isolate.pairs[pair.i,1])
      iid.1    <- as.character(isolate.pairs[pair.i,2])
      fid.2    <- as.character(isolate.pairs[pair.i,3])
      iid.2    <- as.character(isolate.pairs[pair.i,4])
      gender.1 <- pedigree[pedigree[,"fid"] == fid.1 & pedigree[,"iid"] == iid.1,"moi"]
      gender.2 <- pedigree[pedigree[,"fid"] == fid.2 & pedigree[,"iid"] == iid.2,"moi"]
      pair.genotypes   <- cbind(genotypes[,paste(fid.1,iid.1,sep="/")], genotypes[,paste(fid.2,iid.2,sep="/")])
      pop.allele.freqs <- genotypes[,"freq"]
      ibd.estimates <- cbind(fid.1, iid.1, fid.2, iid.2, IBDparameters(pair.genotypes, pop.allele.freqs, gender.1, gender.2))
      colnames(ibd.estimates) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
      ibd.estimates
    }
    ibd.estimates <- rbind(ibd.estimates, ibd.estimates.0)

    # assign new start pair for next subgroup
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      start <- start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group]

    # update progress bar
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)

  # format parameter estimates
  ibd.estimates <- data.frame(ibd.estimates)
  for (i in 1:4)
    ibd.estimates[,i] <- as.character(ibd.estimates[,i])
  for (i in 5:8)
    ibd.estimates[,i] <- as.numeric(as.character(ibd.estimates[,i]))

  return(ibd.estimates)
}
