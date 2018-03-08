#' IBD Posterior Probabilities
#'
#' \code{getIBDposterior()} calculates the posterior probabilities of IBD sharing between pairs of isolates.
#' @param ped.genotypes A list containing 2 objects.
#' See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param parameters A data frame containing meioses and IBD probability estimates for all pairwise combinations of isolates.
#' See the \code{Value} description in \code{\link{getIBDparameters}} for more details on this input.
#' @param number.cores Positive integer. The number of cores used for parallel execution.
#' @param error The genotyping error rate. The default value is 0.001.
#' @return A data frame with the first four columns:
#'  \enumerate{
#' \item Chromosome
#' \item SNP identifiers
#' \item Genetic map distance
#' \item Base-pair position
#' }
#' where each row describes a single SNP. These columns are headed \code{chr, snp_id, pos_M} and \code{pos_bp} respectively.
#' Columns 5 onwards contain the posterior probabilities for each pair of isolates, where a single column corresponds to one pair of isolates.
#' These columns are labeled with merged family IDs and isolate IDs separated by a slash symbol (/).
#' @importFrom foreach "%dopar%"
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile
#' @export
getIBDposterior <- function(ped.genotypes, parameters, number.cores = 1, error = 0.001){

  # check format of input data
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

  # check paramters file is a dataframe with correct fields
  if (!is.data.frame(parameters)) stop ("'parameters' has incorrect format - must be a data.frame")
  if (ncol(parameters) != 8) stop ("'parameters' has incorrect format - must have 8 columns labelled: fid1, iid1, fid2, iid2, m, ibd0, ibd1 and ibd2")
  if (any(colnames(parameters) != c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")))
    stop ("'parameters' has incorrect format - must have 8 columns labelled: fid1, iid1, fid2, iid2, m, ibd0, ibd1 and ibd2")

  # check number.cores
  if (!is.vector(number.cores)) stop ("'number.cores' has incorrect format - must be a vector")
  if (!is.numeric(number.cores)) stop ("'number.cores' has incorrect format - must be numeric")
  if (length(number.cores) != 1) stop ("'number.cores' has incorrect format - must be a single numeric value")
  if (number.cores < 1) stop ("'number.cores' has incorrect format - must be >= 1")

  # check error
  if (!is.vector(error)) stop ("'error' has incorrect format - must be a vector")
  if (!is.numeric(error)) stop ("'error' has incorrect format - must be numeric")
  if (length(error) != 1) stop ("'error' has incorrect format - must be a single numeric value")
  if (error < 0 | error > 1) stop ("'error' has incorrect format - must be between 0 and 1 (inclusive)")

  # create new isolate IDs from PED FIDs and IIDs
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  number.pairs  <- 1:nrow(isolate.pairs)

  # calculate quantiles from the number of pairs, for progress bar only
  pair.quantiles <- unique(round(quantile(number.pairs,probs=seq(0,0.9999,0.01))))
  number.quantiles   <- length(pair.quantiles);

  # create progress bar
  pb <- txtProgressBar(min = 0, max = number.quantiles, style = 3)

  # define number of cores
  doParallel::registerDoParallel(cores=number.cores)

  # create a character vector of unique chromosomes
  chromosome <- unique(as.character(genotypes[,"chr"]))

  # for each subgroup of pairs belonging to the quantiles, get IBD segments (for progress bar)
  start <- 1
  posterior.prob.all <- genotypes[,c("chr", "snp_id", "pos_M","pos_bp")]
  for (quantile.group in 1:number.quantiles) {

    # assign pairs to a subgroup based on which quantile they're in (for progress bar)
    if (number.quantiles == nrow(isolate.pairs))
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      pair.group <- start:(start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles))
      pair.group <- start:nrow(isolate.pairs)

    # get IBD segments for subgroups of pairs
    posterior.prob.0 <- foreach::foreach(pair.i=pair.group, .combine='cbind') %dopar% {

      # select pair
      fid.1    <- as.character(isolate.pairs[pair.i,1])
      iid.1    <- as.character(isolate.pairs[pair.i,2])
      fid.2    <- as.character(isolate.pairs[pair.i,3])
      iid.2    <- as.character(isolate.pairs[pair.i,4])
      gender.1 <- pedigree[pedigree[,"fid"] == fid.1 & pedigree[,"iid"] == iid.1,"moi"]
      gender.2 <- pedigree[pedigree[,"fid"] == fid.2 & pedigree[,"iid"] == iid.2,"moi"]

      # define number of states in model based on genders (MOI)
      if (gender.1 == 2 & gender.2 == 2) {
        number.states = 3
      } else
        number.states = 2

      # get model paramters
      meiosis      <- as.numeric(parameters[parameters[,"fid1"] == fid.1 & parameters[,"fid2"] == fid.2 & parameters[,"iid1"] == iid.1 & parameters[,"iid2"] == iid.2,"m"])
      initial.prob <- as.numeric(parameters[parameters[,"fid1"] == fid.1 & parameters[,"fid2"] == fid.2 & parameters[,"iid1"] == iid.1 & parameters[,"iid2"] == iid.2,c("ibd0","ibd1","ibd2")])

      # change initial probabilities to allow switiching between states.. re-think?? FIXME!!
      if (initial.prob[1] == 1) {
        initial.prob[1] <- 0.999
        initial.prob[2] <- 0.001
      }
      if (initial.prob[2] == 1) {
        initial.prob[1] <- 0.001
        initial.prob[2] <- 0.999
      }

      # create a vector for posterior probabilities for updating
      posterior.prob.sum <- rep(0,nrow(genotypes))

      # for each chromosome, perform IBD analysis
      for (chrom in chromosome) {

        pop.allele.freqs  <- genotypes[genotypes[,"chr"] == chrom,"freq"]
        pair.genotypes    <- cbind(genotypes[genotypes[,"chr"] == chrom,paste(fid.1,iid.1,sep="/")], genotypes[genotypes[,"chr"] == chrom,paste(fid.2,iid.2,sep="/")])
        positions.m  <- genotypes[genotypes[,"chr"] == chrom,"pos_M"]
        positions.bp <- genotypes[genotypes[,"chr"] == chrom,"pos_bp"]
        chromosomes  <- as.character(genotypes[genotypes[,"chr"] == chrom,"chr"])
        markers      <- as.character(genotypes[genotypes[,"chr"] == chrom,"snp_id"])
        number.snps  <- length(positions.m)

        # get posterior probabilities
        gamma <- calculateGamma(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2)
        if(number.states == 2) gamma <- cbind(gamma, rep(0,dim(gamma)[1]))
        posterior.prob <- (gamma[,2]/2 + gamma[,3])
        posterior.prob.sum[genotypes[,"chr"] == chrom] <- posterior.prob.sum[genotypes[,"chr"] == chrom] + posterior.prob

      }
      pair_id <- paste(fid.1, iid.1, fid.2, iid.2, sep="/")
      posterior.prob.sum <- data.frame(posterior.prob.sum)
      colnames(posterior.prob.sum) <- pair_id
      posterior.prob.sum
    }
    posterior.prob.all <- cbind.data.frame(posterior.prob.all, posterior.prob.0)

    # assign new start pair for next subgroup
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      start <- start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group]

    # update progress bar
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)

  return(posterior.prob.all)
}
