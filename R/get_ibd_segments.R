#' IBD Segment Detection
#'
#' \code{getIBDsegments()} detects genomic regions shared IBD between all pairwise combinations of isolates.
#' @param ped.genotypes A list containing 2 objects.
#' See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param parameters A data frame containing meioses and IBD probability estimates for all pairwise combinations of isolates.
#' See the \code{Value} description in \code{\link{getIBDparameters}} for more details on this input.
#' @param number.cores Positive integer. The number of cores used for parallel execution.
#' @param minimum.snps An integer value denoting the minimum number of SNPs in an IBD segment for it to be reported.
#' The default value is 20 SNPs.
#' @param minimum.length.bp The minimum length of a reported IBD segment. The default value is 50,000 bp.
#' @param error The genotyping error rate. The default value is 0.001.
#' @return A data frame with the following columns
#'  \enumerate{
#' \item Family 1 ID
#' \item Isolate 1 ID
#' \item Family 2 ID
#' \item Isolate 2 ID
#' \item Chromosome
#' \item Start SNP
#' \item End SNP
#' \item Start position bp
#' \item End position bp
#' \item Start position M
#' \item End position M
#' \item Number of SNPs
#' \item Length bp
#' \item Length M
#' \item IBD status (1 = one allele shared IBD, 2 = two alleles shared IBD)
#' }
#' where each row describes a unique IBD segment. The data frame is headed \code{fid1, iid1, fid2, iid2, chr, start_snp, end_snp,
#' start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps, length_bp, length_M} and \code{ibd_status} respectively.
#' @importFrom foreach "%dopar%"
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats quantile
#' @export
#' @seealso \code{\link{getGenotypes}} and \code{\link{getIBDparameters}}.
#' @examples
#' \dontrun{
#' # prior to IBD detection, parameter estimates must be estimated.
#' # Assuming this has been done, IBD inference is performed
#' my_ibd <- getIBDsegments(ped.genotypes = png_genotypes,
#'                         parameters = png_parameters,
#'                         number.cores = 1,
#'                         minimum.snps = 20,
#'                         minimum.length.bp = 50000,
#'                         error = 0.001)
#'
#' head(my_ibd)
#' }
getIBDsegments <- function(ped.genotypes, parameters, number.cores = 1, minimum.snps = 20, minimum.length.bp = 50000, error = 0.001){

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

  # check minimum.snps
  if (!is.vector(minimum.snps)) stop ("'minimum.snps' has incorrect format - must be a vector")
  if (!is.numeric(minimum.snps)) stop ("'minimum.snps' has incorrect format - must be numeric")
  if (length(minimum.snps) != 1) stop ("'minimum.snps' has incorrect format - must be a single numeric value")
  if (minimum.snps < 0) stop ("'minimum.snps' has incorrect format - must be >= 0")

  # check minimum.length.bp
  if (!is.vector(minimum.length.bp)) stop ("'minimum.length.bp' has incorrect format - must be a vector")
  if (!is.numeric(minimum.length.bp)) stop ("'minimum.length.bp' has incorrect format - must be numeric")
  if (length(minimum.length.bp) != 1) stop ("'minimum.length.bp' has incorrect format - must be a single numeric value")
  if (minimum.length.bp < 0) stop ("'minimum.length.bp' has incorrect format - must be >= 0")

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

  # create a vector for posterior probabilities for updating
  #posterior.prob.sum <- rep(0,nrow(genotypes))

  # for each subgroup of pairs belonging to the quantiles, get IBD segments (for progress bar)
  start <- 1
  ibd.segments <- NULL
  for (quantile.group in 1:number.quantiles) {

    # assign pairs to a subgroup based on which quantile they're in (for progress bar)
    if (number.quantiles == nrow(isolate.pairs))
      pair.group <- quantile.group
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      pair.group <- start:(start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group] - 1)
    if (number.quantiles < nrow(isolate.pairs) & quantile.group == max(number.quantiles))
      pair.group <- start:nrow(isolate.pairs)

    # get IBD segments for subgroups of pairs
    #ibd.segments <- foreach::foreach(pair=1:nrow(isolate.pairs),  .combine='merge_lists') %dopar% {
    ibd.segments.0 <- foreach::foreach(pair.i=pair.group, .combine='rbind') %dopar% {

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

      # for each chromosome, perform IBD analysis
      ibd.table.2 <- NULL
      for (chrom in chromosome) {

        pop.allele.freqs  <- genotypes[genotypes[,"chr"] == chrom,"freq"]
        pair.genotypes    <- cbind(genotypes[genotypes[,"chr"] == chrom,paste(fid.1,iid.1,sep="/")], genotypes[genotypes[,"chr"] == chrom,paste(fid.2,iid.2,sep="/")])
        positions.m  <- genotypes[genotypes[,"chr"] == chrom,"pos_M"]
        positions.bp <- genotypes[genotypes[,"chr"] == chrom,"pos_bp"]
        chromosomes  <- as.character(genotypes[genotypes[,"chr"] == chrom,"chr"])
        markers      <- as.character(genotypes[genotypes[,"chr"] == chrom,"snp_id"])
        number.snps  <- length(positions.m)

        # get viterbi IBD segments
        viterbi <- calculateViterbi(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2)

        # get posterior probabilities
        #gamma <- calculateGamma(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2)
        #if(number.states == 2) gamma <- cbind(gamma, rep(0,dim(gamma)[1]))
        #posterior.prob <- (gamma[,2]/2 + gamma[,3])
        #posterior.prob.sum[genotypes[,"chr"] == chrom] <- posterior.prob.sum[genotypes[,"chr"] == chrom] + posterior.prob

        # get IBD summary table
        ibd.results <- cbind(fid.1, iid.1, fid.2, iid.2, 1:number.snps, chromosomes, markers, positions.m, positions.bp, viterbi)
        colnames(ibd.results) <- c("fid1","iid1","fid2","iid2","markerNo","chr","marker","pos.m","pos.bp","viterbi")
        ibd.table.1 <- IBDTable(ibd.results)

        # remove IBD segments with less than minimum.snps and less than minimum.length.bp
        # FIXME! can take a while to rbind
        if(length(ibd.table.1) != 0){
          ibd.table.2 <- rbind(ibd.table.2, ibd.table.1[as.numeric(ibd.table.1[,"number.snps"]) >= minimum.snps & as.numeric(ibd.table.1[,"length.bp"]) >= minimum.length.bp,])
        }
      }

      #list(ibd.table.2, posterior.prob.sum)
      ibd.table.2
    }
    ibd.segments <- rbind(ibd.segments, ibd.segments.0)

    # assign new start pair for next subgroup
    if (number.quantiles < nrow(isolate.pairs) & quantile.group != max(number.quantiles))
      start <- start + pair.quantiles[quantile.group+1] - pair.quantiles[quantile.group]

    # update progress bar
    setTxtProgressBar(pb, quantile.group)
  }
  close(pb)

  # format IBD segments
  #ibd.segments[[2]] <- ibd.segments[[2]]/nrow(isolate.pairs)
  rownames(ibd.segments) <- NULL
  ibd.segments <- data.frame(ibd.segments)
  if (nrow(ibd.segments) > 0) {
    for(i in 1:7)
      ibd.segments[,i] <- as.character(ibd.segments[,i])
    for(i in 8:15)
      ibd.segments[,i] <- as.numeric(as.character(ibd.segments[,i]))
    colnames(ibd.segments) <- c("fid1","iid1","fid2","iid2","chr","start_snp","end_snp","start_position_bp",
                                "end_position_bp", "start_position_M", "end_position_M", "number_snps", "length_bp",
                                "length_M", "ibd_status")
    number.pairs.ibd <- length(unique(paste(ibd.segments[,1],ibd.segments[,2],ibd.segments[,3],ibd.segments[,4])))
    cat(paste(number.pairs.ibd,"pairs inferred IBD\n"))
    cat(paste(nrow(ibd.segments),"IBD segments detected\n"))
  } else {
    cat("0 pairs inferred IBD\n")
    cat("0 IBD segments detected\n")
  }

  return(ibd.segments)
}


