#' Between Isolates IBD Segment Detection
#' 
#' \code{getIBDbetween()} detects genomic regions shared IBD between all pairwise combinations of isolates. 
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param parameters A data frame containing meioses and IBD probability estimates for all pairwise combinations of isolates.
#' See the \code{Value} description in \code{\link{between_parameter_estimates}} for more details on this input.
#' @param number.cores Integer. The number of cores used for parallele execution.
#' @param minimum.snps An integer value denoting the minimum number of SNPs in an IBD segment for it to be reported. The default value is 20 SNPs.
#' @param minimum.length.bp The minimum length of a reported IBD segment. The default value is 50,000 bp.
#' @param error The genotyping error rate. The default value is 0.001.
#' @return A data frame with the following columns
#'  \enumerate{
#' \item Family 1 ID
#' \item Isolate 1 ID
#' \item Family 2 ID
#' \item Isolate 2 ID
#' \item Chromosome
#' \item SNP identifier
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
#' @importFrom foreach "%dopar%"
#' @export
getIBDbetween <- function(ped.genotypes, parameters, number.cores = 1, minimum.snps = 20, minimum.length.bp = 50000, error = 0.001){
  
  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 2)
  pedigree  <- ped.genotypes[[1]] 
  genotypes <- ped.genotypes[[2]]
  
  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  
  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) < 8 & nrow(genotypes) <= 1)
    stop ("ped.genotypes has incorrect format")
  colnames(genotypes)[1:5] <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp", "FREQ")
  
  # check paramters file is a dataframe with correct fields
  stopifnot(is.data.frame(parameters))
  if (ncol(parameters) != 8) 
    stop ("parameters has incorrect format")
  colnames(parameters) <- c("fid1", "iid1", "fid2", "iid2", "m", "ibd0", "ibd1", "ibd2")
  
  # check numeric input parameters
  stopifnot(is.numeric(number.cores))
  stopifnot(is.numeric(minimum.snps))
  stopifnot(is.numeric(minimum.length.bp))
  stopifnot(is.numeric(error))
  
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
  chromosome <- unique(as.character(genotypes[,"CHROMOSOME"]))
  
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
    ibd.segments.0 <- foreach::foreach(pair=pair.group, .combine='rbind') %dopar% {
      
      # select pair
      fid.1    <- as.character(isolate.pairs[pair,1])
      iid.1    <- as.character(isolate.pairs[pair,2])
      fid.2    <- as.character(isolate.pairs[pair,3])
      iid.2    <- as.character(isolate.pairs[pair,4])
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
        
        pop.allele.freqs  <- genotypes[genotypes[,"CHROMOSOME"] == chrom,"FREQ"]
        pair.genotypes    <- cbind(genotypes[genotypes[,"CHROMOSOME"] == chrom,paste(fid.1,iid.1,sep="/")], genotypes[genotypes[,"CHROMOSOME"] == chrom,paste(fid.2,iid.2,sep="/")])
        positions.m  <- genotypes[genotypes[,"CHROMOSOME"] == chrom,"POSITION.M"]
        positions.bp <- genotypes[genotypes[,"CHROMOSOME"] == chrom,"POSITION.bp"]
        chromosomes  <- as.character(genotypes[genotypes[,"CHROMOSOME"] == chrom,"CHROMOSOME"])
        markers      <- as.character(genotypes[genotypes[,"CHROMOSOME"] == chrom,"MARKER"])
        number.snps  <- length(positions.m)
        
        # get viterbi IBD segments
        viterbi <- calculateViterbi(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2)

        # get posterior probabilities
        #gamma <- calculateGamma(number.states, initial.prob, meiosis, number.snps, pair.genotypes, pop.allele.freqs, positions.m, error, gender.1, gender.2)
        #if(number.states == 2) gamma <- cbind(gamma, rep(0,dim(gamma)[1]))
        #posterior.prob <- (gamma[,2]/2 + gamma[,3])
        #posterior.prob.sum[genotypes[,"CHROMOSOME"] == chrom] <- posterior.prob.sum[genotypes[,"CHROMOSOME"] == chrom] + posterior.prob
        
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
  for(i in 1:7)
    ibd.segments[,i] <- as.character(ibd.segments[,i])
  for(i in 8:15)
    ibd.segments[,i] <- as.numeric(as.character(ibd.segments[,i]))
  
  number.pairs.ibd <- length(unique(paste(ibd.segments[,1],ibd.segments[,2],ibd.segments[,3],ibd.segments[,4])))
  cat(paste(number.pairs.ibd,"pairs inferred IBD\n"))
  cat(paste(nrow(ibd.segments),"IBD segments detected\n"))
  
  return(ibd.segments)
}


