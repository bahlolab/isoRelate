#' Between Isolates Parameter Estimation
#' 
#' \code{getIBDparameters()} estimates the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD between all pairwise
#' combinations of isolates. 
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param number.cores Integer. The number of cores used for parallele execution.
#' @return A data frame with the following eight columns:
#' \enumerate{
#' \item Family 1 ID
#' \item Isolate 1 ID
#' \item Family 2 ID
#' \item Isolate 2 ID
#' \item The number of meiosis
#' \item Probability of sharing 0 alleles IBD
#' \item Probability of sharing 1 allele IBD
#' \item Probability of sharing 2 alleles IBD
#' }
#' where each row describes a unique pair of isolates. The data frame is headed \code{fid.1, iid1, fid.2, iid2, m, ibd0, ibd1} and \code{ibd2} respectively.
#' @importFrom foreach "%dopar%"
#' @export
getIBDparameters <- function(ped.genotypes, number.cores = 1){
  
  # check input ped and genotypes
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
  
  # check numeric input parameters
  stopifnot(is.numeric(number.cores))
  
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
    ibd.estimates.0 <- foreach::foreach(pair=pair.group, .combine='rbind') %dopar% {
      fid.1    <- as.character(isolate.pairs[pair,1])
      iid.1    <- as.character(isolate.pairs[pair,2])
      fid.2    <- as.character(isolate.pairs[pair,3])
      iid.2    <- as.character(isolate.pairs[pair,4])
      gender.1 <- pedigree[pedigree[,"fid"] == fid.1 & pedigree[,"iid"] == iid.1,"moi"]    
      gender.2 <- pedigree[pedigree[,"fid"] == fid.2 & pedigree[,"iid"] == iid.2,"moi"]     
      pair.genotypes   <- cbind(genotypes[,paste(fid.1,iid.1,sep="/")], genotypes[,paste(fid.2,iid.2,sep="/")])
      pop.allele.freqs <- genotypes[,"FREQ"]
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