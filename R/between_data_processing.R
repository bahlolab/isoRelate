#' IsoRelate Pre-Analysis Data Processing 
#' 
#' \code{getGenotypes()} performs pre-analysis data processing of PLINK formatted unphase haplotype data, 
#' including removal of SNPs and isolates with high proportions of missing data and SNPs with low minor 
#' allele frequencies. It also calculates population allele frequencies for each SNP from either the 
#' input dataset or a specified reference dataet.
#' @param ped.map A list with 2 objects:
#' \enumerate{
#' \item An object of class pedMatrix which contains the PLINK PED information as described in \code{\link{between_load_pedmap}}
#' \item An object of class mapMatrix which contains the PLINK MAP information as described in \code{\link{between_load_pedmap}}
#' }
#' @param maf A numeric value denoting the smallest minor allele frequency allowed in the analysis. The 
#' default value is 0.01.
#' @param isolate.max.missing A numeric value denoting the maximum proportion of missing data allowed for 
#' each isolate. The default value is 0.1.
#' @param snp.max.missing A numeric value denoting the maximum proportion of missing data allowed for each 
#' SNP. The default value is 0.1.
#' @param reference.ped.map An optional list containing reference data used to calculate population allele 
#' frequencies. The list has 2 objects in the same format as the input 
#' \code{ped.map} objects. The default value is \code{NULL} which calculates the population allele 
#' frequencies from the input data.
#' @return A list of two objects: 
#' \enumerate{
#' \item A pedigree containing the isolates that remain after filtering. The pedigree is the first six columns 
#' of the PED file and these columns are
#' headed \code{fid, iid, pid, mid, moi} and \code{aff} respectively.
#' \item A data frame with the first five columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (Morgans, M) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' \item Population allele frequency (type \code{"integer"})
#' } 
#' where each row describes a single marker. These columns are headed \code{CHROMOSOME, MARKER, POSITION.M, POSITION.bp} and \code{FREQ} respectively. 
#' Columns 6 onwards contain the genotype data for each isolate, where a single column corresponds to a single isolate. These columns are 
#' labelled with merged family IDs and isolate IDs separated by a slash symbol (/). 
#' }
#' @export
getGenotypes <- function(ped.map, maf = 0.01, isolate.max.missing = 0.1, snp.max.missing = 0.1, reference.ped.map = NULL){

  # check the input parameters
  
  # check input PED and MAP files
  stopifnot(is.list(ped.map) | length(ped.map) != 2)
  input.ped <- ped.map[[1]]
  input.map <- ped.map[[2]]
  
  # check the PED and MAP files have the same number of SNPs
  n.snps.map <- 2*nrow(input.map)+6
  n.snps.ped <- ncol(input.ped)
  if(is.null(n.snps.ped) | is.null(n.snps.map)) stop()
  if (n.snps.ped != n.snps.map) {
    #stop ("PED and MAP files are not in the correct format")
  }

  # check the MAP file has 4 coloumns
  if (ncol(input.map) != 4)
    stop ("MAP file has incorrect format")
  colnames(input.map) <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp")
  
  # check numeric input parameters
  stopifnot(is.numeric(maf))
  stopifnot(is.numeric(isolate.max.missing))
  stopifnot(is.numeric(snp.max.missing))
  
  # check reference data
  if (!is.null(reference.ped.map)) {
    stopifnot(is.list(reference.ped.map) & length(reference.ped.map) == 2)
    reference.ped <- reference.ped.map[[1]]
    reference.map <- reference.ped.map[[2]]
    
    # check the PED and MAP files have the same number of SNPs
    if (ncol(reference.ped) != (2*nrow(reference.map)+6))
      stop ("reference PED and MAP files are not in the correct format")
    
    # check the MAP file has 4 coloumns
    if (ncol(reference.map) != 4)
      stop ("reference MAP file has incorrect format")
    colnames(reference.map) <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp")
  }
  
  
  # create new isolate IDs from PED FIDs and IIDs
  isolate.names <- paste(input.ped[,1], input.ped[,2], sep="/")
  
  
  # merge input data with reference data
  if (!is.null(reference.ped.map)) {
    input.map.v1      <- cbind(1:nrow(input.map), input.map)
    reference.map.v1  <- cbind(1:nrow(reference.map), reference.map)
    input.map.v1      <- merge(input.map.v1, reference.map.v1, by.x="MARKER", by.y="MARKER")
    input.map.v1      <- input.map.v1[order(input.map.v1[,"1:nrow(input.map)"]),]
    input.ped.columns <- c(1:6, 2*input.map.v1[,"1:nrow(input.map)"] + 5, 2*input.map.v1[,"1:nrow(input.map)"] + 6)
    input.ped.columns <- input.ped.columns[order(input.ped.columns)]
    input.ped.v1      <- input.ped[,input.ped.columns]
    reference.ped.columns  <- c(1:6, 2*input.map.v1[,"1:nrow(reference.map)"] + 5, 2*input.map.v1[,"1:nrow(reference.map)"] + 6)
    reference.ped.columns  <- reference.ped.columns[order(reference.ped.columns)]
    reference.ped.v1       <- reference.ped[,reference.ped.columns]
    input.map.v2           <- input.map.v1[,c("CHROMOSOME.x", "MARKER", "POSITION.M.x", "POSITION.bp.x")]
    colnames(input.map.v2) <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp")
  } else {
    input.map.v2 <- input.map
    input.ped.v1 <- input.ped
  }
  
  
  # call genotypes
  input.matrix        <- as.matrix(input.ped.v1[,7:ncol(input.ped.v1)])
  input.genders       <- input.ped.v1[,5]
  input.genotypes.v0  <- cbind(input.map.v2, haplotypeToGenotype(input.matrix, input.genders))
  if (!is.null(reference.ped.map)) {
    reference.matrix       <- as.matrix(reference.ped.v1[,7:ncol(reference.ped.v1)])
    reference.genders      <- reference.ped.v1[,5]
    reference.genotypes.v0 <- cbind(input.map.v2, haplotypeToGenotype(reference.matrix, reference.genders))
  }
  
  
  # calculate allele frequencies form reference data
  if (is.null(reference.ped.map)) {
    pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(input.genotypes.v0[,5:ncol(input.genotypes.v0)]), input.ped.v1[,5])
    input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,5:ncol(input.genotypes.v0)])
  } else {
    pop.allele.freq    <- calculatePopAlleleFreq(as.matrix(reference.genotypes.v0[,5:ncol(reference.genotypes.v0)]), reference.ped.v1[,5])
    input.genotypes.v1 <- cbind(input.genotypes.v0[,c(1:4)],pop.allele.freq,input.genotypes.v0[,c(5:ncol(input.genotypes.v0))])
  }
  colnames(input.genotypes.v1) <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp", "FREQ", isolate.names)
  cat(paste("Begin filtering of ",length(isolate.names)," isolates and ",nrow(input.genotypes.v1)," SNPs...\n",sep=""))
  
  
  # remove SNPs with low population MAF
  input.genotypes.v2 <- subset(input.genotypes.v1, pop.allele.freq <= (1-maf) & pop.allele.freq >= maf)
  cat(paste(nrow(input.genotypes.v2)," SNPs remain after MAF removal...\n",sep=""))
  
  
  # remove snps with high missingness
  snp.missingness    <- calculateMissingness(as.matrix(t(input.genotypes.v2[,6:ncol(input.genotypes.v2)]))) 
  input.genotypes.v3 <- input.genotypes.v2[snp.missingness <= snp.max.missing,]
  cat(paste(nrow(input.genotypes.v3)," SNPs remain after missingness removal...\n",sep=""))
  
  
  # remove samples with high missingness
  isolate.missingness <- round(calculateMissingness(as.matrix(input.genotypes.v3[,6:ncol(input.genotypes.v3)])),digits=3)
  if (length(isolate.names[isolate.missingness > isolate.max.missing]) > 0) {
    cat(paste("*********\n WARNING: Removing ",isolate.names[isolate.missingness > isolate.max.missing]," because genotype missingness > ",isolate.max.missing*100,"%\n*********\n",sep=""))
    sample.keep        <- input.ped.v1[isolate.missingness <= isolate.max.missing,1:6]
    input.genotypes.v4 <- input.genotypes.v3[,c(1:5, which(isolate.missingness <= isolate.max.missing) + 5)]
    if(nrow(sample.keep) < 1) stop(paste("All isolates removed with missingness > ",isolate.max.missing*100,"%. No isolates remaining.",sep=""))
  } else {
    sample.keep        <- input.ped.v1[,1:6]
    input.genotypes.v4 <- input.genotypes.v3
  }
  colnames(sample.keep) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  cat(paste(ncol(input.genotypes.v4)-5," isolates remain after missingness removal...\n",sep=""))
  
  
  return.genotypes <- list(sample.keep, input.genotypes.v4)
  return(return.genotypes)
}
