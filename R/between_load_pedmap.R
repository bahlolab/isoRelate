#' IsoRelate Load in PLINK PED and MAP files
#' 
#' \code{loadPEDMAP()} reads in a PLINK PED and PLINK MAP file using \code{fread} from the \code{data.table} package.
#' @param ped A character string containing the path and filename of the input PLINK PED file. The PED file is 
#' a white-space (space or tab) delimited file
#' with the first six columns:
#' \enumerate{
#' \item Family ID
#' \item Isolate ID
#' \item Paternal ID
#' \item Maternal ID
#' \item Multiplicity of infection (MOI) (1 = single infection or haploid, 2 = multiple infections or diploid)
#' \item Phenotype (1 = unaffected, 2 = affected, 0 = unknown)
#' }
#' The IDs are alphanumeric: the combination of family and individual ID should uniquely identify an isolate. The phenotype coulmn is not used
#' in IsoRelate analyses however is required for completeness of a standard pedigree. Genotypes (column 7 onwards) should
#' also be white-spaced delimited with the A and B alleles coded as 1 and 2 respectively and missing data coded as 0. All SNPs (whether haploid or not)
#' must have two alleles specified. For haploid chromosomes, genotypes should be specified as homozygous. Either both alleles should be missing (i.e. 0)
#' or neither. No header row should be given.
#' @param map A character string containing the path and filename of the input PLINK MAP file. A MAP file contains exactly four coloums of information:
#' \enumerate{
#' \item Chromosome
#' \item SNP identifier
#' \item Genetic map distance (centi morgans, cM)
#' \item Base-pair position
#' }
#' where each row describes a single marker. Genetic map distance and base-pair positions are expected to be positive values. The MAP file must
#' be ordered by increasing chromosomes and positions. SNP identifiers can contain any characters expect spaces or tabs; also you should avoid
#' * symbols in the names. The MAP file must contain as many markers as are in the PED file. No header row should be given.
#' @return A list of two objects: 
#' \enumerate{
#' \item An object of class \code{pedMatrix} containing the data from \code{ped}.
#' \item An object of class \code{mapMatrix} containing the data from \code{map}. 
#' }
#' @export

loadPEDMAP <- function(ped, map){
  # checks
  stopifnot(is.character(ped) & is.character(map))
  stopifnot(file.exists(ped) & file.exists(map))
  
  # read in ped and map
  input.ped <- data.table::fread(ped, data.table=FALSE)
  input.map <- data.table::fread(map, data.table=FALSE)
  
  # check dimensions of ped and map files
  if (ncol(input.map) != 4) {
    stop(paste(map, "is not in the correct format; requires 4 columns"))
  }
  if (ncol(input.ped) %% 2 != 0) {
    stop(paste(ped, "does not have an even number of columns"))
  }
  if (nrow(input.map) != (ncol(input.ped) - 6)/2) {
    stop(paste("ped has", (ncol(input.ped) - 6)/2, "SNPs and map has", nrow(input.map), "SNPs"))
  }
  if (ncol(input.ped) * nrow(input.ped) > 2^31-1) {
    stop("ped file exceeds R dimension limits..can not perform IBD analysis")
  }
  if (any(input.ped[,7:ncol(input.ped)] < 0) | any(input.ped[,7:ncol(input.ped)] > 2) | any(is.na(input.ped[,7:ncol(input.ped)]))) {
    stop("ped genotypes must be encoded as 0, 1 and 2")
  }
  if (any(is.na(input.map[,3])) | any(is.na(input.map[,4]))) {
    stop("genetic map positions NA for some SNPs")
  }

  # define class and dimension names
  colnames(input.ped) <- NULL
  rownames(input.ped) <- NULL
  colnames(input.map) <- NULL
  rownames(input.map) <- NULL
  
  # define map structure
  input.map[,1] <- as.character(input.map[,1])
  input.map[,2] <- as.character(input.map[,2])
  input.map[,3] <- as.numeric(input.map[,3])
  input.map[,4] <- as.numeric(input.map[,4])

  return(list(input.ped,input.map))
}