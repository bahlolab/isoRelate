#' Between Isolates Binary IBD Matrix
#' 
#' \code{getLocusMatrix()} produces a binary matrix of IBD (1) and non-IBD (0) for each SNP and pair combination.
#' The number of rows is equal to the number of SNPs and the number of columns is equal to the number of pairs.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @return A data frame the first four columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' }
#' where each row describes a unique SNP. These columns are headed \code{CHROMOSOME, MARKER, POSITION.cM} and \code{POSITION.bp} respectively.
#' Columns 5 onwards contain the binary IBD information for each isolate pair, where a single column corresponds to a single pair. 
#' These columns are labelled with merged family IDs and isolate IDs separated by a slash symbol (/). For example famA/ind1/famA/ind2. 
#' @export
getLocusMatrix <- function(ped.genotypes, ibd.segments){
  
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
  
  # check ibd.segments file is a dataframe with correct fields
  stopifnot(is.data.frame(ibd.segments))
  if (ncol(ibd.segments) != 15) 
    stop ("ibd.segments has incorrect format")
  colnames(ibd.segments) <- c("fid1", "ind1", "fid2", "ind2", "chr", "start.snp", "end.snp", "start.position.bp", 
                              "end.position.bp", "start.position.M", "end.position.M", "number.snps", "length.bp",
                              "length.M", "ibd.status")
  
  # create a data frame of pairs; one unique pair per row; columns are fid1, iid1, fid2, iid2 and
  # a unique numeric pair identifier
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  isolate.pairs <- cbind(isolate.pairs, 1:nrow(isolate.pairs))
  colnames(isolate.pairs) <- c("fid1","ind1","fid2","ind2","id")
  
  # merge inferred IBD pairs with their numeric ID
  isolate.id <- merge(isolate.pairs, ibd.segments)
  
  # generate the binary IBD matrix
  ibd.binary.matrix <- IBDMatrix(as.character(genotypes[,"CHROMOSOME"]), genotypes[,"POSITION.bp"], nrow(isolate.pairs),
                              as.numeric(as.character(isolate.id[,"id"])), isolate.id[,"chr"], isolate.id[,"start.position.bp"], 
                              isolate.id[,"end.position.bp"])
  ibd.locus.matrix <- data.frame(genotypes[,c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp")],ibd.binary.matrix)
  colnames(ibd.locus.matrix) <- c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp",
                                  paste(isolate.pairs[,"fid1"],isolate.pairs[,"ind1"],
                                        isolate.pairs[,"fid2"],isolate.pairs[,"ind2"], sep="/"))

  return(data.frame(ibd.locus.matrix))
}


