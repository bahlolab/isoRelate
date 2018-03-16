#' Binary IBD Matrix
#'
#' \code{getIBDmatrix()} produces a binary matrix of IBD (1) and non-IBD (0) results for each SNP and isolate pair combination.
#' Each row identifies a unique SNP while each column identifies a unique isolate pair.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.segments A data frame containing the IBD segments detected by isoRelate.
#' See the \code{Value} description in \code{\link{getIBDsegments}} for more details on this input.
#' @return A data frame with the first four columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' }
#' where each row describes a unique SNP. Columns 1-4 are headed \code{chr, snp_id, pos_M} and \code{pos_bp} respectively.
#' Columns 5 onwards contain the binary IBD information for each isolate pair, where a single column corresponds to a single pair.
#' These columns are labeled with merged family IDs and isolate IDs separated by a slash symbol (/). For example fid1/iid1/fid2/iid2.
#' @export
#' @seealso \code{\link{getGenotypes}}, \code{\link{getIBDsegments}}, \code{\link{getIBDproportion}}, \code{\link{getIBDiR}}.
#' @examples
#' # generate a binary IBD matrix
#' my_matrix <- getIBDmatrix(ped.genotypes = png_genotypes,
#'                           ibd.segments = png_ibd)
getIBDmatrix <- function(ped.genotypes, ibd.segments){

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

  # check ibd.segments file is a dataframe with correct fields
  if (!is.data.frame(ibd.segments)) stop ("'ibd.segments' has incorrect format - 'ibd.segments' is not a data.frame")
  if (ncol(ibd.segments) != 15) stop (paste("'ibd.segments' has incorrect format - 'ibd.segments' must have 15 columns: fid1, iid1, fid2, iid2, chr",
        "start_snp, end_snp, start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps",
        "length_bp, length_M, ibd_status"))
  if (any(colnames(ibd.segments) != c("fid1","iid1","fid2","iid2","chr","start_snp","end_snp","start_position_bp","end_position_bp",
                          "start_position_M", "end_position_M", "number_snps", "length_bp", "length_M", "ibd_status")))
    stop (paste("'ibd.segments' has incorrect format - 'ibd.segments' must have columns: fid1, iid1, fid2, iid2, chr",
                "start_snp, end_snp, start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps",
                "length_bp, length_M, ibd_status"))


  # create a data frame of pairs; one unique pair per row; columns are fid1, iid1, fid2, iid2 and
  # a unique numeric pair identifier
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  isolate.pairs <- cbind(isolate.pairs, 1:nrow(isolate.pairs))
  colnames(isolate.pairs) <- c("fid1","iid1","fid2","iid2","id")

  # merge inferred IBD pairs with their numeric ID
  isolate.id <- merge(isolate.pairs, ibd.segments)

  # generate the binary IBD matrix
  ibd.binary.matrix <- IBDMatrix(as.character(genotypes[,"chr"]), genotypes[,"pos_bp"], nrow(isolate.pairs),
                              as.numeric(as.character(isolate.id[,"id"])), isolate.id[,"chr"], isolate.id[,"start_position_bp"],
                              isolate.id[,"end_position_bp"])
  ibd.locus.matrix <- data.frame(genotypes[,c("chr","snp_id","pos_M","pos_bp")],ibd.binary.matrix)
  colnames(ibd.locus.matrix) <- c("chr","snp_id","pos_M","pos_bp",
                                  paste(isolate.pairs[,"fid1"],isolate.pairs[,"iid1"],
                                        isolate.pairs[,"fid2"],isolate.pairs[,"iid2"], sep="/"))
  rownames(ibd.locus.matrix) <- NULL

  return(data.frame(ibd.locus.matrix))
}


