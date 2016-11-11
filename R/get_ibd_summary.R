#' IBD Segment Summary
#'
#' \code{getIBDsummary()} prints a brief summary of the detected IBD segments to the console.
#' @param ped.genotypes A list containing 2 objects.
#' See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.segments A data frame containing the IBD segments detected by isoRelate.
#' See the \code{Value} description in \code{link{getIBDsegments}} for more details on this input.
#' @export
getIBDsummary <- function(ped.genotypes, ibd.segments) {

  # check format of genotype data
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) stop ("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) stop ("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  if (!is.data.frame(pedigree)) stop ("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6) stop ("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", "moi", "aff")))
    stop ("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")

  # check ibd segment
  if (!is.data.frame(ibd.segments)) stop ("'ibd.segments' has incorrect format - must be a data.frame")
  if (ncol(ibd.segments) != 15 | any(colnames(ibd.segments) != c("fid1","iid1","fid2","iid2","chr","start_snp","end_snp","start_position_bp",
                                                             "end_position_bp", "start_position_M", "end_position_M", "number_snps", "length_bp",
                                                             "length_M", "ibd_status")))
    stop("'ibd.segments' has incorrect format - must have 15 columns")
  if (nrow(ibd.segments) == 0) stop ("no IBD segments detected")

  total.ibd.seg <- nrow(ibd.segments)
  ave.length.bp <- round(mean(ibd.segments[,"length_bp"]),0)

  num.isolates <- nrow(ped.genotypes[[1]])
  num.pairs <- num.isolates*(num.isolates-1)/2
  num.pairs.ibd <- length(unique(paste(ibd.segments[,"fid1"],ibd.segments[,"iid1"],
                                       ibd.segments[,"fid2"],ibd.segments[,"iid2"])))
  num.isolates.ibd <- length(unique(c(paste(ibd.segments[,"fid1"],ibd.segments[,"iid1"]),
                                      paste(ibd.segments[,"fid2"],ibd.segments[,"iid2"]))))

  cat(num.isolates.ibd,"/",num.isolates,"isolates IBD\n")
  cat(num.pairs.ibd,"/",num.pairs,"pairs of isolates IBD\n")
  cat("total IBD segments detected =",total.ibd.seg,"\n")
  cat("average length of segments (bp) =",ave.length.bp,"\n")

  # per chromosome
  if (length(unique(ibd.segments[,"chr"])) > 1) {
    for (chr in sort(unique(ibd.segments[,"chr"]))) {
      ibd.chr <- ibd.segments[ibd.segments[,"chr"] == chr,]
      cat(nrow(ibd.chr),"IBD segments detected on chromosome",chr,"\n")
    }
  }
}












