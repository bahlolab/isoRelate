#' Between Isolates IBD Plotting
#' 
#' \code{plotIBDpairwise()} produces figures of inferred IBD segments. Regions of IBD are represented by coloured rectangles, where 
#' \code{IBD=1} (one allele shared IBD) and \code{IBD=2} (two alleles shared IBD) regions are distingued by different colours. 
#' The user can choose which chromosomes they wish to plot as well as zoom in on specific regions of interest.
#' @param ibd.segments A data frame with inferred IBD tracts. See the \code{Value} description in \code{\link{between_ibd_segments}} 
#' for more details. This data frame can be alterted prior to plotting, for example subsetting by specific combinations of isolate pairs.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_processing}} for more details on this input.
#' @param chromosomes A vector of chromosomes to plot. Can be either numeric or character.
#' @param zoom A data frame containing three columns:
#' \enumerate{
#' \item Chromosome
#' \item Start location bp
#' \item End location bp
#' }
#' More than one region can be specified per chromosome however they will be plotted in spearate figures.
#' @param ibd.colours A vector with two values. These can be either numeric or the name of colours to be used in the IBD plot (must be vaild colour
#' names).
#' @param minimum.length.bp The minimum length of an IBD segment in bp to be included in the summary figure. The default value is 0 bp.
#' @param pairs.per.page The maximum number of pairwise results to plot on each page of the PDF.
#' @param iid.only Use only the isolate IDs on the figure legend of the plot. If using, the isolate IDs must be unique for each isolate.
#' \code{iid.only} and \code{fid.only} are mutually exclusive.
#' @param fid.only Use only the family IDs on the figure legend of the plot. If using, the family IDs must be unique for each isolate.
#' @export

plotIBDpairwise <- function(ped.genotypes, ibd.segments, chromosomes = NULL, zoom = NULL, ibd.colours = c("dodgerblue","darkgoldenrod1"),
                                        minimum.length.bp = 0, pairs.per.page = 50, iid.only = FALSE, fid.only = FALSE){
  
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
  
  # check chromosomes
  if (!is.null(chromosomes))
    stopifnot(is.vector(chromosomes))
  
  # check zoom file is a dataframe (if not NULL) with 3 columns
  if (!is.null(zoom)){
    stopifnot(is.data.frame(zoom))
    stopifnot(ncol(zoom) != 3) 
  }
  
  # check chromosome colours
  stopifnot(is.vector(ibd.colours))
  if (length(ibd.colours) > 2){
    warning("ibd.colours contains more than 2 colours. Using first 2 colours only.")
    ibd.colours <- ibd.colours[1:2]
  }
    
  # check numeric input parameters
  stopifnot(is.numeric(minimum.length.bp))
  stopifnot(is.numeric(pairs.per.page))
  
  # check logical input parameters
  stopifnot(is.logical(fid.only))
  stopifnot(is.logical(iid.only))
  
  # minimum IBD tract length to plot
  ibd.segments.1 <- ibd.segments[ibd.segments[,"length.bp"] >= minimum.length.bp,]
  if (nrow(ibd.segments.1) == 0) 
    stop(paste("No IBD tracts with length greater than ",minimum.length.bp," bp",sep=""))
  
  # plot figure
  if (is.null(zoom)) {
    plotIBD(ibd.segments.1, ped.genotypes, chromosomes, ibd.colours, pairs.per.page, iid.only, fid.only)
  } else {
    plotIBDzoom(ibd.segments.1, ped.genotypes, ibd.colours, pairs.per.page, iid.only, fid.only, zoom)
  }
}