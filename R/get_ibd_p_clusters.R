#' Genome Cluster Networks
#'
#' \code{getIBDpclusters()} produces a network of clusters of isolates that share a minimum proportion of genome IBD.
#' Isolates that do not share a minimum proportion IBD are not included in the network or output.
#' The networks are created using R package \code{igraph}.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.segments A data frame containing the IBD segments detected by isoRelate.
#' See the \code{Value} description in \code{\link{getIBDsegments}} for more details on this input.
#' @param prop Numeric value between (0,1].
#' The minimum proportion of genome shared IBD between a pair of isolates in order for the pair to be included in the network.
#' For example, if \code{prop=1} then two isolates will be included if they are IBD over the entire genome, i.e. identical isolates,
#' whereas \code{prop=0.5} will include isolates that share at least 50\% of their genome IBD.
#' The default is \code{prop=1}.
#' @param hi.clust Logical. Whether to perform hierarchical clustering using the \code{fastgreedy.community} approach in the \code{igraph} package.
#' @return A list of three objects named \code{clusters}, \code{i.network} and \code{hi.clust}:
#' \enumerate{
#' \item A list where each object contains the names of isolates that form a disjoint cluster in the network.
#' If hierarchical clustering has been performed then the clusters may not be disjoint.
#' \item An \code{igraph} network used in the construction of network plots. See \url{ http://igraph.org/r/} doe more details.
#' \item Logical. Whether or not hierarchical clustering has been performed.
#' }
#' @importFrom igraph graph.data.frame E E<- V V<- layout_with_fr fastgreedy.community components
#' @export
#' @seealso \code{\link{getGenotypes}}, \code{\link{getIBDsegments}} and \code{\link{getIBDiclusters}}.
getIBDpclusters <- function(ped.genotypes, ibd.segments, prop=1, hi.clust = FALSE){

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
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", "moi", "aff")))
    stop ("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")

  # check ibd.segments file is a dataframe with correct fields
  if (!is.data.frame(ibd.segments)) stop ("'ibd.segments' has incorrect format - 'ibd.segments' is not a data.frame")
  if (ncol(ibd.segments) != 15) stop (paste("'ibd.segments' has incorrect format - 'ibd.segments' must have 15 columns: fid1, iid1, fid2, iid2, chr",
                                            "start_snp, end_snp, start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps",
                                            "length_bp, length_M, ibd_status"))
  if (nrow(ibd.segments) == 0) stop ("no IBD segments detected")
  if (any(colnames(ibd.segments) != c("fid1","iid1","fid2","iid2","chr","start_snp","end_snp","start_position_bp","end_position_bp",
                                      "start_position_M", "end_position_M", "number_snps", "length_bp", "length_M", "ibd_status")))
    stop (paste("'ibd.segments' has incorrect format - 'ibd.segments' must have columns: fid1, iid1, fid2, iid2, chr",
                "start_snp, end_snp, start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps",
                "length_bp, length_M, ibd_status"))

  # check the proportion of genome IBD
  if (!is.vector(prop)) stop ("'prop' has incorrect format - must be a numeric vector")
  if (!is.numeric(prop)) stop ("'prop' must be a numeric value between 0 and 1 (inclusive)")
  if (length(prop) != 1) stop ("'prop' must be a single numeric value between 0 and 1 (inclusive)")
  if (prop > 1 | prop <= 0) stop ("'prop' must be a numeric value between (0,1]")

  # check hierarchical clustering paramter
  if (!is.vector(hi.clust)) stop ("'hi.clust' has incorrect format - must be a logical vector")
  if (!is.logical(hi.clust)) stop ("'hi.clust' has incorrect format - must be a single logical value")
  if (length(hi.clust) != 1) stop ("'hi.clust' has incorrect format - must be a single logical value")

  # calculate the length of genome
  genome.length <- 0
  for(chr in unique(as.character(genotypes[,"chr"]))) {
    genotypes.0 <- genotypes[genotypes[,"chr"] == chr,]
    genome.length <- genome.length + max(genotypes.0[,"pos_bp"]) - min(genotypes.0[,"pos_bp"])
  }

  # find pairs haring large amounts of genome
  ibd.pairs <- paste(ibd.segments[,"fid1"],ibd.segments[,"iid1"],ibd.segments[,"fid2"],ibd.segments[,"iid2"],sep="/")
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  isolate.pairs <- paste(isolate.pairs[,1],isolate.pairs[,2],isolate.pairs[,3],isolate.pairs[,4],sep="/")
  highly.related <- NULL
  for(i in 1:length(unique(ibd.pairs))) {
    ibd.segments.0 <- ibd.segments[ibd.pairs == unique(ibd.pairs)[i],]
    if(sum(ibd.segments.0[,"length_bp"])/genome.length >= prop)
      highly.related <- c(highly.related, unique(ibd.pairs)[i])
  }
  ibd.interval <- ibd.segments[ibd.pairs %in% highly.related,]
  if (nrow(ibd.interval) == 0) stop ("no pairs share >= ",prop*100,"% of their genome IBD")

  # create a subset data frame
  ibd.melt <- data.frame(from = paste(ibd.interval[,"fid1"],ibd.interval[,"iid1"],sep="/"),
                         to = paste(ibd.interval[,"fid2"],ibd.interval[,"iid2"],sep="/"),
                         value = ibd.interval[,"length_bp"])
  ibd.melt <- ibd.melt[!duplicated(paste(ibd.melt[,"from"], ibd.melt[,"to"], sep="/")),]

  # get nodes and their original name (from pedigree)
  nodes <- unique(c(as.character(ibd.melt[,1]),as.character(ibd.melt[,2])))

  # create an igraph dataframe
  i.network <- igraph::graph.data.frame(ibd.melt, nodes, directed=F)

  # get clusters
  if (hi.clust) {
    # hierarchical clustering
    clusters.h <- igraph::fastgreedy.community(i.network)

    # create a list of clusters
    cluster.list <- list()
    for (clust in 1:length(clusters.h)) {
      cluster.list[[clust]] <- clusters.h[[clust]]
    }
  } else {
    # disjoint clustering
    clusters <- igraph::components(i.network)
    cluster.belong <- clusters[[1]]

    # create a list of each cluster
    cluster.list <- list()
    for (clust in 1:clusters[[3]]) {
      cluster.list[[clust]] <- names(cluster.belong[cluster.belong == clust])
    }
  }

  # print summary info
  cluster.list.1 <- cluster.list[order(sapply(cluster.list,length), decreasing=T)]
  clusterSummary(cluster.list.1)
  return.list <- list(cluster.list.1, i.network, hi.clust)
  names(return.list) <- c("clusters","i.network","hi.clust")
  return(return.list)
}



