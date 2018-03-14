#' Interval Cluster Networks
#'
#' \code{getIBDiclusters()} produces a network of clusters of isolates that have been inferred IBD over a specified interval.
#' Isolates that are not IBD over the interval are not included in the network or output. The networks are created
#' using the R package \code{igraph}.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.segments A data frame containing the IBD segments detected by isoRelate.
#' See the \code{Value} description in \code{\link{getIBDsegments}} for more details on this input.
#' @param interval A vector of length 3 containing the region to identify clusters over. This vector should contain the
#' chromosome ID, the start of the interval in base-pairs and the end of the interval in base-pairs; in this order respectively.
#' @param prop Numeric value between 0 and 1 (inclusive).
#' The minimum proportion of an interval (in base-pairs) shared IBD between a pair of isolates in order for the pair to be included in the network.
#' For example, if \code{prop=1} then two isolates will be included if they are IBD over the entire interval,
#' whereas \code{prop=0.5} will include isolates that are IBD over at least 50\% of the interval.
#' The default is \code{prop=0}, which includes isolates if they have an IBD segment that overlaps the interval in any way.
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
#' @seealso \code{\link{getGenotypes}}, \code{\link{getIBDsegments}} and \code{\link{getIBDpclusters}}.
getIBDiclusters <- function(ped.genotypes, ibd.segments, interval = NULL, prop=0, hi.clust = FALSE){

  # check format of input data
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) stop ("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) stop ("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  if (!is.data.frame(pedigree)) stop ("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6) stop ("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
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

  # check if region specified
  if (is.null(interval)) stop ("must specify 'interval'")
  if (!is.vector(interval)) stop ("'interval' has incorrect format - must be a vector of length 3")
  if (length(interval) != 3) stop ("'interval' has incorrect format - must be a vector of length 3")
  interval.chr   <- as.character(interval[1])
  interval.start <- as.numeric(interval[2])
  interval.stop  <- as.numeric(interval[3])
  if (!(interval.chr %in% ibd.segments[,"chr"])) stop (paste("no pairs IBD on chromosome =",interval.chr))
  if (interval.start > interval.stop) stop(paste("interval start=",interval.start," is greater than interval end=",interval.stop,sep=""))

  # check the proportion of interval IBD over
  if (!is.vector(prop)) stop ("'prop' has incorrect format - must be a numeric vector")
  if (!is.numeric(prop)) stop ("'prop' must be a numeric value between 0 and 1 (inclusive)")
  if (length(prop) != 1) stop ("'prop' must be a single numeric value between 0 and 1 (inclusive)")
  if (prop > 1 | prop < 0) stop ("'prop' must be a numeric value between 0 and 1 (inclusive)")

  # check hierarchical clustering paramter
  if (!is.vector(hi.clust)) stop ("'hi.clust' has incorrect format - must be a logical vector")
  if (!is.logical(hi.clust)) stop ("'hi.clust' has incorrect format - must be a single logical value")
  if (length(hi.clust) != 1) stop ("'hi.clust' has incorrect format - must be a single logical value")

  # for each IBD segment, find segments overlapping interval
  ibd.chr <- ibd.segments[ibd.segments[,"chr"] == interval.chr,]
  ibd.overlap <- rep(0, nrow(ibd.chr))
  for(i in 1:nrow(ibd.chr)){
    my.overlap <- getOverlap(ibd.chr[i,c("start_position_bp","end_position_bp")], c(interval.start, interval.stop))
    if(my.overlap[2] > my.overlap[1])
      if((my.overlap[2] - my.overlap[1]) >= (as.numeric(interval.stop) - as.numeric(interval.start))*prop)
        ibd.overlap[i] <- 1
  }
  ibd.interval <- ibd.chr[ibd.overlap == 1,]
  if (nrow(ibd.interval) == 0) stop ("no IBD segments detected over interval")

  # remove duplicate isolates (from IBD=1 and IBD=2 inferred over same interval)
  ibd.pairs <- paste(ibd.interval[,"fid1"],ibd.interval[,"iid1"],ibd.interval[,"fid2"],ibd.interval[,"iid2"],sep="/")
  ibd.interval <- ibd.interval[!duplicated(ibd.pairs),]

  # create a subset data frame
  ibd.melt <- data.frame(from = paste(ibd.interval[,"fid1"],ibd.interval[,"iid1"],sep="/"),
                         to = paste(ibd.interval[,"fid2"],ibd.interval[,"iid2"],sep="/"),
                         value = ibd.interval[,"length_bp"])

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


#' A function that prints summary information of clusters
#' @param cluster.list a list of isolate groups where unique elements correspond to isolates in unique clusters
clusterSummary <- function(cluster.list) {
  cat(paste("Number of clusters = ",length(cluster.list),"\n",sep=""))
  cat(paste("Maximum cluster length = ",length(cluster.list[[1]]),"\n",sep=""))
  cat(paste("Minimum cluster length = ",length(cluster.list[[length(cluster.list)]]),"\n",sep=""))
}

