#' Between Isolates Clusters
#' 
#' \code{getClusters()} identifies clusters of isolates that share IBD over part of the entire specified region.
#' Each cluster represents a disjoint network where each isolate in the network is conneceted to all other isolates
#' in the network, either directly or through a series of connected isoaltes.
#' @param ped.genotypes a list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param locus.matrix a data frame containing the binary IBD information for each SNP and each pair. 
#' See the returned \code{value} in \code{\link{between_locus_matrix}} for more details.
#' @param snp a character string of length 1 giving the name of a SNP to find clusters of isolates who are IBD here.
#' Note that \code{snpName} and \code{interval} are mutually exclusive and only one option should be specified.
#' @param interval a vector of length 3 containing the region to identify clusters over. This vector should contain the
#' chromosome ID, the start of the interval in bp and the end of the interval in bp; in this order respectively. Note that 
#' \code{snpName} and \code{interval} are mutually exclusive and only one option should be specified.
#' @param prop numeric value between 0 and 1. The minimum proportion of an interval (in bp) shared IBD between a pair of isolates in order for the pair 
#' to be directly connected in a cluster. For example, if \code{prop=1} then two isolates will only be connected if they are IBD over the
#' entire specified region, whereas if \code{prop=0.5} then two isolates will be connected if they are IBD over at least half of the interval.
#' The default is \code{prop=0}, which connects two isolates if they are IBD at all over the interval.
#' @return A list with two objects:
#' \enumerate{
#' \item a list of disjoint clusters that contain the names of isoaltes that are part of the same network.
#' \item the binary IBD matrix for the specified region, where 1 represents IBD sharing by pairs over the 
#' specified interval and 0 otherwise.
#' }
#' @importFrom igraph graph.data.frame E E<- V V<- layout_with_fr fastgreedy.community components
#' @export
getClusters <- function(ped.genotypes, locus.matrix, snp=NULL, interval=NULL, prop=0, hclust=FALSE){
  
  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 2)
  pedigree <- ped.genotypes[[1]] 
  
  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  
  # check locus matrix input
  if (ncol(locus.matrix) < 5)
    stop ("locus.matrix has incorrect format")
  colnames(locus.matrix)[1:4] <- c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp")
  
  # check if region specified
  if (is.null(snp) & is.null(interval))
    stop("must specify 'snp' or 'interval'")
  
  # if interval specified
  if (!is.null(interval) & length(interval) == 3) {
    chr   <- as.character(interval[1])
    start <- as.numeric(interval[2])
    stop  <- as.numeric(interval[3])
    stopifnot(chr %in% locus.matrix[,"CHROMOSOME"])
    if(start > stop) 
      stop(paste("interval start=",interval[2]," is greater than interval end=",interval[3],sep=""))
    locus.interval <- locus.matrix[locus.matrix[,"CHROMOSOME"] == chr & locus.matrix[,"POSITION.bp"] >= start &
                                     locus.matrix[,"POSITION.bp"] <= stop,]
    if(nrow(locus.interval) == 0)
      stop("no SNPs in interval")
  } else {
    stopifnot(length(snp) == 1)
    stopifnot(snp %in% locus.matrix[,"MARKER"])
    locus.interval <- locus.matrix[locus.matrix[,"MARKER"] == snp,]
  } 
  
  # vector of binary IBD across all pairs
  # 1 indicates some amount of sharing between pairs over region
  binary.ibd.matrix <- matrix(nrow=nrow(pedigree), ncol=nrow(pedigree))
  if (nrow(locus.interval) > 1) {
    # function to determine whether a pair shares a portion of the specified region IBD
    fun <- function(pair.binary.ibd, positions, prop) {
      if (prop == 0) {
        if (any(pair.binary.ibd == 1)) {
          return(1) 
        } else 
          return(0)
      }
      interval.length <- max(positions) - min(positions)
      if (all(pair.binary.ibd == 0)) {
        pair.ibd.length <- 0 
      } else 
        pair.ibd.length <- max(positions[pair.binary.ibd==1]) - min(positions[pair.binary.ibd==1])
      if ((pair.ibd.length/interval.length) >= prop) {
        return(1) 
      } else 
        return(0)
    }
    ibd.vector <- apply(locus.interval[,5:ncol(locus.interval)], 2, fun, positions=locus.interval[,"POSITION.bp"], prop=prop)
  } else {
    ibd.vector <- NULL
    for (i in 5:ncol(locus.interval)) {
      if (locus.interval[,i] == 1){
        ibd.vector <- c(ibd.vector, 1) 
      } else  
        ibd.vector <- c(ibd.vector, 0)
    }
  }
  
  # turn vector into a square matrix of binary IBD
  j <- 1
  for (i in 1:(nrow(pedigree)-1)) {
    binary.ibd.matrix[i,(i+1):nrow(pedigree)] <- ibd.vector[j:(j + nrow(pedigree) - i - 1)]
    binary.ibd.matrix[(i+1):nrow(pedigree),i] <- ibd.vector[j:(j + nrow(pedigree) - i - 1)]
    binary.ibd.matrix[i,i] <- 0
    j <- j + nrow(pedigree) - i
  }
  binary.ibd.matrix[i+1,i+1] <- 0
  
  # name columns
  isolate.names <- paste(pedigree[,1],pedigree[,2],sep="/")
  colnames(binary.ibd.matrix) <- isolate.names
  rownames(binary.ibd.matrix) <- isolate.names
  
  # create an igraph data frame with clusters
  # lower triangle of -1
  binary.ibd.matrix[lower.tri(binary.ibd.matrix)] <- -1
  
  # reshape the data using the function melt
  links <- melt(binary.ibd.matrix)
  colnames(links) <- c("from","to","value")
  
  # select only pairs with IBD inferred
  links <- links[links[,3] == 1,]
  
  # get nodes and their original name (from pedigree)
  nodes <- unique(c(as.character(links[,1]),as.character(links[,2])))
  
  # create an igraph dataframe
  net <- graph.data.frame(links, nodes, directed=F)
  
  # get clusters
  if (hclust) {
    # hierarchical clustering
    clusters.h <- fastgreedy.community(net)
    
    # create a list of clusters
    cluster.list <- list()
    for (clust in 1:length(clusters.h)) {
      cluster.list[[clust]] <- clusters.h[[clust]]
    }
  } else {
    # disjoint clustering
    clusters <- components(net)
    cluster.belong <- clusters[[1]]
    
    # create a list of each cluster
    cluster.list <- list()
    for (clust in 1:clusters[[3]]) {
      cluster.list[[clust]] <- names(cluster.belong[cluster.belong == clust])
    }
  }
  
  # print summary info
  cluster.list.1 <- cluster.list[order(sapply(cluster.list,length),decreasing=T)]
  clusterSummary(cluster.list.1)
  
  #return(list(cluster.list.1, binary.ibd.matrix))
  return(list(cluster.list.1, net, hclust))
}


#' A function that prints summary information of clusters 
#' @param cluster.list a list of isolate groups where unique elements correspond to isolates in unique clusters
clusterSummary <- function(cluster.list) {
  cat(paste("Number of clusters = ",length(cluster.list),"\n",sep=""))
  cat(paste("Maximum cluster length = ",length(cluster.list[[1]]),"\n",sep=""))
  cat(paste("Minimum cluster length = ",length(cluster.list[[length(cluster.list)]]),"\n",sep=""))
}

