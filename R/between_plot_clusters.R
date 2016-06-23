#' Between Isolates Clusters
#' 
#' \code{getClusters()} identifies clusters of isolates that share IBD over part of the entire specified region.
#' Each cluster represents a disjoint network where each isolate in the network is conneceted to all other isolates
#' in the network, either directly or through a series of connected isoaltes.
#' @param ped.genotypes a list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param clusters a data frame containing the binary IBD information for each SNP and each pair. 
#' See the returned \code{value} in \code{\link{between_locus_matrix}} for more details.
#' @param groups a character string of length 1 giving the name of a SNP to find clusters of isolates who are IBD here.
#' Note that \code{snpName} and \code{interval} are mutually exclusive and only one option should be specified.
#' @param groups a vector of length 3 containing the region to identify clusters over. This vector should contain the
#' chromosome ID, the start of the interval in bp and the end of the interval in bp; in this order respectively. Note that 
#' \code{snpName} and \code{interval} are mutually exclusive and only one option should be specified.
#' @param vertex.size numeric value between 0 and 1. The minimum proportion of an interval (in bp) shared IBD between a pair of isolates in order for the pair 
#' to be directly connected in a cluster. For example, if \code{prop=1} then two isolates will only be connected if they are IBD over the
#' entire specified region, whereas if \code{prop=0.5} then two isolates will be connected if they are IBD over at least half of the interval.
#' The default is \code{prop=0}, which connects two isolates if they are IBD at all over the interval.
#' @param vertex.name logical. Whether to add isolate names to the nodes. Default is \code{vertex.name=FALSE}.
#' @param add.legend logical. Whether to plot a legend. This only applies when isolates have MOI > 1 or \code{groups} is specified.
#' @param vertex.frame.color character string. A single colour that outlines vertices. Default is \code{vertex.fame.color="white"}.
#' @param ... other plotting parameters allowed by \code{igraph}.
#' @return A list of disjoint clusters that contain the names of isoaltes that are part of the same network.
#' The last object in the list is the binary IBD matrix for the specified region, where 1 represents IBD sharing by
#' pairs over the specified interval and 0 otherwise.
#' @importFrom igraph graph.data.frame is.igraph E E<- V V<- layout_with_fr fastgreedy.community components
#' @export
plotClusters <- function(ped.genotypes, clusters, groups = NULL, vertex.size = 4, vertex.name = FALSE, 
                         vertex.frame.color = "white", add.legend = TRUE, ...){
  
  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 2)
  pedigree <- ped.genotypes[[1]] 
  
  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  
  # check cluster format
  stopifnot(is.list(clusters) | length(clusters) == 3)
  if (!is.list(clusters[[1]]) | !is.igraph(clusters[[2]]) | !is.logical(clusters[[3]]))
    stop("clusters has incorrect format")
  
  # check null arguments
  if (!is.null(groups))
    stopifnot(is.data.frame(groups))
  
  # check numeric arguments
  stopifnot(is.numeric(vertex.size))
  
  # check logical
  stopifnot(is.logical(vertex.name))
  stopifnot(is.logical(add.legend))
  
  # check characters
  stopifnot(is.character(vertex.frame.color) | is.numeric(vertex.frame.color))
  
  # get igraph network
  net <- clusters[[2]]
  
  # node names
  nodes <- V(net)$name
  if (vertex.name) {
    node.names <- nodes
  } else {
    node.names <- ""
  }
  
  # get moi - order of nodes is not that in pedigree
  isolate.names <- paste(pedigree[,1], pedigree[,2], sep="/")
  moi <- NULL
  for (node.0 in nodes) {
    moi <- c(moi, pedigree[isolate.names == node.0,"moi"])
  }
  node.shape <- NULL
  node.shape[moi == 1] <- "circle"
  node.shape[moi == 2] <- "square"
  
  # legend parameters
  legend.labels <- NULL 
  legend.col    <- NULL

  # assign colour to groups (id, group1, group2)
  if (!is.null(groups)) {
    stopifnot(ncol(groups) > 1)
    
    # check the isolates in igraph belong to a group
    if (!is.null(groups)) {
      group.names <- paste(groups[,1],groups[,2],sep="/")
      if (!all(nodes %in% group.names)) 
        stop("groups is missing information for some isoaltes")
    }
    
    # subset group by nodes in igraph
    groups.0 <- NULL
    for (node.0 in nodes) {
      groups.0 <- rbind(groups.0, groups[group.names == node.0,])
    }
    
    # assign colours
    if (ncol(groups.0) > 4){
      cat("using first 4 columns of groups")
      groups.0 <- groups.0[,1:4]
    }
    if (ncol(groups.0) == 4) {
      groups.1    <- as.character(unique(groups.0[,3])) # groups
      col.group.1 <- getColourPaletteMajor(length(groups.1)) # unique group colour
      col.g      <- NULL
      for (i in 1:length(groups.1)) {
        sub_groups <- paste(groups.0[groups.0[,3] == groups.1[i],3],
                            groups.0[groups.0[,3] == groups.1[i],4], sep=", ") # sub groups
        groups.2    <- unique(sub_groups)
        col.group.2 <- getColourPaletteMinor(ol.group.1[i],length(groups.2))
        for (j in 1:length(groups.2)) {
          col.g[paste(groups.0[,3], groups.0[,4],sep=", ") == groups.2[j]] <- col.group.2[j]
        }
        legend.labels <- c(legend.labels, groups.2)
        legend.col    <- c(legend.col, col.group.2[1:length(groups.2)])
      }
      groups.0 <- cbind(groups.0, col.g)
    }
    if (ncol(groups.0) == 3) {
      groups.1    <- as.character(unique(groups.0[,3]))
      col.group.1 <- getColourPaletteMajor(length(groups.1))
      col.g      <- NULL
      for (i in 1:length(groups.1)) {
        col.g[groups.0[,3] == groups.1[i]] <- col.group.1[i]
      }
      groups.0 <- cbind(groups.0, col.g)
      legend.labels <- c(legend.labels, groups.1)
      legend.col    <- c(legend.col, col.group.1)
    }
  } else{
    groups.0 <- data.frame(col.g=rep(colours(1),length(nodes)))
  }
  groups.0[,"col.g"] <- as.character(groups.0[,"col.g"])
  
  
  # plot networks
  if (clusters[[3]]) {
    # hierachical clustering
    clusters.h <- fastgreedy.community(net)
    plot(clusters.h, net, 
         col=groups.0[,"col.g"], 
         vertex.size=vertex.size,
         vertex.shape=node.shape,
         vertex.label=node.names,
         vertex.frame.color=vertex.frame.color,
         edge.color="gray60",
         mark.border="white",
         mark.col="gray94",
         ...)
  } else {
    plot(net, 
         vertex.color=groups.0[,"col.g"], 
         vertex.size=vertex.size,
         vertex.shape=node.shape,
         vertex.label=node.names,
         vertex.frame.color=vertex.frame.color,
         edge.color="gray60",
         ...)
  }
  
  # legend
  if(any(moi == 2) & add.legend){
    legend(x=-1.5, y=-0.65, c("MOI=1","MOI>1"), pch=c(21,22), col="#777777", 
           pt.bg="gray72", pt.cex=2, cex=.8, bty="n", ncol=1)
  }
  if(!is.null(legend.labels) & !is.null(legend.col) & add.legend) { 
    legend(x=-1.5, y=-0.9, legend.labels, pch=21,col="#777777", pt.bg=legend.col, 
           pt.cex=2, cex=.8, bty="n", ncol=1)
  }
}