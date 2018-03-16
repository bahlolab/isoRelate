#' Plot Cluster Networks
#'
#' \code{plotIBDclusters()} Produces a figure of an isoRelate cluster network, where unique isolates are represented by vertices and a line is
#' drawn between two vertices if the isolates have been inferred IBD via the criteria specified in either \code{getIBDiclusters} or
#' \code{getIBDpclusters}.The networks are created using the R package \code{igraph}.
#'
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param clusters A named list of three objects containing network information.
#' See the \code{Value} description in either \code{\link{getIBDiclusters}} or \code{\link{getIBDpclusters}} for more details on this input.
#' @param groups A data frame with 3 columns of information:
#' \enumerate{
#' \item Family ID
#' \item Isolate ID
#' \item Group ID
#' }
#' Group ID, for example, can be the geographic regions where the isolates were collected.
#' If \code{groups} is specified then each isolate in the pedigree must belong to a group.
#' Vertices in the network will be colored according to group allocation.
#' The default is \code{groups=NULL} and all vertices will have the same color.
#' @param vertex.color A vector of characters or numeric values of the vertex colors in the network.
#' If \code{groups} is specified then \code{vertex.color} should contain the same number of colors as unique groups.
#' @param vertex.frame.color Character string or numeric value. A single color that will be used as the vertex border. Default is \code{vertex.fame.color="white"}.
#' @param vertex.size Numeric value indicating the size of the vertices in the network. Default is \code{vertex.size=4}.
#' @param vertex.name Logical. Whether to add isolate names to the vertices. Default is \code{vertex.name=FALSE}.
#' @param edge.color Character string or numeric value. A single color to be used for all edges. Default is \code{edge.color="gray60"}.
#' @param edge.width Numeric. A single value indicating the width of the edges. Default is \code{edge.width=0.8}.
#' @param mark.border Character string or numeric value. A single color to be used for all borders in hierarchical clustering groups. Default is \code{mark.border="white"}.
#' @param mark.col Character string or numeric value. A single color to be used to fill hierarchical clustering groupings. Default is \code{mark.col="gray94"}.
#' @param add.legend Logical. Whether to include a legend in the plot. Default is \code{add.legend=TRUE}.
#' @param legend.x Numerical. A single value indicating the x-coordinate of the legend, with default \code{legend.x=-1.5}.
#' @param legend.y Numerical. A single value indicating the y-coordinate of the legend, with default \code{legend.y=-0.25}.
#' @param layout A matrix containing the x and y coordinates of the vertices, generated using the Fruchterman-Reingold force-directed layout.
#' @param return.layout Logical. Whether or not to return the layout matrix (vertex positions) in the network.
#' This layout can be used as the input for the parameter \code{layout} to avoid different network configurations each time \code{plotClusters()}
#' is run on the same network.
#' @importFrom igraph graph.data.frame is.igraph E E<- V V<- layout_with_fr fastgreedy.community components
#' @importFrom ggnetwork ggnetwork
#' @importFrom graphics plot legend
#' @export
#' @seealso \code{\link{getGenotypes}}, \code{\link{getIBDpclusters}} and \code{\link{getIBDiclusters}}.
#' @examples
#' # generate the isolates who are IBD over the Plasmodium falciparum CRT gene
#' my_i_clusters <- getIBDiclusters(ped.genotypes = png_genotypes,
#'                                  ibd.segments = png_ibd,
#'                                  interval = c("Pf3D7_07_v3", 403222, 406317),
#'                                  prop=0,
#'                                  hi.clust = FALSE)
#'
#' str(my_i_clusters)
#'
#' # creating a stratification dataset
#' my_groups <- png_genotypes[[1]][,1:3]
#' my_groups[1:10,"pid"] <- "a"
#' my_groups[11:25,"pid"] <- "b"
#' my_groups[26:38,"pid"] <- "c"
#'
#' # plot the network of clusters
#' plotIBDclusters(ped.genotypes = png_genotypes,
#'                 clusters = my_i_clusters,
#'                 groups = my_groups,
#'                 vertex.color = NULL,
#'                 vertex.frame.color = "white",
#'                 vertex.size = 4,
#'                 vertex.name = FALSE,
#'                 edge.color = "gray60",
#'                 edge.width = 0.8,
#'                 mark.border = "white",
#'                 mark.col = "gray94",
#'                 add.legend = TRUE,
#'                 legend.x = -1.5,
#'                 legend.y = -0.25,
#'                 layout = NULL,
#'                 return.layout = FALSE)
#'
#' # generate the isolates who share at least than 90% of their genome IBD
#' my_p_clusters <- getIBDpclusters(ped.genotypes = png_genotypes,
#'                                  ibd.segments = png_ibd,
#'                                  prop=0.9,
#'                                  hi.clust = FALSE)
#'
#' # plot the network of clusters
#' plotIBDclusters(ped.genotypes = png_genotypes,
#'                 clusters = my_p_clusters,
#'                 groups = my_groups,
#'                 vertex.color = NULL,
#'                 vertex.frame.color = "white",
#'                 vertex.size = 4,
#'                 vertex.name = FALSE,
#'                 edge.color = "gray60",
#'                 edge.width = 0.8,
#'                 mark.border = "white",
#'                 mark.col = "gray94",
#'                 add.legend = TRUE,
#'                 legend.x = -1.5,
#'                 legend.y = -0.25,
#'                 layout = NULL,
#'                 return.layout = FALSE)
plotIBDclusters <- function(ped.genotypes, clusters, groups = NULL, vertex.color = NULL, vertex.frame.color = "white",
                         vertex.size = 4, vertex.name = FALSE, edge.color = "gray60", edge.width = 0.8, mark.border = "white",
                         mark.col = "gray94", add.legend = TRUE, legend.x = -1.5, legend.y = -0.25, layout = NULL, return.layout = FALSE){

  # check format of input data
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) stop ("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) stop ("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  if (!is.data.frame(pedigree)) stop ("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6) stop ("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", "moi", "aff")))
    stop ("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")

  # check cluster format
  if (!is.list(clusters) | length(clusters) != 3) stop ("'clusters' must be a list containing 3 objects: 'clusters', 'i.network' and 'hi.clust'")
  if (any(names(clusters) != c("clusters", "i.network", "hi.clust"))) stop ("'clusters' must be a list containing 3 objects: 'clusters', 'i.network' and 'hi.clust'")
  if (!is.list(clusters[[1]])) stop ("'clusters[[1]]' must be a list containing the network 'clusters'")
  if (!is.igraph(clusters[[2]])) stop ("'clusters[[2]]' must be an igraph network")
  if (!is.logical(clusters[[3]])) stop ("'clusters[[3]]' logical")
  i.network <- clusters[[2]]
  hi.clust <- clusters[[3]]

  # check null arguments
  if (!is.null(groups)) {
    if (!is.data.frame(groups)) stop ("'groups' is not a data.frame")
    if (ncol(groups) != 3) stop ("'groups' has incorrect format - must have 3 columns: fid, iid, group")
    colnames(groups)[1:2] <- c("fid","iid")
  }

  # check vertex colours
  if (!is.null(vertex.color)) {
    if (!is.vector(vertex.color)) stop ("'vertex.color' has incorrect format - must be a vector")
    if (!is.character(vertex.color) & !is.numeric(vertex.color)) stop ("'vertex.color' must be of type 'character' or 'numeric'")
    if (!all(areColors(vertex.color))) stop ("some 'vertex.color' are not valid colors")
  }
  if (!is.vector(vertex.frame.color)) stop ("'vertex.frame.color' has incorrect format - must be a vector")
  if (!is.character(vertex.frame.color) & !is.numeric(vertex.frame.color)) stop ("'vertex.frame.color' must be of type 'character' or 'numeric'")
  if (length(vertex.frame.color) != 1) stop ("'vertex.frame.color' must be of length 1")
  if (!all(areColors(vertex.frame.color))) stop ("'vertex.frame.color' is not a valid colur")

  # check numeric arguments
  if (!is.vector(vertex.size)) stop ("'vertex.size' has incorrect format - must be a vector")
  if (!is.numeric(vertex.size)) stop ("'vertex.size' has incorrect format - must be a single numeric value")
  if (length(vertex.size) != 1) stop ("'vertex.size' has incorrect format - must be a single numeric value")

  # check vertex.name
  if (!is.vector(vertex.name)) stop ("'vertex.name' has incorrect format - must be a vector")
  if (!is.logical(vertex.name)) stop ("'vertex.name' has incorrect format - must be a single logical value")
  if (length(vertex.name) != 1) stop ("'vertex.name' has incorrect format - must be a single logical value")

  # check add.legend
  if (!is.vector(add.legend)) stop ("'add.legend' has incorrect format - must be a vector")
  if (!is.logical(add.legend)) stop ("'add.legend' has incorrect format - must be a single logical value")
  if (length(add.legend) != 1) stop ("'add.legend' has incorrect format - must be a single logical value")

  # check edge colours
  if (!is.vector(edge.color)) stop ("'edge.color' has incorrect format - must be a vector")
  if (!is.character(edge.color) & !is.numeric(edge.color)) stop ("'edge.color' must be of type 'character' or 'numeric'")
  if (length(edge.color) != 1) stop ("'edge.color' must be of length 1")
  if (!all(areColors(edge.color))) stop ("'edge.color' is not a valid color")

  # check edge width
  if (!is.vector(edge.width)) stop ("'edge.width' has incorrect format - must be a vector")
  if (!is.numeric(edge.width) | length(edge.width) != 1) stop ("'edge.width' has incorrect format - must be a single numeric value")
  if (length(edge.width) != 1) stop ("'edge.width' has incorrect format - must be a single numeric value")

  # check edge colours
  if (!is.vector(mark.border)) stop ("'mark.border' has incorrect format - must be a vector")
  if (!is.character(mark.border) & !is.numeric(mark.border)) stop ("'mark.border' must be of type 'character' or 'numeric'")
  if (!all(areColors(mark.border))) stop ("'edge.color' is not a valid color")
  if (!is.vector(mark.col)) stop ("'mark.col' has incorrect format - must be a vector")
  if (!is.character(mark.col) & !is.numeric(mark.col)) stop ("'mark.col' must be of type 'character' or 'numeric'")
  #if (length(mark.col) != 1) stop ("'mark.col' must be of length 1")
  if (!all(areColors(mark.col))) stop ("'mark.col' is not a valid color")

  # check legend xy coordinates
  if (!is.vector(legend.x)) stop ("'legend.x' has incorrect format - must be a vector")
  if (!is.numeric(legend.x) | length(legend.x) != 1) stop ("'legend.x' has incorrect format - must be a single numeric value")
  if (!is.vector(legend.y)) stop ("'legend.y' has incorrect format - must be a vector")
  if (!is.numeric(legend.y) | length(legend.y) != 1) stop ("'legend.y' has incorrect format - must be a single numeric value")

  # check layout
  if (!is.null(layout)) {
    if (!is.matrix(layout)) stop("'layout' has incorrect format - must be a data.frame")
    if (ncol(layout) != 2) stop("'layout' has incorrect format - must have 2 columns")
  }

  # check return.layout
  if (!is.vector(return.layout)) stop ("'return.layout' has incorrect format - must be a vector")
  if (!is.logical(return.layout)) stop ("'return.layout' is not logical")

  # node names
  nodes <- V(i.network)$name
  if (vertex.name) {
    node.names <- nodes
  } else {
    node.names <- ""
  }
  if (length(nodes) <= 1) stop ("cannot plot network with a single vertex")

  # check layout has same number of nodes
  if (!is.null(layout)) {
    if (nrow(layout) != length(nodes)) stop ("'layout' has incorrect format - must have same number of vertices as in 'clusters'")
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
  #node.shape[moi == 1] <- "fcircle"
  #node.shape[moi == 2] <- "fsquare"

  # legend parameters
  legend.labels <- NULL
  legend.col    <- NULL

  # assign colour to groups (id, group1, group2)
  if (!is.null(groups)) {
    # check the isolates in igraph belong to a group
    group.names <- paste(groups[,1],groups[,2],sep="/")
    if (!all(nodes %in% group.names)) stop ("'groups' has incorrect format - some isolates 'ped.genotypes' are missing from 'groups'")

    # subset group by nodes in igraph - order is important here
    groups.0 <- NULL
    for (node.0 in nodes) {
      groups.0 <- rbind(groups.0, groups[group.names == node.0,])
    }

    # assign colours
    groups.1 <- as.character(unique(groups.0[,3]))
    if (!is.null(vertex.color)) {
      if (length(vertex.color) == length(groups.1)) {
        col.group.1 <- vertex.color
      }
      if (length(vertex.color) < length(groups.1)) {
        stop (paste0("'vertex.color' has incorrect format - must have at least ",length(groups.1)," colors specified"))
      }
      if (length(vertex.color) > length(groups.1)) {
        warning (paste0("'vertex.color' has ",length(vertex.color)," colors specified and requires ",
                        length(groups.1)," colors. Using first ",length(groups.1), " colors."))
        col.group.1 <- vertex.color[1:length(groups.1)]
      }
    } else {
      col.group.1 <- getColourPaletteMajor(length(groups.1))
    }
    col.g <- NULL
    for (i in 1:length(groups.1)) {
      col.g[groups.0[,3] == groups.1[i]] <- col.group.1[i]
    }
    groups.0 <- cbind(groups.0, col.g)
    legend.labels <- c(legend.labels, groups.1)
    legend.col    <- c(legend.col, col.group.1)
  } else{
    if (!is.null(vertex.color)) {
      if (length(vertex.color) > 1)
        warning (paste0("'vertex.color' has ",length(vertex.color)," colors specified and requires 1 color. Using first color only."))
      groups.0 <- data.frame(col.g=rep(vertex.color[1],length(nodes)))
    } else
      groups.0 <- data.frame(col.g=rep("grey",length(nodes)))
  }
  groups.0[,"col.g"] <- as.character(groups.0[,"col.g"])

  # get coordinates of large network (Fruchterman-Reingold force-directed layout)
  if (is.null(layout)) {
    my.network <- ggnetwork::ggnetwork(i.network)
    my.network.i.network <- my.network[!duplicated(my.network[,c(1,2,4)]),]
    layout <- as.matrix(my.network.i.network[,1:2], ncol=2)
  }

  # plot networks
  if (clusters[[3]]) {
    # hierarchical clustering
    clusters.h <- igraph::fastgreedy.community(i.network)
    plot(clusters.h, i.network,
         col=groups.0[,"col.g"],
         vertex.size=vertex.size,
         vertex.shape=node.shape,
         vertex.label=node.names,
         vertex.frame.color=vertex.frame.color,
         edge.width=edge.width,
         edge.color=edge.color,
         mark.border=mark.border,
         mark.col=mark.col,
         layout=layout)
  } else {
    plot(i.network,
         vertex.color=groups.0[,"col.g"],
         vertex.size=vertex.size,
         vertex.shape=node.shape,
         vertex.label=node.names,
         vertex.frame.color=vertex.frame.color,
         edge.width=edge.width,
         edge.color=edge.color,
         layout=layout)
  }

  # legend
  if(add.legend & any(moi == 2) & !is.null(legend.labels) & !is.null(legend.col)){
    legend(x=legend.x, y=legend.y, c(legend.labels, "MOI=1", "MOI>1"), pch=c(rep(21, length(legend.labels)), 21, 22),
           col="gray85", pt.bg=c(legend.col, "gray72", "gray72"), pt.cex=2, cex=.8, bty="n", ncol=1)
  } else if(add.legend & all(moi == 1) & !is.null(legend.labels) & !is.null(legend.col)){
    legend(x=legend.x, y=legend.y, legend.labels, pch=21,
           col="#777777", pt.bg=legend.col, pt.cex=2, cex=.8, bty="n", ncol=1)
  } else if(add.legend & any(moi == 2) & is.null(legend.labels) & is.null(legend.col)){
    legend(x=legend.x, y=legend.y, c("MOI=1","MOI>1"), pch=c(21,22), col="#777777",
           pt.bg="gray72", pt.cex=2, cex=.8, bty="n", ncol=1)
  }

  if (return.layout) return(as.matrix(layout))
}

