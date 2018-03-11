#' Plot iR Statistics
#'
#' \code{plotIBDiR()} plots the -log10 (p-values) used to assess the significance of excess IBD sharing.
#'
#' @param ibd.iR A data frame containing the iR summary statistics for each SNP.
#' See the returned \code{Value} in \code{\link{getIBDiR}} for more details.
#' If multiple subpopulations are specified (column name "subpop") then the iR statistics for each subpopulation will be plotted
#' on a separate facet in the figure (see \url{http://docs.ggplot2.org/current/facet_grid.html} on faceting).
#' If there are many subpopulations (>8) it may be better to plot subsets of the subpopulations as apposed to all subpopulations in a single figure.
#' If multiple populations (column name "pop") are specified then only the first population will be included in the figure.
#' Genomic locations of annotation genes can be included in the figure and specific intervals can be highlighted.
#' @param interval A vector of length 3 containing the genomic locations of a specific region to plot. This vector should contain the
#' chromosome ID, the start of the interval in base-pairs and the end of the interval in base-pairs; in this order respectively.
#' The default is \code{interval=NULL} which will plot iR statistics over all chromosomes in \code{ibd.iR}.
#' @param annotation.genes A data frame containing information on annotation genes to be included in the figure.
#' This data frame must have at least 5 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{"numeric"} or \code{"integer"})
#' \item Gene name (type \code{"character"})
#' \item Start location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item End location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item Gene strand (+ or -) (type \code{"character"})
#' }
#' \code{annotation.genes} must contain the following headers \code{chr, name, start, end} and \code{strand}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels. The default is {annotation.genes=NULL}.
#' @param annotation.genes.color A vector of characters or numeric values containing the two colors according to gene stand (positive (+) or negative (-))
#' @param highlight.genes A data frame containing information of genes or regions to highlight.
#' The data frame must have at least 4 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{"numeric"} or \code{"integer"})
#' \item Gene name (type \code{"character"})
#' \item Start location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' \item End location of the gene in base-pairs (type \code{"numeric"} or \code{"integer"})
#' }
#' \code{highlight.genes} should contain the following headers \code{chr, name, start} and \code{end}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels. The default is {highlight.genes=NULL}.
#' @param highlight.genes.labels Logical. Whether to include gene names as labels in the figure. The default is \code{highlight.genes.labels=FALSE}.
#' @param highlight.genes.color Character string or numeric value. A single color that will be used to highlight a region/gene. The default is \code{highlight.genes.color=NULL}.
#' @param highlight.genes.alpha Numeric. A single value between 0 and 1 indicating the gene color transparency. The default is \code{highlight.genes.alpha=0.1}.
#' @param point.size Numeric. The size of the points in the figures. The default is \code{point.size=1}.
#' @param point.color A vector of characters or numeric values denoting the color of points to be plotted.
#' If there are multiple subpopulations then the number of colors specified should equal the number of subpopulations.
#' The default is \code{point.color=NULL} which will use isoRelate default colors.
#' @param add.rug Logical. Whether to include SNP positions as a rug in the figure. The default is \code{add.rug=FALSE}
#' @param plot.title A character string of a title to be added to the figure The default is \code{plot.title=NULL} which does not add a title to the plot.
#' @param add.legend Logical. Whether a legend containing subpopulation information should be plotted. The default is \code{add.legend=FALSE}.
#' @param facet.label Logical. Whether to include facet labels if multiple subpopulations (column name "subpop") are specified.
#' @param facet.scales A character string of either \code{"fixed"}, \code{"free"}, \code{"free_x"} or \code{"free_y"} specifying the facet axis-scales.
#' See \url{http://docs.ggplot2.org/current/facet_grid.htmlfacet.scales} for more information on this parameter. The default is \code{facet.scales="fixed"}
#' @import ggplot2
#' @export
plotIBDiR <- function(ibd.iR, interval = NULL, annotation.genes = NULL, annotation.genes.color = NULL,
                   highlight.genes = NULL, highlight.genes.labels = TRUE, highlight.genes.color = NULL, highlight.genes.alpha = 0.1,
                   point.size = 1, point.color = NULL, add.rug = FALSE, plot.title = NULL, add.legend = FALSE, facet.label = TRUE,
                   facet.scales = "fixed"){

  # check locus matrix
  if (!is.data.frame(ibd.iR)) stop ("'ibd.iR' has incorrect format - must be a data.frame")
  if (ncol(ibd.iR) != 8)
    stop ("'ibd.iR' has incorrect format - must be a data.frame with 8 columns: chr, snp_id, pos_M, pos_bp, pop, subpop, iR and log10_pvalue")
  if (any(colnames(ibd.iR) != c("chr", "snp_id", "pos_M", "pos_bp", "pop", "subpop", "iR", "log10_pvalue")))
    stop ("'ibd.iR' has incorrect format - must be a data.frame with 8 columns: chr, snp_id, pos_M, pos_bp, pop, subpop, iR and log10_pvalue")
  pops <- as.character(unique(ibd.iR[,"pop"]))
  if (length(unique(pops)) > 1) warning ("can only plot one population at a time - plotting pop=",unique(pops)[1]," only")
  ibd.iR <- ibd.iR[ibd.iR[,"pop"] == unique(pops)[1],]
  subpops <- as.character(unique(ibd.iR[,"subpop"]))

  # check interval
  if (!is.null(interval)) {
    if (!is.vector(interval)) stop ("'interval' has incorrect format - must be a vector of length 3")
    if (length(interval) != 3) stop ("'interval' has incorrect format - must be a vector of length 3")
    interval.chr   <- as.character(interval[1])
    interval.start <- as.numeric(interval[2])
    interval.stop  <- as.numeric(interval[3])
    if (interval.start > interval.stop) stop(paste("interval start=",interval.start," is greater than interval end=",interval.stop,sep=""))
  }

  # check annotation genes
  if (!is.null(annotation.genes)) {
    if (!is.data.frame(annotation.genes)) stop ("'annotation.genes' has incorrect format - must be a data.frame")
    if (!all(c("name","strand","chr","start","end") %in% colnames(annotation.genes)))
      stop ("'annotation.genes' has incorrect format - must have columns: name, strand, chr, start and end")
    if (is.factor(annotation.genes[,"name"])) annotation.genes[,"name"] <- as.character(annotation.genes[,"name"])
    if (is.factor(annotation.genes[,"chr"])) annotation.genes[,"chr"] <- as.character(annotation.genes[,"chr"])
    if (!is.numeric(annotation.genes[,"start"])) annotation.genes[,"start"] <- as.numeric(as.character(annotation.genes[,"start"]))
    if (!is.numeric(annotation.genes[,"end"])) annotation.genes[,"end"] <- as.numeric(as.character(annotation.genes[,"end"]))
  }

  # check annotation.genes.color
  if (!is.null(annotation.genes.color)) {
    if (!is.vector(annotation.genes.color)) stop ("'annotation.genes.color' has incorrect format - must be a vector")
    if (!is.character(annotation.genes.color) & !is.numeric(annotation.genes.color)) stop ("'annotation.genes.color' must be of type 'character' or 'numeric'")
    if (length(annotation.genes.color) < 2) stop ("'annotation.genes.color' has incorrect format - must specify 2 colors")
    if (length(annotation.genes.color) > 2) {
      warning ("'annotation.genes.color' has ",length(annotation.genes.color)," colors specified and requires 2. Using first 2 colors only")
      annotation.genes.color <- annotation.genes.color[1:2]
    }
    if (!all(areColors(annotation.genes.color))) stop ("some 'annotation.genes.color' are not valid colors")
  } else
    annotation.genes.color <- c("gold","firebrick1")

  # check highlight genes
  if (!is.null(highlight.genes)) {
    if (!is.data.frame(highlight.genes)) stop ("'highlight.genes' has incorrect format - must be a data.frame")
    if (!all(c("name","chr","start","end") %in% colnames(highlight.genes)))
      stop ("'highlight.genes' has incorrect format - must have columns: name, chr, start and end")
    if (is.factor(highlight.genes[,"name"])) highlight.genes[,"name"] <- as.character(highlight.genes[,"name"])
    if (is.factor(highlight.genes[,"chr"])) highlight.genes[,"chr"] <- as.character(highlight.genes[,"chr"])
    if (!is.numeric(highlight.genes[,"start"])) highlight.genes[,"start"] <- as.numeric(as.character(highlight.genes[,"start"]))
    if (!is.numeric(highlight.genes[,"end"])) highlight.genes[,"end"] <- as.numeric(as.character(highlight.genes[,"end"]))
  }

  # check highlight.genes.labels
  if (!is.vector(highlight.genes.labels)) stop ("'highlight.genes.labels' has incorrect format - must be a single logical value")
  if (!is.logical(highlight.genes.labels)) stop ("'highlight.genes.labels' has incorrect format - must be a single logical value")
  if (length(highlight.genes.labels) != 1) stop ("'highlight.genes.labels' has incorrect format - must be a single logical value")

  # check highlight.genes.color
  if (!is.null(highlight.genes.color)) {
    if (!is.vector(highlight.genes.color)) stop ("'highlight.genes.color' has incorrect format - must be a vector")
    if (!is.character(highlight.genes.color) & !is.numeric(highlight.genes.color)) stop ("'highlight.genes.color' must be of type 'character' or 'numeric'")
    if (length(highlight.genes.color) < 1) stop ("'highlight.genes.color' has incorrect format - must specify 1 color")
    if (length(highlight.genes.color) > 2) {
      warning ("'highlight.genes.color' has ",length(highlight.genes.color)," colors specified and requires 1. Using first color only")
      highlight.genes.color <- highlight.genes.color[1:2]
    }
    if (!all(areColors(highlight.genes.color))) stop ("'highlight.genes.color' is not a valid color")
  } else
    highlight.genes.color <- "gray40"

  # check highlight.genes.alpha
  if (!is.vector(highlight.genes.alpha)) stop ("'highlight.genes.alpha' has incorrect format - must be a single numeric value")
  if (!is.numeric(highlight.genes.alpha)) stop ("'highlight.genes.alpha' has incorrect format - must be a single numeric value")
  if (length(highlight.genes.alpha) != 1) stop ("'highlight.genes.alpha' has incorrect format - must be a single numeric value")
  if (highlight.genes.alpha > 1 | highlight.genes.alpha <= 0) stop ("'highlight.genes.alpha' has incorrect format - must be a single numeric value between (0,1]")

  # check rug
  if (!is.vector(add.rug)) stop ("'add.rug' has incorrect format - must be a logical vector")
  if (!is.logical(add.rug)) stop ("'add.rug' has incorrect format - must be a logical vector")
  if (length(add.rug) != 1) stop ("'add.rug' has incorrect format - must be a single logical value")

  # check title
  if (!is.null(plot.title)) {
    if (!is.vector(plot.title)) stop ("'plot.title' has incorrect format - must be a character vector")
    if (!is.character(plot.title)) stop ("'plot.title' has incorrect format - must be a character vector")
  }

  # check legend
  if (!is.vector(add.legend)) stop ("'add.legend' has incorrect format - must be a logical vector")
  if (!is.logical(add.legend)) stop ("'add.legend' has incorrect format - must be a logical vector")
  if (length(add.legend) != 1) stop ("'add.legend' has incorrect format - must be a single logical value")

  # check legend
  if (!is.vector(facet.label)) stop ("'facet.label' has incorrect format - must be a logical vector")
  if (!is.logical(facet.label)) stop ("'facet.label' has incorrect format - must be a logical vector")
  if (length(facet.label) != 1) stop ("'facet.label' has incorrect format - must be a single logical value")

  # check facet.scales
  if (!is.vector(facet.scales)) stop ("'facet.scales' has incorrect format - must be a character vector")
  if (!is.character(facet.scales)) stop ("'facet.scales' has incorrect format - must be a character vector")
  if (!(facet.scales %in% c("fixed","free","free_y","free_x"))) stop ("'facet.scales' has incorrect format - must be either fixed, free, free_y or free_x")

  # point.size
  if (!is.vector(point.size)) stop ("'point.size' has incorrect format - must be a numeric vector")
  if (!is.numeric(point.size)) stop ("'point.size' has incorrect format - must be a numeric vector")
  if (length(point.size) != 1) stop ("'point.size' has incorrect format - must be a single numeric value")
  if (point.size <= 0) stop ("'point.size' has incorrect format - must be a single numeric value > 0")

  # check point colours
  if (!is.null(point.color)) {
    if (!is.vector(point.color)) stop ("'point.color' has incorrect format - must be a vector")
    if (!is.character(point.color) & !is.numeric(point.color)) stop ("'point.color' must be of type 'character' or 'numeric'")
    if (!all(areColors(point.color))) stop ("some 'point.color' are not valid colors")

    # exact amount of colors specified:
    if (length(point.color) == length(subpops))
      col.line <- point.color
    # too few colours specified:
    if (length(point.color) < length(subpops)) {
      stop (paste0("'point.color' has incorrect format - must have at least ",length(subpops)," colors specified"))
    }
    # too many colours specified:
    if (length(point.color) > length(subpops)) {
      if (length(subpops) == 1) {
        warning ("'point.color' has ",length(point.color)," colors specified and requires 1 color. Using the first color only.")
      } else
        warning (paste0("'point.color' has ",length(point.color)," colors specified and requires either ",
                        length(pops)," or ", length(subpops)," colors. Using the first ",
                        length(subpops)," colors only."))
      col.line <- point.color[1:length(subpops)]
    }
  } else {
    col.line <- getColourPaletteMajor(length(subpops))
  }

  # check same chromosomes in each data.frame
  if (is.null(interval)) {
    chromosomes <- as.character(unique(ibd.iR[,"chr"]))
  } else {
    if (!(interval.chr %in% unique(ibd.iR[,"chr"])))
      stop("'interval' chromosome is not in 'ibd.iR'")
    chromosomes <- interval.chr
  }

  # subset annotation genes by interval
  if (!is.null(annotation.genes)) {
    if (!is.null(interval)) {
      annotation.genes.chr <- annotation.genes[annotation.genes[,"chr"] == interval.chr,]
      annotation.genes.interval <- rep(0, nrow(annotation.genes.chr))
      if (nrow(annotation.genes.chr) > 0) {
        for(g in 1:nrow(annotation.genes.chr)){
          gene.overlap <- getOverlap(annotation.genes.chr[g,c("start","end")],c(interval.start,interval.stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            annotation.genes.interval[g] <- 1
        }
        annotation.genes.overlap <- annotation.genes.chr[annotation.genes.interval == 1,]
      } else
        annotation.genes.overlap <- data.frame()
    } else {
      annotation.genes.overlap <- annotation.genes[annotation.genes[,"chr"] %in% chromosomes,]
    }
    if (nrow(annotation.genes.overlap) == 0) {
      warning("no 'annotation.genes' in interval")
    }
  }

  # subset highlight genes by interval
  if (!is.null(highlight.genes)) {
    if (!is.null(interval)) {
      highlight.genes.chr <- highlight.genes[highlight.genes[,"chr"] == interval.chr,]
      highlight.genes.interval <- rep(0, nrow(highlight.genes.chr))
      if (nrow(highlight.genes.chr) > 0) {
        for(g in 1:nrow(highlight.genes.chr)){
          gene.overlap <- getOverlap(highlight.genes.chr[g,c("start","end")],c(interval.start,interval.stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            highlight.genes.interval[g] <- 1
        }
        highlight.genes.overlap <- highlight.genes.chr[highlight.genes.interval == 1,]
      } else
        highlight.genes.overlap <- data.frame()
    } else {
      highlight.genes.overlap <- highlight.genes[highlight.genes[,"chr"] %in% chromosomes,]
    }
    if (nrow(highlight.genes.overlap) == 0) {
      warning("no 'highlight.genes' in interval")
    }
  }


  # define continuous plot positions if interval not specified
  if (is.null(interval) | length(chromosomes) > 1) {
    chrstart <- 0
    chradd   <- NULL
    labelpos <- NULL
    # change genomic positions for each
    for (i in 1:length(chromosomes)) {
      maxpos <- max(ibd.iR[ibd.iR[,"chr"] == chromosomes[i],"pos_bp"])
      minpos <- min(ibd.iR[ibd.iR[,"chr"] == chromosomes[i],"pos_bp"])
      ibd.iR[ibd.iR[,"chr"] == chromosomes[i],"pos_bp"] <- ibd.iR[ibd.iR[,"chr"] == chromosomes[i],"pos_bp"] + chrstart
      labelpos[i] <- (maxpos - minpos + 2*chrstart)/2
      chradd[i]   <- chrstart
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.overlap) != 0){
          annotation.genes.overlap[annotation.genes.overlap[,"chr"] == chromosomes[i],"start"] <- annotation.genes.overlap[annotation.genes.overlap[,"chr"] == chromosomes[i],"start"] + chrstart
          annotation.genes.overlap[annotation.genes.overlap[,"chr"] == chromosomes[i],"end"] <- annotation.genes.overlap[annotation.genes.overlap[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      if (!is.null(highlight.genes)) {
        if (nrow(highlight.genes.overlap) != 0){
          highlight.genes.overlap[highlight.genes.overlap[,"chr"] == chromosomes[i],"start"] <- highlight.genes.overlap[highlight.genes.overlap[,"chr"] == chromosomes[i],"start"] + chrstart
          highlight.genes.overlap[highlight.genes.overlap[,"chr"] == chromosomes[i],"end"] <- highlight.genes.overlap[highlight.genes.overlap[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      chrstart <- chrstart + maxpos
    }
    chradd <- chradd[-1]
  } else {
    ibd.iR.temp <- ibd.iR
    ibd.iR <- ibd.iR.temp[ibd.iR.temp[,"chr"] == interval.chr & ibd.iR.temp[,"pos_bp"] >= interval.start &
                                                  ibd.iR.temp[,"pos_bp"] <= interval.stop,]
    if (nrow(ibd.iR) == 0) stop("no SNPs in 'interval'")
    pop.subpop.2 <- as.character(unique(paste(ibd.iR[,"pop"],ibd.iR[,"subpop"],sep="/")))
    if (length(pop.subpop.2) != length(subpops)) stop ("some subpops have no SNPs over interval")
  }

  # remove subpops with all NAs
  not.na <- NULL
  for (i in 1:length(subpops)) {
    if (all(!is.na(ibd.iR[ibd.iR[,"subpop"] == subpops[i],"iR"]))) {
      not.na <- c(not.na, subpops[i])
    }
  }
  if (length(not.na) == 0 & is.null(not.na))
    stop ("all subpops have NA values - cannot plot")
  if (length(not.na) < length(subpops)) {
    warning ("subpops cannot be plotting due to NA values: ", paste(subpops[!(subpops %in% not.na)], collapse=", "))
    ibd.iR <- ibd.iR[ibd.iR[,"subpop"] %in% not.na,]
  }
  col.line <- col.line[1:length(not.na)]

  # check format of data frame
  ibd.iR[,"pos_bp"] <- as.numeric(ibd.iR[,"pos_bp"])
  ibd.iR[,"pop"] <- as.character(ibd.iR[,"pop"])
  ibd.iR[,"subpop"] <- as.character(ibd.iR[,"subpop"])
  ibd.iR[,"log10_pvalue"] <- as.numeric(ibd.iR[,"log10_pvalue"])

  # plot:

  # setting up ggplot
  ggp <- ggplot()
  ggp <- ggp + geom_point(data = ibd.iR, aes_string("pos_bp", "log10_pvalue", color = "subpop"), size=point.size)
  ggp <- ggp + theme_bw()
  ggp <- ggp + ylab("-log10(P-value)")
  ggp <- ggp + theme(panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              legend.title = element_blank())

  if (add.rug)
    ggp <- ggp + geom_rug(data=ibd.iR, aes_string(x = "pos_bp"), size = 0.1, colour = "gray30")
  if (!is.null(plot.title))
    ggp <- ggp + ggtitle(plot.title)
  if (!add.legend)
    ggp <- ggp + theme(legend.position = "none")

  # add facet
  if (length(subpops) > 1)
    ggp <- ggp + facet_grid(subpop~., scales=facet.scales)

  # facet labels
  if (!facet.label)
    ggp <- ggp + theme(strip.text.y = element_blank())

  # line colours
  ggp <- ggp + scale_colour_manual(values = col.line)

  # interval:
  if (is.null(interval)) {
    ggp <- ggp + xlab("Chromosome")
    ggp <- ggp + geom_vline(xintercept = chradd, colour = "gray87", linetype = "longdash", size = 0.4)
    ggp <- ggp + scale_x_continuous(breaks = labelpos, labels = chromosomes)
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  } else {
    ggp <- ggp + xlab(paste("Chromosome", chromosomes))
  }

  # annotation.genes:
  if (!is.null(annotation.genes)) {
    if(nrow(annotation.genes.overlap) != 0) {
      #min.y <- -0.05*max(ibd.iR[,"log10_pvalue"])
      #max.y <- -0.01*max(ibd.iR[,"log10_pvalue"])
      gene.hight <- 0.05*(max(ibd.iR[,"log10_pvalue"]) - min(ibd.iR[,"log10_pvalue"]))
      max.y <- min(ibd.iR[,"log10_pvalue"]) - gene.hight*0.5
      min.y <- max.y - gene.hight
      pos.strand <- annotation.genes.overlap[annotation.genes.overlap[,"strand"] == "+",]
      neg.strand <- annotation.genes.overlap[annotation.genes.overlap[,"strand"] != "+",]
      ggp <- ggp + geom_rect(data=pos.strand, aes_string(xmin = "start", xmax = "end"), ymin = min.y, ymax = max.y, alpha = 0.9, fill = annotation.genes.color[1])
      ggp <- ggp + geom_rect(data=neg.strand, aes_string(xmin = "start", xmax = "end"), ymin = min.y, ymax = max.y, alpha = 0.9, fill = annotation.genes.color[2])
      ggp <- ggp + ylim(min.y, max(ibd.iR[,"log10_pvalue"])) # overrides facets="free"
    }
  }

  # highlight.genes:
  if (!is.null(highlight.genes)) {
    if(nrow(highlight.genes.overlap) != 0) {
      lab.pos <- min(ibd.iR[,"log10_pvalue"])*0.05
      ggp <- ggp + geom_rect(data=highlight.genes.overlap, aes_string(xmin = "start", xmax = "end"), ymin = -Inf, ymax = Inf, fill = highlight.genes.color, alpha = highlight.genes.alpha)
      ggp <- ggp + geom_vline(data=highlight.genes.overlap, aes_string(xintercept = "start"), colour = highlight.genes.color, linetype = "solid", alpha = highlight.genes.alpha)
      if (highlight.genes.labels)
        ggp <- ggp + geom_text(data=highlight.genes.overlap, aes_string(x = "start", label = "name"), y = lab.pos, colour = "gray20", angle = 90, hjust = -0.1, vjust = -0.2, size = 3, alpha = 0.6)
    }
  }

  print(ggp)
}
