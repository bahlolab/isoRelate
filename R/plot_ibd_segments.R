#' Plot IBD Segments
#'
#' \code{plotIBDsegments()} plots IBD segments for pairs across the genome. IBD segments are depicted by colored blocks.
#'
#' @param ped.genotypes A list containing 3 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.segments A data frame containing the IBD segments inferred from pairs of isolates
#' See the returned \code{value} in \code{\link{getIBDsegments}} for more details.
#' @param interval A vector of length 3 containing the genomic locations of a specific region to plot. This vector should contain the
#' chromosome ID, the start of the interval in base-pairs and the end of the interval in base-pairs; in this order respectively.
#' The default is \code{interval=NULL} which will plot the segments over all chromosomes.
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
#' @param annotation.genes.color A vector of characters or numeric values containing the two colors representing gene stand (positive (+) or negative (-))
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
#' @param segment.height A numeric value giving the hight of IBD segment blocks, such that 0 < segment.height <= 1. The default is \code{segment.height=0.5}.
#' @param segment.color A vector of characters or numeric values denoting the color of the segments to be plotted.
#' Two colors must be specified, one for segments with 1 allele IBD and one for segments with 2 alleles IBD.
#' @param number.per.page A numeric value indicating the maximum number of IBD pairs to plot in a single graphics window. The default is
#' \code{number.per.page=NULL} which will plot all IBD pairs in a single window. This may not be ideal when there are many IBD pairs. If
#' \code{number.per.page} is set, it is recommended to plot the output to a file as opposed to the R plotting window.
#' @param fid.label Logical. If \code{fid.label=TRUE}, family IDs will be included in the y-axis labels; otherwise family IDs will be omitted.
#' The default is \code{add.fid.name=TRUE}.
#' @param iid.label Logical. If \code{iid.label=TRUE}, isolate IDs will be included in the y-axis labels; otherwise isolate IDs will be omitted.
#' The default is \code{add.iid.name=TRUE}.
#' @param ylabel.size A numeric value indicating the size of the y-axis labels if drawn. The default is \code{ylabel.size=9}.
#' @param add.rug Logical. Whether to include SNP positions as a rug in the figure. The default is \code{add.rug=FALSE}
#' @param plot.title A character string of a title to be added to the figure The default is \code{plot.title=NULL} which does not add a title to the plot.
#' @param add.legend Logical. If \code{add.legend=TRUE}, a legend specifying the IBD status (1 allele IBD or 2 alleles IBD) will be included. The default is \code{add.legend=TRUE}.
#' @import ggplot2
#' @export
#' @seealso \code{\link{getGenotypes}} and \code{\link{getIBDsegments}}.
plotIBDsegments <- function (ped.genotypes, ibd.segments, interval = NULL, annotation.genes = NULL, annotation.genes.color = NULL,
                             highlight.genes = NULL, highlight.genes.labels = TRUE, highlight.genes.color = NULL,
                             highlight.genes.alpha = 0.1, segment.height = 0.5, segment.color = NULL, number.per.page = NULL,
                             fid.label = TRUE, iid.label = TRUE, ylabel.size = 9, add.rug = FALSE, plot.title = NULL, add.legend = TRUE) {

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
  if (nrow(ibd.segments) == 0) stop ("no IBD segments detected")
  if (any(colnames(ibd.segments) != c("fid1","iid1","fid2","iid2","chr","start_snp","end_snp","start_position_bp","end_position_bp",
                                      "start_position_M", "end_position_M", "number_snps", "length_bp", "length_M", "ibd_status")))
    stop (paste("'ibd.segments' has incorrect format - 'ibd.segments' must have columns: fid1, iid1, fid2, iid2, chr",
                "start_snp, end_snp, start_position_bp, end_position_bp, start_position_M, end_position_M, number_snps",
                "length_bp, length_M, ibd_status"))

  # check if region specified
  if (!is.null(interval)) {
    if (!is.vector(interval)) stop ("'interval' has incorrect format - must be a vector of length 3")
    if (length(interval) != 3) stop ("'interval' has incorrect format - must be a vector of length 3")
    interval.chr   <- as.character(interval[1])
    interval.start <- as.numeric(interval[2])
    interval.stop  <- as.numeric(interval[3])
    if (!(interval.chr %in% ibd.segments[,"chr"])) stop (paste("no pairs IBD on chromosome =",interval.chr))
    if (!(interval.chr %in% genotypes[,"chr"])) stop (paste("chromosome =",interval.chr,"is not in 'ped.genotypes'"))
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

  # check segment.height
  if (!is.vector(segment.height)) stop ("'segment.height' has incorrect format - must be a single numeric value")
  if (!is.numeric(segment.height)) stop ("'segment.height' has incorrect format - must be a single numeric value")
  if (length(segment.height) != 1) stop ("'segment.height' has incorrect format - must be a single numeric value")
  if (segment.height > 1 | segment.height <= 0) stop ("'segment.height' has incorrect format - must be a single numeric value between (0,1]")

  # check number.per.page
  if (!is.null(number.per.page)) {
    if (!is.vector(number.per.page)) stop ("'number.per.page' has incorrect format - must be a single numeric value")
    if (!is.numeric(number.per.page)) stop ("'number.per.page' has incorrect format - must be a single numeric value")
    if (length(number.per.page) != 1) stop ("'number.per.page' has incorrect format - must be a single numeric value")
    if (number.per.page <= 0) stop ("'number.per.page' has incorrect format - must be a single numeric value > 0")
  }

  # check fid.label
  if (!is.vector(fid.label)) stop ("'fid.label' has incorrect format - must be a single logical value")
  if (!is.logical(fid.label)) stop ("'fid.label' has incorrect format - must be a single logical value")
  if (length(fid.label) != 1) stop ("'fid.label' has incorrect format - must be a single logical value")

  # check iid.label
  if (!is.vector(iid.label)) stop ("'iid.label' has incorrect format - must be a single logical value")
  if (!is.logical(iid.label)) stop ("'iid.label' has incorrect format - must be a single logical value")
  if (length(iid.label) != 1) stop ("'iid.label' has incorrect format - must be a single logical value")

  # check ylabel.size
  if (!is.vector(ylabel.size)) stop ("'ylabel.size' has incorrect format - must be a single numeric value")
  if (!is.numeric(ylabel.size)) stop ("'ylabel.size' has incorrect format - must be a single numeric value")
  if (length(ylabel.size) != 1) stop ("'ylabel.size' has incorrect format - must be a single numeric value")

  # check rug
  if (!is.vector(add.rug)) stop ("'add.rug' has incorrect format - must be a single logical value")
  if (!is.logical(add.rug)) stop ("'add.rug' has incorrect format - must be a single logical value")
  if (length(add.rug) != 1) stop ("'add.rug' has incorrect format - must be a single logical value")

  # check title
  if (!is.null(plot.title)) {
    if (!is.vector(plot.title)) stop ("'plot.title' has incorrect format - must be a character vector")
    if (!is.character(plot.title)) stop ("'plot.title' has incorrect format - must be a character vector")
  }

  # check legend
  if (!is.vector(add.legend)) stop ("'add.legend' has incorrect format - must be a single logical value")
  if (!is.logical(add.legend)) stop ("'add.legend' has incorrect format - must be a single logical value")
  if (length(add.legend) != 1) stop ("'add.legend' has incorrect format - must be a single logical value")

  # check segment.color
  if (!is.null(segment.color)) {
    if (!is.vector(segment.color)) stop ("'segment.color' has incorrect format - must be a vector")
    if (!is.character(segment.color) & !is.numeric(segment.color)) stop ("'segment.color' must be of type 'character' or 'numeric'")
    if (length(segment.color) < 2) stop ("'segment.color' has incorrect format - must specify 2 colors")
    if (length(segment.color) > 2) {
      warning ("'segment.color' has ",length(segment.color)," colors specified and requires 2. Using first 2 colors only")
      segment.color <- segment.color[1:2]
    }
    if (!all(areColors(segment.color))) stop ("some 'segment.color' are not valid colors")
  } else
    segment.color <- c("#69B4FF","#99DD55")

  # check chromosomes in data.frame
  if (is.null(interval)) {
    chromosomes <- as.character(unique(genotypes[,"chr"]))
    chromosomes.ibd <- as.character(unique(ibd.segments[,"chr"]))
    if (!all(chromosomes.ibd %in% chromosomes)) stop ("different chromosomes in 'ibd.segments' and 'ped.genotypes'")
  } else {
    if (!(interval.chr %in% unique(genotypes[,"chr"])))
      stop("'interval' chromosome is not in 'ped.genotypes'")
    if (!(interval.chr %in% unique(ibd.segments[,"chr"])))
      stop("'interval' chromosome is not in 'ibd.segments'")
    chromosomes <- interval.chr
  }

  # subset IBD segments by interval
  if (!is.null(interval)) {
    ibd.chr <- ibd.segments[ibd.segments[,"chr"] == interval.chr,]
    ibd.overlap <- rep(0, nrow(ibd.chr))
    for(i in 1:nrow(ibd.chr)){
      my.overlap <- getOverlap(ibd.chr[i,c("start_position_bp","end_position_bp")], c(interval.start, interval.stop))
      if(my.overlap[2] > my.overlap[1])
        ibd.overlap[i] <- 1
    }
    ibd.interval <- ibd.chr[ibd.overlap == 1,]
    if (nrow(ibd.interval) == 0) stop ("no IBD segments detected over interval")
  } else
    ibd.interval <- ibd.segments

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
    newpos   <- NULL
    chradd   <- NULL
    labelpos <- NULL
    # change genomic positions for each
    for (i in 1:length(chromosomes)) {
      maxpos <- max(genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"])
      minpos <- min(genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"])
      newpos <- c(newpos, genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"] + chrstart)
      ibd.interval[ibd.interval[,"chr"] == chromosomes[i],"start_position_bp"] <- ibd.interval[ibd.interval[,"chr"] == chromosomes[i],"start_position_bp"] + chrstart
      ibd.interval[ibd.interval[,"chr"] == chromosomes[i],"end_position_bp"] <- ibd.interval[ibd.interval[,"chr"] == chromosomes[i],"end_position_bp"] + chrstart
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
    newpos <- genotypes[genotypes[,"chr"] == interval.chr & genotypes[,"pos_bp"] >= interval.start & genotypes[,"pos_bp"] <= interval.stop,"pos_bp"]
    if (length(newpos) == 0) warning ("no SNPs in 'interval'")
  }

  # for each unique group pair, get numeric values for pairs
  ibd.interval$ibd.pairs <- paste(ibd.interval[,"fid1"],ibd.interval[,"iid1"],ibd.interval[,"fid2"],ibd.interval[,"iid2"],sep="/")
  unique.pairs <- unique(ibd.interval[,"ibd.pairs"])
  pair.id <- 1:length(unique.pairs)
  unique.pair.id <- data.frame(unique.pairs, pair.id)
  ibd.interval.2 <- merge(ibd.interval,unique.pair.id,by.x="ibd.pairs",by.y="unique.pairs")

  # set number per page
  if (!is.null(number.per.page)) {
    page.start <- 1
    for (i in seq(1, length(unique.pairs), number.per.page)) {
      if (i == max(seq(1, length(unique.pairs), number.per.page))) {
        page.pairs <- unique.pairs[page.start:length(unique.pairs)]
      } else
        page.pairs <- unique.pairs[page.start:(page.start+number.per.page-1)]
      ibd.interval.2[ibd.interval.2[,"ibd.pairs"] %in% page.pairs, "page.num"] <- i
      for (j in 1:length(page.pairs)) {
        ibd.interval.2[ibd.interval.2[,"ibd.pairs"] == page.pairs[j], "pair.id"] <- j
      }
      page.start <- page.start + number.per.page
    }
  } else {
    ibd.interval.2[,"page.num"] <- 1
  }

  # plotting segments:

  for (i in unique(ibd.interval.2[,"page.num"])) {
    ibd.page <- ibd.interval.2[ibd.interval.2[,"page.num"] == i,]

    # base plot
    ggp <- ggplot()
    ggp <- ggp + theme_bw()
    ggp <- ggp + theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       strip.background = element_rect(fill = "white",color = "white"),
                       legend.title = element_blank(),
                       strip.text.y = element_text(angle = 360),
                       axis.title.y = element_blank(),
                       axis.text.y = element_text(size=ylabel.size))
    if (add.legend) {
      ggp <- ggp + geom_rect(data=ibd.page, aes_(xmin = ~start_position_bp, xmax = ~end_position_bp, ymin = ~pair.id, ymax = ~pair.id+segment.height,
                                                fill = ~as.factor(ibd_status), alpha=0.8))
      ggp <- ggp + scale_fill_manual("", values = c(segment.color[1], segment.color[2]), labels=c("IBD = 1", "IBD = 2"))
    } else
      ggp <- ggp + geom_rect(data=ibd.page, aes_(xmin = ~start_position_bp, xmax = ~end_position_bp, ymin = ~pair.id, ymax = ~pair.id+segment.height),
                                            fill = ifelse(ibd.page[,"ibd_status"] == 1, segment.color[1], segment.color[2]), alpha=0.8)
    if (add.rug & length(newpos) != 0)
      ggp <- ggp + geom_rug(aes(x = newpos), size = 0.1, colour = "gray30")
    if (!is.null(plot.title))
      ggp <- ggp + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))

    # interval:
    if (is.null(interval) & length(chromosomes) > 1) {
      ggp <- ggp + xlab("Chromosome")
      ggp <- ggp + geom_vline(xintercept = chradd, colour = "gray87", linetype = "longdash", size = 0.4)
      ggp <- ggp + scale_x_continuous(breaks = labelpos, labels = chromosomes)
      ggp <- ggp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    } else {
      ggp <- ggp + xlab(paste("Chromosome", chromosomes))
      if (!is.null(interval))
        ggp <- ggp + coord_cartesian(xlim = c(interval.start, interval.stop))
    }

    # annotation.genes:
    if (!is.null(annotation.genes)) {
      if(nrow(annotation.genes.overlap) != 0) {
        #min.y <- -0.05*max(ibd.page[,"pair.id"])
        #max.y <- -0.01*max(ibd.page[,"pair.id"])
        gene.hight <- 0.05*(max(ibd.page[,"pair.id"]) - min(ibd.page[,"pair.id"]))
        max.y <- min(ibd.page[,"pair.id"]) - gene.hight*0.5
        min.y <- max.y - gene.hight
        pos.strand <- annotation.genes.overlap[annotation.genes.overlap[,"strand"] == "+",]
        neg.strand <- annotation.genes.overlap[annotation.genes.overlap[,"strand"] != "+",]
        ggp <- ggp + geom_rect(data=pos.strand, aes_string(xmin = "start", xmax = "end"), ymin = min.y, ymax = max.y, alpha = 0.9, fill = annotation.genes.color[1])
        ggp <- ggp + geom_rect(data=neg.strand, aes_string(xmin = "start", xmax = "end"), ymin = min.y, ymax = max.y, alpha = 0.9, fill = annotation.genes.color[2])
        #ggp <- ggp + ylim(min.y, max(ibd.page[,"pair.id"])) # overrides facets="free"
      }
    }

    # highlight.genes:
    if (!is.null(highlight.genes)) {
      if(nrow(highlight.genes.overlap) != 0) {
        lab.pos <- min(ibd.page[,"pair.id"])*0.05
        ggp <- ggp + geom_rect(data=highlight.genes.overlap, aes_string(xmin = "start", xmax = "end"), ymin = -Inf, ymax = Inf, fill = highlight.genes.color, alpha = highlight.genes.alpha)
        ggp <- ggp + geom_vline(data=highlight.genes.overlap, aes_string(xintercept = "start"), colour = highlight.genes.color, linetype = "solid", alpha = highlight.genes.alpha)
        if (highlight.genes.labels)
          ggp <- ggp + geom_text(data=highlight.genes.overlap, aes_string(x = "start", label = "name"), y = lab.pos, colour = "gray20", angle = 90, hjust = -0.1, vjust = -0.2, size = 3, alpha = 0.6)
      }
    }

    # y-axis labels & limits
    y.labs <- ibd.page[!duplicated(ibd.page[,"pair.id"]),]
    y.breaks <- y.labs[,"pair.id"] + segment.height/2
    if (fid.label & iid.label){
      y.labels <- paste(y.labs[,"fid1"], y.labs[,"iid1"], y.labs[,"fid2"], y.labs[,"iid2"],sep="/")
    }
    if (fid.label & !iid.label){
      y.labels <- paste(y.labs[,"fid1"], y.labs[,"fid2"],sep="/")
    }
    if (!fid.label & iid.label){
      y.labels <- paste(y.labs[,"iid1"], y.labs[,"iid2"],sep="/")
    }
    if (!fid.label & !iid.label){
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.overlap) > 0)
          ggp <- ggp + ylim(min.y, max(ibd.page[,"pair.id"]))
      }
      ggp <- ggp + theme(axis.text.y=element_blank(), axis.ticks=element_blank())
    } else {
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.overlap) > 0) {
          ggp <- ggp + scale_y_continuous(breaks = y.breaks, labels = y.labels, limits=c(min.y, max(ibd.page[,"pair.id"])))
        } else
          ggp <- ggp + scale_y_continuous(breaks = y.breaks, labels = y.labels)
      } else
        ggp <- ggp + scale_y_continuous(breaks = y.breaks, labels = y.labels)
    }


    print(ggp)
  }

}
