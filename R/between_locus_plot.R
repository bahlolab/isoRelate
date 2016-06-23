#' Between Isolates Loci Sum Plot
#' 
#' \code{plotIBDglobal()} plots the proportion of pairs IBD for each SNP
#' @param locus.proportions a data frame containing the binary IBD information for each SNP and each pair. 
#' See the returned \code{value} in \code{\link{between_locus_matrix}} for more details.
#' @param interval a vector of length 3 containing the region to identify clusters over. This vector should contain the
#' chromosome ID, the start of the interval in bp and the end of the interval in bp; in this order respectively. Note that 
#' \code{snpName} and \code{interval} are mutually exclusive and only one option should be specified.
#' @param genes a data frame with at least 5 columns of information:
#' \enumerate{
#' \item chr
#' \item name
#' \item start
#' \item end
#' \item strand
#' }
#' This data frame does not have to be in a specific order, however it must contain all of the above information with identical labels.
#' @importFrom ggplot2 ggplot
#' @export
plotIBDglobal <- function(locus.proportions, interval = NULL, genes = NULL){

  # check locus matrix input
  if (ncol(locus.proportions) < 5)
    stop ("locus.proportions has incorrect format")
  colnames(locus.proportions)[1:4] <- c("CHROMOSOME","MARKER","POSITION.M","POSITION.bp")
  
  # if interval specified
  if (!is.null(interval) & length(interval) == 3) {
    chr   <- as.character(interval[1])
    start <- as.numeric(interval[2])
    stop  <- as.numeric(interval[3])
    stopifnot(chr %in% locus.proportions[,"CHROMOSOME"])
    if(start > stop) 
      stop(paste("interval start=",interval[2]," is greater than interval end=",interval[3],sep=""))
    locus.interval <- locus.proportions[locus.proportions[,"CHROMOSOME"] == chr & locus.proportions[,"POSITION.bp"] >= start &
                                     locus.proportions[,"POSITION.bp"] <= stop,]
    if(nrow(locus.interval) == 0)
      stop("no SNPs in interval")
  } else {
    locus.interval <- locus.proportions
  }
  
  # check genes
  if (!is.null(genes)) {
    stopifnot(is.data.frame(genes))
    stopifnot(c("name","strand","chr","start","end") %in% colnames(genes))
    
    # subset genes if interval specified
    if (!is.null(interval)) {
      genes.0 <- genes[genes[,"chr"] == chr,]
      # function to find genes that overlap interval
      fun.overlap <- function(region.1, region.2) {
        a <- max(region.1[1], region.2[1])
        b <- min(region.1[2], region.2[2])
        return(c(a,b))
      }
      genes.1 <- NULL
      for(g in 1:nrow(genes.0)){
        gene.overlap <- fun.overlap(genes.0[g,c("start","end")],c(start,stop))
        if (gene.overlap[2] - gene.overlap[1] > 0)
          genes.1 <- rbind(genes.1, genes.0[g,])
      }
      genes.1 <- data.frame(genes.1)
    } else {
      genes.1 <- genes
    }
    if(nrow(genes.1) == 0) {
      cat("no genes in interval")
    } else {
      genes.1[,"start"] <- as.numeric(genes.1[,"start"])
      genes.1[,"end"] <- as.numeric(genes.1[,"end"])
    }
  }
  
  # chromosome names
  chromosomes <- as.character(unique(locus.interval[,"CHROMOSOME"]))
  
  # genome length
  genome.length <- 0
  for(chrom in chromosomes){
    positions.chrom <- locus.interval[locus.interval[,"CHROMOSOME"] == chrom,"POSITION.bp"]
    genome.length <- genome.length + (max(positions.chrom) - min(positions.chrom))
  }
  
  # define new plot positions
  if (is.null(interval)) {
    chrstart <- 0
    chradd   <- 0
    labelpos <- NULL
    newpos   <- NULL
    for(i in 1:length(chromosomes)){
      maxpos   <- max(locus.interval[locus.interval[,"CHROMOSOME"] == chromosomes[i],"POSITION.bp"])
      minpos   <- min(locus.interval[locus.interval[,"CHROMOSOME"] == chromosomes[i],"POSITION.bp"])
      newpos   <- c(newpos, locus.interval[locus.interval[,"CHROMOSOME"] == chromosomes[i],"POSITION.bp"] + chrstart)
      labelpos[i] <- (maxpos - minpos + 2*chrstart)/2
      chradd[i]   <- chrstart
      if (!is.null(genes)) {
        if (nrow(genes.1) != 0){
          genes.1[genes.1[,"chr"] == chromosomes[i],"start"] <- genes.1[genes.1[,"chr"] == chromosomes[i],"start"] + chrstart
          genes.1[genes.1[,"chr"] == chromosomes[i],"end"] <- genes.1[genes.1[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      chrstart <- chrstart + maxpos
    }
    chradd <- chradd[-1]
  }
  
  # create dataframes for plotting
  if (is.null(interval)) {
    locus.df <- data.frame(pos=newpos,locus.proportions[,5:ncol(locus.proportions)])
    #colnames(locus.df) <- c("pos", colnames(locus.prop))
    locus.df.melt <- melt(locus.df,id="pos")
  } else {
    locus.df <- data.frame(pos=locus.interval[,"POSITION.bp"],locus.proportions[,5])
    #colnames(locus.df) <- c("pos", colnames(locus.prop))
    locus.df.melt <- melt(locus.df,id="pos")
  }
  number.groups <- ncol(locus.df) - 1
  
  # calculate y-axis limit - -10% of max prop if reference genes supplied
  if (!is.null(genes)) {
    y.min <- -max(locus.proportions[,5:ncol(locus.proportions)])*0.1
  } else
    y.min <- 0
  y.max <- max(locus.proportions[,5:ncol(locus.proportions)])
  
  # setting up ggplot
  # groups:
  if (number.groups == 1) {
    ggp <- ggplot2::ggplot(locus.df.melt, ggplot2::aes(pos, value))
  } else {
    ggp <- ggplot2::ggplot(locus.df.melt, ggplot2::aes(pos, value, col=variable)) +
      # separate plots by groups
      ggplot2::facet_grid(variable~.) +
      # get group colours
      ggplot2::scale_colour_manual(values=getColourPaletteMajor(number.groups))
  }
  
  # parameters in all plots
  ggp <- ggp +
    # line plot
    ggplot2::geom_line() + 
    # y axis limits
    ggplot2::ylim(y.min,y.max) +
    # remove gray colour from background
    ggplot2::theme_bw() + 
    # remove grid lines
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major=ggplot2::element_blank()) +
    # remove box around facet (group labels) and transform
    ggplot2::theme(strip.background = ggplot2::element_rect(fill="white",color="white"),
          strip.text.y = ggplot2::element_text(angle=360)) +
    # y axis labels
    ggplot2::ylab("Proportion of Pairs IBD") + 
    # remove legend
    ggplot2::theme(legend.position = "none")
  
  # interval:
  if (is.null(interval)) {
    # x label
    ggp <- ggp + ggplot2::xlab("Chromosome") +
      # add vertical lines to separate chromosomes
      ggplot2::geom_vline(xintercept = chradd, colour="gray87", linetype = "longdash") +
      # change x axis text and transform
      ggplot2::scale_x_continuous(breaks=labelpos, labels=chromosomes) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,vjust=0.5))
  } else {
    # x label
    ggp <- ggp + ggplot2::xlab(paste("Chromosome",chromosomes)) +
      ggplot2::geom_rug(sides="b", color="black")
  }
  
  # genes:
  if (!is.null(genes)) {
    if(nrow(genes.1) != 0) {
      genes.new <- data.frame(genes.1,y1=y.min,y2=0)
      ggp <- ggp + ggplot2::geom_rect(data=genes.new, ggplot2::aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=strand), inherit.aes=FALSE)
      #ggp <- ggp + geom_rect(data=genes.new, aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=ifelse(strand=="+", "black", "red")), inherit.aes=FALSE)
    }
  }
  
  plot(ggp)
}