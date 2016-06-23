#' Internal Function
#' 
#' \code{plotIBD()} plots all IBD segments for all chromosomes in one figure. A maximum of 50 pairs are plotted in one figure and remaining
#' pairs are plotted on subsequent pdf pages with a maximum of 50 pairs in one figure.
#' @param ibd.segments A data frame with inferred IBD tracts. See the \code{Value} description in \code{\link{between_ibd_segments}} 
#' for more details 
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_processing}} for more details on this input.
#' @param chrom A vector of chromosomes to plot. Can be either numeric or character.
#' @param ibd.colours A vector with two values. These can be either numeric or the name of colours to be used in the IBD plot (must be vaild colour
#' names).
#' @param pairs.per.page Integer. The maximum number of pairwise results to plot on each page of the PDF.
#' @param iid.only Use only the isolate IDs on the figure legend of the plot. If using, the isolate IDs must be unique for each isolate.
#' \code{iid.only} and \code{fid.only} are mutually exclusive.
#' @param fid.only Use only the family IDs on the figure legend of the plot. If using, the family IDs must be unique for each isolate.

plotIBD <- function(ibd.segments, ped.genotypes, chrom = NULL, ibd.colours = c("dodgerblue","darkgoldenrod1"), 
                     pairs.per.page = 50, iid.only = FALSE, fid.only = FALSE){
  
  # subset genotype data
  genotypes <- ped.genotypes[[2]]
  
  # select chomosomes for plotting
  if(!is.null(chrom)){
    chromosomes <- unlist(strsplit(chrom,","))
  } else {
    chromosomes <- unique(as.character(genotypes[,"CHROMOSOME"]))
  }
  
  # find the length of each chromosome
  chrom.length <- NULL
  for(i in chromosomes){
    chrom.length <- c(chrom.length, max(genotypes[genotypes[,"CHROMOSOME"] == i, "POSITION.bp"]))
  }
  names(chrom.length) <- chromosomes
  
  # create a pseudo position for each SNP
  newCH  <- NULL
  SNPpos <- NULL
  allPos <- NULL
  newPos <- 0
  xaxlab <- NULL
  xaxPos <- NULL
  for (i in chromosomes) {
    CH=ibd.segments[ibd.segments$chr==i,]
    if (nrow(CH) > 0) {
      CH[,"start.position.bp"] <- CH[,"start.position.bp"] + newPos
      CH[,"end.position.bp"]   <- CH[,"end.position.bp"] + newPos
      newCH <- rbind(newCH, CH)
      xaxlab <- c(xaxlab, newPos + chrom.length[names(chrom.length)==i]/2)
      xaxPos <- c(xaxPos, i)
      SNPpos <- c(SNPpos, genotypes[genotypes[,"CHROMOSOME"] == i,"POSITION.bp"] + newPos)
      newPos <- newPos + chrom.length[names(chrom.length)==i]
      allPos <- c(allPos, newPos)
    } 
  }
  
  if (nrow(newCH) > 0) {
    # y-axis labels
    if (iid.only == TRUE) {
      newCH$SampleID <- paste(newCH$ind1,newCH$ind2,sep="  /  ")
    } else if (fid.only == TRUE) {
      newCH$SampleID <- paste(newCH$fid1,newCH$fid2,sep="  /  ")
    } else {
      newCH$Fam1ID   <- paste(newCH$fid1,newCH$ind1,sep="_")
      newCH$Fam2ID   <- paste(newCH$fid2,newCH$ind2,sep="_")
      newCH$SampleID <- paste(newCH$Fam1ID,newCH$Fam2ID,sep="  /  ")
    }
    
    # create pair IDs
    Samples 			       <- as.data.frame(unique(newCH$SampleID))
    colnames(Samples)[1] <- "SampleID"
    #Samples$index       <- seq(1,nrow(Samples))
    Samples$index        <- seq(nrow(Samples),1)
    newCH 				       <- merge(newCH,Samples,by.x="SampleID",by.y="SampleID")
    ymar=8
    
    for (k in seq(1,nrow(Samples),pairs.per.page)) {
      # select samples to plot
      StartSample <- k
      EndSample   <- min(k+pairs.per.page-1,nrow(Samples))
      ch          <- newCH[newCH$index>=StartSample & newCH$index <= EndSample,]
      
      # assign IBD segment colour
      for (j in 1:nrow(ch)) {
        if (ch$ibd_status[j]==1) { 
          ch$Col[j] <- ibd.colours[1] 
        } else 
          ch$Col[j] <- ibd.colours[2]
      }
      
      # define plot
      par(yaxt = "n", xaxt="n",las=0, cex.axis=1,mar=c(2,ymar,2,2.1))
      plot(c(ch$start.position.bp,ch$end.position.bp),c(ch$index,ch$index), type = "n", xlab = "",ylab="",
           ylim=c(min(ch$index-0.3),max(ch$index+0.3)),xlim=c(0,newPos))
      title(expression("IBD = 1" * phantom(" and IBD = 2")),col.main=ibd.colours[1])
      title(expression(phantom("IBD = 1 and ") * "IBD = 2"),col.main=ibd.colours[2])
      title(expression(phantom("IBD = 1 ") * "and " * phantom("IBD = 2"),col.main="black")) 
      
      # plotting IBD segments
      rect(ch$start.position.bp,ch$index-0.3,ch$end.position.bp,ch$index+0.3,border=ch$Col,col=ch$Col)
      par(las=2,yaxt="s",cex.axis=0.5)
      samples=unique(data.frame(ch$SampleID,ch$index))      
      axis(2, at = samples$ch.index, labels =  samples$ch.SampleID) #Y label for sample ID
      par(las=0,xaxt="s",cex.axis=0.6)
      axis(1, at = xaxlab, labels = xaxPos,tcl=-0.01)
      abline(v=allPos[-length(allPos)],lty=2)
      axis(1, at=SNPpos,tick=TRUE,labels=FALSE)
      box()
      cat("plotting ",nrow(Samples),"samples, sample",StartSample,"to sample", EndSample, " ... \n")
    }
  }
}