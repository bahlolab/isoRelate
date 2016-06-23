#' Internal Function
#' 
#' \code{overlappingSegment()} returns the maximum start of two regions and the minimum end of two regions.
#' @param a A vector of containing the start and end positions of an interval in bp
#' @param b A vector of containing the start and end positions of an interval in bp

overlappingSegment <- function(a, b){
  x <- max(a[1], b[1])
  y <- min(a[2], b[2])
  return(c(x,y))
}

#' Internal Function
#' 
#' \code{plotIBDzoom()} plots IBD segments in user specified regions.
#' @param ibd.segments A data frame with inferred IBD tracts. See the \code{Value} description in \code{\link{between_ibd_segments}} 
#' for more details 
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{between_processing}} for more details on this input.
#' @param ibd.colours A vector with two values. These can be either numeric or the name of colours to be used in the IBD plot (must be vaild colour
#' names).
#' @param pairs.per.page Integer. The maximum number of pairwise results to plot on each page of the PDF.
#' @param iid.only Use only the isolate IDs on the figure legend of the plot. If using, the isolate IDs must be unique for each isolate.
#' \code{iid.only} and \code{fid.only} are mutually exclusive.
#' @param fid.only Use only the family IDs on the figure legend of the plot. If using, the family IDs must be unique for each isolate.
#' @param zoom A data frame containing three columns:
#' \enumerate{
#' \item Chromosome
#' \item Start location bp
#' \item End location bp
#' }
#' More than one region can be specified per chromosome however they will be plotted in spearate figures.

plotIBDzoom <- function(ibd.segments, ped.genotypes, ibd.colours = c("dodgerblue","darkgoldenrod1"), 
                          pairs.per.page = 50, iid.only = FALSE, fid.only = FALSE, zoom){
  
  # subset genotype data
  genotypes <- ped.genotypes[[2]]
  
  for (i in 1:nrow(zoom)) {
    # subset IBD segments by chromosome
    ibdSegments_chr <- ibd.segments[ibd.segments[,"chr"] == zoom[i,1],]
    
    # find IBD segments that overlap zoomed region
    ibdSegments_chr$overlap0 <- 0
    for (j in 1:nrow(ibdSegments_chr)) {
      overlap0 <- overlappingSegment(ibdSegments_chr[j,c("start.position.bp","end.position.bp")],zoom[i,2:3])
      if (overlap0[2] - overlap0[1] > 0) {
        ibdSegments_chr$overlap0[j] <- 1 
        if (ibdSegments_chr$start.position.bp[j] < zoom[i,2]) 
          ibdSegments_chr$start.position.bp[j] = zoom[i,2]
        if (ibdSegments_chr$end.position.bp[j] > zoom[i,3])
          ibdSegments_chr$end.position.bp[j] = zoom[i,3]
      }
    }
    
    CH     <- ibdSegments_chr[ibdSegments_chr[,"overlap0"] == 1,]
    SNPpos <- genotypes[genotypes[,"CHROMOSOME"] == zoom[i,1], "POSITION.bp"]
    
    if (nrow(CH) > 0) {
      # y-axis labels
      if (iid.only == TRUE) {
        CH$SampleID <- paste(CH$ind1, CH$ind2,sep="  /  ")
      } else if (fid.only == TRUE) {
        CH$SampleID <- paste(CH$fid1, CH$fid2,sep="  /  ")
      } else {
        CH$Fam1ID   <- paste(CH$fid1, CH$ind1,sep="_")
        CH$Fam2ID   <- paste(CH$fid2, CH$ind2,sep="_")
        CH$SampleID <- paste(CH$Fam1ID, CH$Fam2ID,sep="  /  ")
      }
      
      # create pair IDs
      Samples   		       <- as.data.frame(unique(CH$SampleID))
      colnames(Samples)[1] <- "SampleID"
      Samples$index        <- seq(1,nrow(Samples))
      CH 			             <- merge(CH,Samples,by.x="SampleID",by.y="SampleID")
      
      for (k in seq(1, nrow(Samples), pairs.per.page)) {
        # select samples to plot
        StartSample <- k
        EndSample   <- min(k+pairs.per.page-1,nrow(Samples))
        ch          <- CH[CH$index>=StartSample & CH$index <= EndSample,]
        cat("plotting Chr ",zoom[i,1],",",nrow(Samples),"samples, sample",StartSample,"to sample", EndSample, " ... ")
        
        # assign IBD segment colour
        for (j in 1:nrow(ch)) {
          if (ch$ibd_status[j]==1){ ch$Col[j] <- ibd.colours[1] } else ch$Col[j] <- ibd.colours[2]
        }
        
        # define plot
        par(xaxt = "s", yaxt = "n",las=0, cex.axis=1,mar=c(4,8,2,1))
        plot(c(ch$start.position.bp,ch$end.position.bp),c(ch$index,ch$index), type = "n", xlab = "",ylab="",tcl=-0.8,
             ylim=c(min(ch$index-0.3),max(ch$index+0.3)),xlim=c(zoom[i,2],zoom[i,3]))
        title(expression("IBD = 1" * phantom(" and IBD = 2")),col.main=ibd.colours[1])
        title(expression(phantom("IBD = 1 and ") * "IBD = 2"),col.main=ibd.colours[2])
        title(expression(phantom("IBD = 1 ") * "and " * phantom("IBD = 2"),col.main="black")) 
        
        # plotting IBD segments
        rect(ch$start.position.bp,ch$index-0.3,ch$end.position.bp,ch$index+0.3,border=ch$Col,col=ch$Col,density=30)
        par(las=2,yaxt="s",cex.axis=0.6)
        samples=unique(data.frame(ch$SampleID,ch$index))      
        axis(2, at = samples$ch.index, labels =  samples$ch.SampleID) #Y label for sample ID
        mtext(paste("Chromosome ",zoom[i,1],sep=""),1, line = 3, las = 0)
        axis(1, at=SNPpos,tick=TRUE,labels=FALSE)
        box()
        cat("\n")
      }
    } else cat(paste("No IBD segments overlapping chromosome=",zoom[i,1],":",zoom[i,2],"-",zoom[i,3],sep=""))
  }
}
