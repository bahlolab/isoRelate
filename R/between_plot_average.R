#' Between Isolates Average Sharing Plot
#' 
#' \code{plotIBDaverage()} plots the proportion of pairs IBD for each SNP
#' @param locus.matrix A list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @export
plotIBDaverage <- function(ped.genotypes, ibd.segments, groups = NULL, number.cores = 1){
  
  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 2)
  pedigree  <- ped.genotypes[[1]] 
  genotypes <- ped.genotypes[[2]]
  
  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("ped.genotypes has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "moi", "aff")
  #pedigree <- checkClass(pedigree, c(rep("c",4),rep("n",2)))
  
  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) < 8 & nrow(genotypes) <= 1)
    stop ("ped.genotypes has incorrect format")
  colnames(genotypes)[1:5] <- c("CHROMOSOME", "MARKER", "POSITION.M","POSITION.bp", "FREQ")
  
  # check ibd.segments file is a dataframe with correct fields
  stopifnot(is.data.frame(ibd.segments))
  if (ncol(ibd.segments) != 15) 
    stop ("ibd.segments has incorrect format")
  colnames(ibd.segments) <- c("fid1", "ind1", "fid2", "ind2", "chr", "start.snp", "end.snp", "start.position.bp", 
                              "end.position.bp", "start.position.M", "end.position.M", "number.snps", "length.bp",
                              "length.M", "ibd.status")
  #ibd.segments <- checkClass(ibd.segments, c(rep("c",7),rep("n",8)))
  
  # check groups
  if (!is.null(groups)) {
    stopifnot(is.data.frame(groups))
    stopifnot(ncol(groups) > 2)
    if (ncol(groups) > 4){
      cat("using first 4 columns of groups")
      groups <- groups[,1:4]
    }
    colnames(groups)[1:2] <- c("fid","iid")
    #groups <- checkClass(groups, c("c","c","f"))
  }
    
  # check numeric values
  stopifnot(is.numeric(number.cores))
    
  # chromosome names
  chromosomes <- as.character(unique(genotypes[,"CHROMOSOME"]))
  
  # genome length
  genome.length <- 0
  for(chrom in chromosomes){
    positions.chrom <- genotypes[genotypes[,"CHROMOSOME"] == chrom,"POSITION.bp"]
    genome.length <- genome.length + (max(positions.chrom) - min(positions.chrom))
  }
  
  # create new isolate IDs from PED FIDs and IIDs
  isolate.pairs <- isolatePairs(pedigree[,1], pedigree[,2])
  isolate.pairs <- paste(isolate.pairs[,1],isolate.pairs[,2],isolate.pairs[,3],isolate.pairs[,4],sep="/")
  
  # check the isolates belong to a group
  if (!is.null(groups)) {
    group.names <- paste(groups[,"fid"],groups[,"iid"],sep="/")
    isolate.names <- paste(pedigree[,"fid"],pedigree[,"iid"],sep="/")
    if (!all(group.names %in% isolate.names)) 
      stop("'groups' is missing information for some isoaltes")
    
    # assign number ID to each isoalte
    pedigree.0 <- data.frame(num.id=1:nrow(pedigree), pedigree)
    
    # merge pedigree with proups by IDs
    pedigree.group <- merge(pedigree.0, groups, by=c("fid", "iid"))
    
    # reorder merged pedigree by numberic IDs
    pedigree.group <- pedigree.group[order(pedigree.group[,"num.id"]),]
    
    # get isolates group pairs - dataframe with 2 coloumns
    group.pairs <- groupPairs(as.character(pedigree.group[,8]))
    
    # reorder pairs groups
    groups.unique <- as.character(unique(pedigree.group[,8]))
    groups.unique.pairs <- groupPairs(groups.unique)
    if (nrow(groups.unique.pairs) > 0) {
      for (i in 1:nrow(groups.unique.pairs)){
        group.pairs[group.pairs[,1] == groups.unique.pairs[i,2] & group.pairs[,2] == groups.unique.pairs[i,1],] <- groups.unique.pairs[i,]
      }
    }
    
    # combine group pairs into a single vector
    group.pairs.1 <- paste(group.pairs[,1],group.pairs[,2],sep="\n")
    group.pairs.1 <- as.factor(group.pairs.1)
  } else {
    group.pairs.1 <- as.factor(rep("NA",length(isolate.pairs)))
  }
  
  # create a dataframe with isolate pairs, group ids, and average value = 0
  ave.pairs <- data.frame(isolate.pairs = isolate.pairs, group.pairs = group.pairs, number.ibd = 0, 
                          ave.length.bp = 0, ave.length.M = 0, ave.proportion = 0)
  pair.id <- 1:nrow(ave.pairs)
  
  # create pairs from ibd.segments
  ibd.pairs <- paste(ibd.segments[,"fid1"],ibd.segments[,"ind1"],ibd.segments[,"fid2"],ibd.segments[,"ind2"],sep="/")
  
  # define number of cores
  doParallel::registerDoParallel(cores=number.cores)
  
  # for each pair with ibd inferred, calculate the number of ibd segments and the size of them
  average.ibd <- foreach::foreach(pair=unique(ibd.pairs), .combine='rbind') %dopar% {
    pair.ibd.segments <- ibd.segments[ibd.pairs == pair,]
    number.ibd <- nrow(pair.ibd.segments)
    ave.length.bp <- mean(pair.ibd.segments[,"length.bp"])
    ave.length.M  <- round(mean(pair.ibd.segments[,"length.M"]),2)
    ave.proportion <- sum(pair.ibd.segments[,"length.bp"])/genome.length
    cbind(pair, number.ibd, ave.length.bp, ave.length.M,ave.proportion)
  }
  ave.pairs[pair.id[ave.pairs[,"isolate.pairs"] %in% unique(ibd.pairs)],c("isolate.pairs",
    "number.ibd", "ave.length.bp", "ave.length.M", "ave.proportion")] <- average.ibd
  
  # plot pair colours by groups
  par(mfrow=c(2,1))
  if (!is.null(groups)) {
    # for each original group, choose a colour
    new.group.pairs <- unique(paste(group.pairs[,1],group.pairs[,2],sep="/"))
    number.groups <- length(new.group.pairs)

    # outline and fill colour of boxplots
    out.colour  <- getColourPaletteMajor(number.groups)
    fill.colour <- NULL
    for (i in 1:number.groups) {
      fill.colour.pallete <- colorRampPalette(c(out.colour[i], "white"))
      fill.colour[i] <- fill.colour.pallete(3)[2]
    }
    
    # plot average kb
    ave.pairs.kb <- as.numeric(ave.pairs[,"ave.length.bp"])/1000
    boxplot(ave.pairs.kb~group.pairs.1,col=fill.colour,border=out.colour,axes=FALSE)
    grid(nx=NA, ny=NULL) 
    par(new=TRUE)
    boxplot(ave.pairs.kb~group.pairs.1,col=fill.colour,border=out.colour,las=2, ylab="Average Length (kb)")
    
    # plot number of IBD segments
    number.ibd.0 <- as.numeric(ave.pairs[,"number.ibd"])
    boxplot(number.ibd.0~group.pairs.1,col=fill.colour,border=out.colour,axes=FALSE)
    grid(nx=NA, ny=NULL) 
    par(new=TRUE)
    boxplot(number.ibd.0~group.pairs.1,col=fill.colour,border=out.colour,las=2, ylab="Number of IBD Segments")

  } else {
    
    # plot average kb
    ave.pairs.kb <- as.numeric(ave.pairs[,"ave.length.bp"])/1000
    hist(ave.pairs.kb,axes=FALSE,xlab="",ylab="",main="")
    grid(nx=NULL, ny=NULL) 
    par(new=TRUE)
    hist(ave.pairs.kb,las=2, xlab="Average Length (kb)",main="")
    
    # plot number of IBD segments
    number.ibd.0 <- as.numeric(ave.pairs[,"number.ibd"])
    hist(number.ibd.0,axes=FALSE,xlab="",ylab="",main="")
    grid(nx=NULL, ny=NULL) 
    par(new=TRUE)
    hist(number.ibd.0, xlab="Number of IBD Segments",main="")
  }
  #return(ave.pairs)
}
