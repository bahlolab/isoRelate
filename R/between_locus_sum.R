#' Between Isolates Proportion of Pairs IBD
#' 
#' \code{getLocusProportion()} calculates the number of IBD segments inferred at each SNP over all pairs.
#' @param ped.genotypes a list containing 2 objects. See the \code{Value} description in \code{\link{between_data_processing}} for more details on this input.
#' Note the family IDs and isolate IDs in obeject 1 of this list must match the family IDs and isolate IDs in the header of object 2 of this list.
#' @param locus_matrix A data frame containing the binary IBD information for each SNP and each pair. 
#' See the returned \code{value} in \code{\link{between_locus_matrix}} for more details.
#' @param groups a data frame with a least 3 coloumns of information:
#' \enumerate{
#' \item Family ID
#' \item Isolate ID
#' \item Group ID 1
#' }
#' Group ID, for example, can be geographic regions where the isolates were collected. A fourth column can also be included with a second level of group IDs.
#' If \code{groups} is specified, each isolate in the pedigree should belong to a group. If more than two levels of groups are specified, only the first
#' two levels will be used, i.e. columns 1 to 4 in \code{groups} will be used. The default is \code{groups=NULL} for no groups.
#' @return A data frame the following 5 columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' \item Number of IBD segments (type \code{"integer"})
#' }
#' where each row describes a unique SNP. The data frame is headed 
#' \code{CHROMOSOME, MARKER, POSITION.cM, POSITION.bp} and \code{prop} respectively.
#' @export
getLocusProportion <- function(ped.genotypes, locus.matrix, groups = NULL){
  
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
  
  # check groups
  if (!is.null(groups)) {
    stopifnot(is.data.frame(groups))
    stopifnot(ncol(groups) > 2)
    if (ncol(groups) > 4){
      cat("using first 4 columns of groups")
      groups <- groups[,1:4]
    }
    colnames(groups)[1:2] <- c("fid","iid")
    
    # check isolates belong to a group
    group.names <- paste(groups[,"fid"],groups[,"iid"],sep="/")
    isolate.names <- paste(pedigree[,"fid"],pedigree[,"iid"],sep="/")
    if (!all(isolate.names %in% group.names)) 
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
    
    # get number of pairwise analyses for each group
    number.pairs.groups <- NULL
    for (i in unique(paste(group.pairs[,1],group.pairs[,2],sep="/"))) {
      id_1 <- unlist(strsplit(i,"/"))[1]
      id_2 <- unlist(strsplit(i,"/"))[2]
      number.pairs.groups[group.pairs[,1] == id_1 & group.pairs[,2] == id_2] <- paste("npairs=",
                                                                                      dim(group.pairs[group.pairs[,1] == id_1 & group.pairs[,2] == id_2,])[1],
                                                                                      sep="")
    }
    
    # combine group pairs into a single vector
    group.pairs.1 <- paste(group.pairs[,1],group.pairs[,2],number.pairs.groups,sep="\n")
  }
  
  # calculate proportion IBD at each SNP
  locus.pairs <- locus.matrix[,5:ncol(locus.matrix)]
  if (!is.null(groups)) {
    locus.prop <- NULL
    for (g in unique(group.pairs.1)) {
      locus.pairs.g <- locus.pairs[,group.pairs.1 == g]
      locus.prop <- cbind(locus.prop, rowMeans(locus.pairs.g))
    }
    colnames(locus.prop) <- unique(group.pairs.1)
  } else {
    locus.prop <- rowMeans(locus.pairs)
  }
  
  return.locus.prop <- cbind(locus.matrix[,1:4],locus.prop)
  if(is.null(groups))
    colnames(return.locus.prop) <- "prop"
  
  return(data.frame(return.locus.prop))
}


