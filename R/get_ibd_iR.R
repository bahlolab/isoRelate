#' Selection Significance Statistic
#'
#' \code{getIBDiR()} calculates a summary statistic for each SNP that can be used to assess the significance of excess IBD sharing at genomic loci,
#' thus identifying regions under positive selection.
#' First relatedness between isolates and SNP allele frequencies are accounted for, then normalization procedures are applied where we assume
#' our transformed summary statistic follows a chi-squared distribution with 1 degree of freedom.
#' This allows the calculation of –log10 (P-values) which we denote as the iR statistic.
#' SNPs with iR values greater than some threshold (i.e. –log10 (P-values) > –log10 (0.05)) provide evidence of positive selection.
#' \code{getIBDiR} can return \code{NA} iR statistics for a number of reasons, including trying to generate iR statistics when there are no IBD pairs
#' or when all pairs are IBD, or when only several isolates are analyzed.
#' @param ped.genotypes A list containing 2 objects. See the \code{Value} description in \code{\link{getGenotypes}} for more details on this input.
#' @param ibd.matrix A data frame containing the binary IBD information for each SNP and each pair.
#' See the returned \code{Value} in \code{\link{getIBDmatrix}} for more details.
#' @param groups A data frame with 3 columns of information:
#' \enumerate{
#' \item Family ID
#' \item Isolate ID
#' \item Group ID
#' }
#' where IBD proportions are calculated for
#' \enumerate{
#' \item all pairs of isolates within the same group
#' \item all pairwise-group comparisons where isolates belong to different groups
#' }
#' Group ID, for example, can be the geographic regions where the isolates were collected.
#' The default is \code{groups=NULL} and IBD proportions will be calculated over all pairs.
#' @return A data frame the following 7 columns:
#' \enumerate{
#' \item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
#' \item SNP identifiers (type \code{"character"})
#' \item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
#' \item Base-pair position (type \code{"integer"})
#' \item Population (type \code{"character"} or \code{"numeric"})
#' \item Subpopulation (type \code{"character"} or \code{"numeric"})
#' \item iR statistic (type \code{"numeric"})
#' }
#' where each row describes a unique SNP.
#' The column 'Population' is filled with ones by default, while 'Subpopulation' contains the group ID’s from  \code{groups},
#' where the proportion of pairs IBD has been calculated for all pairs of isolates belonging to the same group as well as all pairs of
#' isolates where each isolate belongs to a different group.
#' If \code{groups=NULL} then 'Subpopulation' will be filled with zeros also.
#' The population columns have been included for plotting purposes.
#' The data frame is headed \code{chr, snp_id, pos_M, pos_bp, pop, subpop} and \code{iR} respectively.
#' @export
getIBDiR <- function(ped.genotypes, ibd.matrix, groups = NULL){

  # check format of input data
  if (!is.list(ped.genotypes) | length(ped.genotypes) != 2) stop ("'ped.genotypes' must be a named list containing 2 objects: 'pedigree' and 'genotypes'")
  if (any(names(ped.genotypes) != c("pedigree", "genotypes"))) stop ("'ped.genotypes' must be a named list containing 'pedigree' and 'genotypes'")
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  if (!is.data.frame(pedigree)) stop ("'ped.genotypes' has incorrect format - 'pedigree' is not a data.frame")
  if (!is.data.frame(genotypes)) stop ("'ped.genotypes' has incorrect format - 'genotypes' is not a data.frame")

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6) stop ("'ped.genotypes' has incorrect format - 'pedigree' must have 6 columns: fid, iid, pid, mid, moi and aff")
  if (nrow(pedigree) < 3) stop ("cannot calculate iR statistic with only 2 isolates") # can't get iR with single pair
  if (any(colnames(pedigree) != c("fid", "iid", "pid", "mid", "moi", "aff")))
    stop ("'ped.genotypes' has incorrect format - 'pedigree' must have columns labelled: fid, iid, pid, mid, moi and aff")
  number.isolates <- nrow(pedigree)
  number.pairs <- (number.isolates * (number.isolates - 1))/2

  # check there are ped.genotypes and pairs to perform analysis
  if (ncol(genotypes) <= 6) stop ("'ped.genotypes' has incorrect format - minimum of 2 isolates required for analysis")
  if (nrow(genotypes) <= 1) stop ("'ped.genotypes' has incorrect format - too few SNPs for analysis")
  if (any(colnames(genotypes)[1:5] != c("chr", "snp_id", "pos_M","pos_bp", "freq")))
    stop("'ped.genotypes' has incorrect format - 'genotypes' must have columns labelled: chr, snp_id, pos_M, pos_bp and freq")
  if (ncol(genotypes) != (number.isolates + 5)) stop ("'ped.genotypes' has incorrect format - some isolates without genotype data")

  # check locus matrix input
  if (ncol(ibd.matrix) < 5) stop ("'ibd.matrix' has incorrect format - must have minimum 5 columns")
  if (any(colnames(ibd.matrix)[1:4] != c("chr", "snp_id", "pos_M", "pos_bp")))
    stop ("'ibd.matrix' has incorrect format - must have first 4 columns labelled: chr, snp_id, pos_M, pos_bp")
  if (ncol(ibd.matrix) != (number.pairs + 4)) stop ("'ibd.matrix' has incorrect format - somes pairs are missing")

  # check groups
  if (!is.null(groups)) {
    if (!is.data.frame(groups)) stop ("'groups' has incorrect format - must be a data.frame")
    if (ncol(groups) != 3) stop ("'groups' has incorrect format - must have 3 columns: fid, iid, group")
    colnames(groups)[1:2] <- c("fid","iid")

    # check isolates belong to a group
    group.names <- paste(groups[,"fid"],groups[,"iid"],sep="/")
    isolate.names <- paste(pedigree[,"fid"],pedigree[,"iid"],sep="/")
    if (!all(isolate.names %in% group.names)) stop ("'groups' has incorrect format - some isolates 'ped.genotypes' are missing from 'groups'")

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
        change.pair <- which(group.pairs[,1] == groups.unique.pairs[i,2] & group.pairs[,2] == groups.unique.pairs[i,1])
        group.pairs[change.pair,1] <- groups.unique.pairs[i,1]
        group.pairs[change.pair,2] <- groups.unique.pairs[i,2]
      }
    }

    # get number of pairwise analyses for each group
    npairs.groups <- NULL
    group.pairs.1 <- paste(group.pairs[,1],group.pairs[,2],sep="/")
    for (i in unique(group.pairs.1)) {
      id_1 <- unlist(strsplit(i,"/"))[1]
      id_2 <- unlist(strsplit(i,"/"))[2]
      npairs <- dim(group.pairs[group.pairs[,1] == id_1 & group.pairs[,2] == id_2,])[1]
      if (is.null(npairs)) npairs <- 1
      npairs.groups <- c(npairs.groups, npairs)
    }
  }

  # calculate iR at each SNP
  locus.pairs <- ibd.matrix[,5:ncol(ibd.matrix)]
  if (!is.null(groups)) {
    locus.prop <- NULL
    na.iR <- NULL
    locus.prop.colnames <- NULL
    for (g in 1:length(unique(group.pairs.1))) {
      if (npairs.groups[g] == 1) {
        locus.g <- rep(NA, nrow(ibd.matrix))
      } else {
        locus.pairs.g <- locus.pairs[,group.pairs.1 == unique(group.pairs.1)[g]]
        locus.g <- iRfunction(locus.pairs.g, genotypes[,"freq"])
      }
      if (all(is.na(locus.g)))
        na.iR <- c(na.iR, g)
      locus.prop <- cbind(locus.prop, locus.g)
      locus.prop.colnames <- c(locus.prop.colnames, unique(group.pairs.1)[g])
    }
    # groups with NA iR
    if (length(na.iR) > 0)
      warning ("NA values generated for: ",paste(unique(group.pairs.1)[na.iR], collapse=", "))

    # name columns
    if (length(unique(group.pairs.1)) != 1)
      colnames(locus.prop) <- unique(group.pairs.1)
  } else {
    locus.prop <- iRfunction(locus.pairs, genotypes[,"freq"])
    if (all(is.na(locus.prop)))
      warning ("NA values generated")
  }

  # melt data.frame if groups
  return.locus.prop <- cbind(ibd.matrix[,1:4], pop=rep(1, nrow(ibd.matrix)), locus.prop)
  if (!is.null(groups) & ncol(return.locus.prop) > 6) {
    return.locus.prop <- data.table::melt(return.locus.prop, id.vars=colnames(return.locus.prop)[1:5])
  } else {
    return.locus.prop <- data.frame(ibd.matrix[,1:4], pop=rep(1, nrow(ibd.matrix)), subpop=1, locus.prop)
  }
  colnames(return.locus.prop) <- c("chr", "snp_id", "pos_M", "pos_bp", "pop", "subpop", "iR")

  return(data.frame(return.locus.prop))
}


