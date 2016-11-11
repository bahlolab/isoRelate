#' Plasmodium Falciparum Gene Annotation Dataset
#'
#'
#' Gene annotations for the reference genome 3D7 were downloaded from
#' \url{http://www.plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/},
#' release PlasmoDB-28, last modified 23/03/2016.
#'
#' @format A data frame with 6 columns of information
#'  \describe{
#' \item{chr}{Chromosomes}
#' \item{start}{Base-pair positions of start of genes}
#' \item{end}{Base-pair positions of end of genes}
#' \item{strand}{The positive or negaitve gene strand}
#' \item{name}{Gene name, commonly abbreviated}
#' \item{gene_id}{Gene id}
#' }
"annotation_genes"
