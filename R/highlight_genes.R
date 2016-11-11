#' Plasmodium Falciparum Gene Highlight Dataset
#'
#'
#' Gene annotations for 31 commonly studied Plasmodium falciparum genes,
#' including antimalarial drug resistance genes, vaccine candidates and var genes.
#' Gene annotations are for the reference genome 3D7 were downloaded from
#' \url{http://www.plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/},
#' release PlasmoDB-28, last modified 23/03/2016.
#'
#' @format A data frame with 5 columns of information
#'  \describe{
#' \item{chr}{Chromosomes}
#' \item{start}{Base-pair positions of start of genes}
#' \item{end}{Base-pair positions of end of genes}
#' \item{name}{Gene name, commonly abbreviated}
#' \item{gene_id}{Gene id}
#' }
"highlight_genes"
