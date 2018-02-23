#include <Rcpp.h>
using namespace Rcpp;

//' Call Genotypes from Haplotype Data
//'
//' \code{haplotypeToGenotype} transforms PLINK haplotype data into genotype data of the form -1, 0, 1 and 2 representing missing genotypes,
//' homozygous reference, heterozygous and homozygous alternative respectively. Haploid isolates are coded as diploid although will not have
//' heterozygous genotypes.
//'
//' @param haplotypes An integer matrix of haplotype data in PLINK format. A allele is denoted 1, B allele is denoted 2 and missing data is denoted 0
//' @param moi An integer vector of multiplicity of infection (MOI) estimates for each isoalte. Isolate MOI estimates should be ordered such that value \code{n} of
//' \code{moi} corresponds to column \code{n} of \code{haplotypes}.
//' @return A matrix with genotype calls where columns correspond to isolates and rows correspond to SNPs
// [[Rcpp::export]]
IntegerMatrix haplotypeToGenotype(IntegerMatrix haplotypes, IntegerVector moi){
  IntegerMatrix genotypes(haplotypes.ncol()/2, haplotypes.nrow());
  int number_isolates = haplotypes.nrow();
  int number_snps = haplotypes.ncol();

  for (int i = 0; i < number_isolates; i++) {
    for (int j = 0; j < number_snps; j=j+2) {
      if(haplotypes(i,j) == 1 && haplotypes(i,j+1) == 1) genotypes(j/2,i) = 0;
      if(haplotypes(i,j) == 2 && haplotypes(i,j+1) == 2) genotypes(j/2,i) = 2;
      if(haplotypes(i,j) == 0 && haplotypes(i,j+1) == 0) genotypes(j/2,i) = -1;
      if((moi[i] == 1) && ((haplotypes(i,j) == 1 && haplotypes(i,j+1) == 2) || (haplotypes(i,j) == 2 && haplotypes(i,j+1) == 1))) genotypes(j/2,i) = -1;
      if((moi[i] == 2) && ((haplotypes(i,j) == 1 && haplotypes(i,j+1) == 2) || (haplotypes(i,j) == 2 && haplotypes(i,j+1) == 1))) genotypes(j/2,i) = 1;
      if(haplotypes(i,j) != 0 && haplotypes(i,j) != 1 && haplotypes(i,j) != 2) genotypes(j/2,i) = -1;
      if(haplotypes(i,j+1) != 0 && haplotypes(i,j+1) != 1 && haplotypes(i,j+1) != 2) genotypes(j/2,i) = -1;
    }
  }
  return genotypes;
}

