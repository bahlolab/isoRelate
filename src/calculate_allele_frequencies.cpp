#include <Rcpp.h>
using namespace Rcpp;


//' Calculate Allele Frequencies for SNPs from Genotype Data
//' 
//' \code{calculatePopAlleleFreq} calculates reference allele frequencies for each SNP given genotype data.
//' 
//' @param genotypes An integer matrix of genotype data of the form -1, 0, 1 and 2 representing missing genotypes, homozygous reference, 
//' heterozygous and homozygous alternative respectively. Each column of \code{genotypes} represents a unique isolate and each row of \code{genotypes}
//' represents a unique SNP.
//' @param moi An integer vector of multiplicity of infection (MOI) estimates for each isoalte. Isolate MOI estimates should be ordered such that value \code{n} of 
//' \code{moi} corresponds to column \code{n} of \code{genotypes}.
//' @export
// [[Rcpp::export]]
NumericVector calculatePopAlleleFreq(IntegerMatrix genotypes, IntegerVector moi) {
  NumericVector pop_allele_freqs(genotypes.nrow());
  int number_isolates = genotypes.ncol();
  int number_snps = genotypes.nrow();
  
  for (int t = 0; t < number_snps; t++) {
    double A = 0, B = 0; 
    for (int i = 0; i < number_isolates; i++) {
      if (genotypes(t,i) == 0 && moi[i] == 1)  A += 1;
      if (genotypes(t,i) == 0 && moi[i] == 2)  A += 2;
      if (genotypes(t,i) == 1 && moi[i] == 2){ A += 1; B += 1; }
      if (genotypes(t,i) == 2 && moi[i] == 1)  B += 1;
      if (genotypes(t,i) == 2 && moi[i] == 2)  B += 2;
    }
    if (A + B == 0) pop_allele_freqs[t] = -1; else pop_allele_freqs[t] = A/(A+B);
  }
  return pop_allele_freqs;
}
