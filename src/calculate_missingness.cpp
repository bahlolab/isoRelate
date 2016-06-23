#include <Rcpp.h>
using namespace Rcpp;

//' Calculate Missingness Proportions
//' 
//' Calculates the proportion of missing data for each SNPs or each isolate where missing values are denoted by -1. 
//' Missing values are calculated for each column of \code{genotypes} (where columns are isolates and rows are SNPs),
//' however \code{genotypes} can be transposed to calculate missingness proportions for SNPs.
//' 
//' @param genotypes An integer matrix of genotype data of the form -1, 0, 1 and 2 representing missing genotypes, homozygous reference, heterozygous and homozygous alternative respectively.
//' @return A vector of length \code{n} where \code{n} is the number of columns in \code{genotypes}. 
//' @export
// [[Rcpp::export]]
NumericVector calculateMissingness(IntegerMatrix genotypes) {
  NumericVector proportion_missing(genotypes.ncol());
  int number_snps = genotypes.nrow();
  int number_isolates = genotypes.ncol();
  double number_snps_1 = genotypes.nrow();
  
  for (int i = 0; i < number_isolates; i++) {
    double number_missing = 0.0;
    for (int j = 0; j < number_snps; j++) {
      if(genotypes(j,i) == -1 ) number_missing += 1;
    }
    proportion_missing[i] = number_missing/number_snps_1;
  }
  return proportion_missing;
}
