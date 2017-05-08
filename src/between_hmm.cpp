#include <Rcpp.h>
using namespace Rcpp;

// Below are all the functions used in the hidden Markov model
// to detect IBD between isoaltes.


//----- round decimal -----

//' Round digits to specified decimal places
//' @param number A number to round
//' @param number The number of digits to round to
// [[Rcpp::export]]
double roundDecimal(double number, int digits){
  double new_number;
  new_number = ceil( number * pow(10,digits) - 0.4999999 )/ pow(10,digits);
  if(new_number == -0) return 0;
  else return new_number;
}


//----- emission probabilities -----

//' The emission probabilities for 2 haploid chromosomes
//' @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
//' @param genotype_1 The genotype for isolate 1 from the pair for SNP i
//' @param genotype_2 The genotype for isolate 2 from the pair for SNP i
//' @param ibd The IBD state
// [[Rcpp::export]]
double emissionProbHH(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  double alt_allele_freq = 1 - pop_allele_freq;
  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,2), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal( pop_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal( alt_allele_freq, 6);

  return 0;
}

//' The emission probabilities for 1 haploid chromosome and 1 diploid chromosome
//' @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
//' @param genotype_1 The genotype for isolate 1 from the pair for SNP i
//' @param genotype_2 The genotype for isolate 2 from the pair for SNP i
//' @param ibd The IBD state
//' @param male_column The haploid isolate from the pair. Either 1 or 2
//' @param female_column The diploid isolate from the pair. Either 1 or 2
// [[Rcpp::export]]
double emissionProbHD(double pop_allele_freq, int genotype_1, int genotype_2, int ibd, int male_column, int female_column) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  int geno_male, geno_female;

  if(male_column == 1) { geno_male = genotype_1; geno_female = genotype_2; }
  if(male_column == 2) { geno_male = genotype_2; geno_female = genotype_1; }

  if(geno_female == 0 && geno_male == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,3), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 0 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 0) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(geno_female == 2 && geno_male == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,3), 6);

  if(geno_female == 0 && geno_male == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(geno_female == 0 && geno_male == 2 && ibd == 1) return 0;
  if(geno_female == 1 && geno_male == 0 && ibd == 1) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 1 && geno_male == 2 && ibd == 1) return roundDecimal( pop_allele_freq*alt_allele_freq, 6);
  if(geno_female == 2 && geno_male == 0 && ibd == 1) return 0;
  if(geno_female == 2 && geno_male == 2 && ibd == 1) return roundDecimal( pow(alt_allele_freq,2), 6);

  return 0;
}

//' The emission probabilities for 2 diploid chromosomes
//' @param pop_allele_freq The population allele frequency for SNP i. This corresponds to the reference allele
//' @param genotype_1 The genotype for isolate 1 from the pair for SNP i
//' @param genotype_2 The genotype for isolate 2 from the pair for SNP i
//' @param ibd The IBD state
// [[Rcpp::export]]
double emissionProbDD(double pop_allele_freq, int genotype_1, int genotype_2, int ibd) {
  if(pop_allele_freq > 1) { pop_allele_freq = 1; }
  double alt_allele_freq = 1 - pop_allele_freq;

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,4), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 0) return roundDecimal( 2*pow(pop_allele_freq,3)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 0) return roundDecimal( 4*pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 0) return roundDecimal( pow(pop_allele_freq,2)*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 0) return roundDecimal( 2*pop_allele_freq*pow(alt_allele_freq,3), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 0) return roundDecimal( pow(alt_allele_freq,4), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,3), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 1) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 1) return roundDecimal( pow(pop_allele_freq,2)*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 1) return roundDecimal( (pow(pop_allele_freq,2)*alt_allele_freq + pop_allele_freq*pow(alt_allele_freq,2)), 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 1) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 1) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 1) return roundDecimal( pop_allele_freq*pow(alt_allele_freq,2), 6);
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 1) return roundDecimal( pow(alt_allele_freq,3), 6);

  if(genotype_1 == 0 && genotype_2 == 0 && ibd == 2) return roundDecimal( pow(pop_allele_freq,2), 6);
  if(genotype_1 == 0 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 0 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 1 && genotype_2 == 1 && ibd == 2) return roundDecimal( 2*pop_allele_freq*alt_allele_freq, 6);
  if(genotype_1 == 1 && genotype_2 == 2 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 0 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 1 && ibd == 2) return 0;
  if(genotype_1 == 2 && genotype_2 == 2 && ibd == 2) return roundDecimal( pow(alt_allele_freq,2), 6);

  return 0;
}


//----- transition probabilities -----

//' The transition probabilities for 2 haploid chromosomes
//' @param omega_0 The probability of sharing 0 alleles IBD
//' @param meiosis The number of meiosis separating the two isoaltes
//' @param dist_cM The genetic map distance (cM) between SNP i and SNP j
//' @param ibd_current The IBD state of SNP j
//' @param ibd_previous The IBD state of SNP i
// [[Rcpp::export]]
double transitionProbHH(double omega_0, int meiosis, double dist_cM, int ibd_current, int ibd_previous) {
  double dist_M = dist_cM;
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}

//' The transition probabilities for 1 haploid and 1 diploid chromosome
//' @param omega_0 The probability of sharing 0 alleles IBD
//' @param meiosis The number of meiosis separating the two isoaltes
//' @param dist_cM The genetic map distance (cM) between SNP i and SNP j
//' @param ibd_current The IBD state of SNP j
//' @param ibd_previous The IBD state of SNP i
// [[Rcpp::export]]
double transitionProbHD(double omega_0, int meiosis, double dist_cM, int ibd_current, int ibd_previous) {
  double dist_M = dist_cM;
  double omega_1 = 1 - omega_0;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (omega_0 + omega_1 * exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( (omega_0 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( (omega_1 * (1.0 - exp(-alpha * dist_M))), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( (omega_1 + omega_0 * exp(-alpha * dist_M)), 6);

  return 0;
}

//' The transition probabilities for 2 diploid chromosomes
//' @param omega_0 The probability of sharing 0 alleles IBD
//' @param meiosis The number of meiosis separating the two isoaltes
//' @param dist_cM The genetic map distance (cM) between SNP i and SNP j
//' @param ibd_current The IBD state of SNP j
//' @param ibd_previous The IBD state of SNP i
// [[Rcpp::export]]
double transitionProbDD(double omega_0, double omega_1, double omega_2, int meiosis, double dist_cM, int ibd_current, int ibd_previous) {
  double dist_M = dist_cM;
  double theta = 0.5 * (1.0 - exp(-2.0 * dist_M));
  double alpha = -meiosis * log(1.0 - theta);

  double T02 = (exp(-alpha * omega_1 * dist_M) * omega_2)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_0 * omega_1)/(omega_1 - 1) + omega_2;
  double T20 = (exp(-alpha * omega_1 * dist_M) * omega_0)/(omega_1 - 1) + exp(-alpha * dist_M) * omega_1 + (exp(-alpha * dist_M) * omega_2 * omega_1)/(omega_1 - 1) + omega_0;

  if(ibd_previous == 0 && ibd_current == 0) return roundDecimal( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T02), 6);
  if(ibd_previous == 0 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 0 && ibd_current == 2) return roundDecimal( T02, 6);
  if(ibd_previous == 1 && ibd_current == 0) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_0), 6);
  if(ibd_previous == 1 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1 + exp(-alpha * dist_M)), 6);
  if(ibd_previous == 1 && ibd_current == 2) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_2), 6);
  if(ibd_previous == 2 && ibd_current == 0) return roundDecimal( T20, 6);
  if(ibd_previous == 2 && ibd_current == 1) return roundDecimal( ((1 - exp(-alpha * dist_M)) * omega_1), 6);
  if(ibd_previous == 2 && ibd_current == 2) return roundDecimal( (1 - (1 - exp(-alpha * dist_M)) * omega_1 - T20), 6);

  return 0;
}


//----- genotype error probabilities -----

//' The genotyping error probability for 1 haploid chromosome
//' @param truth The true genotype
//' @param observed The observed genotype
//' @param error The genotype error rate
// [[Rcpp::export]]
double genotypeErrorH(int truth, int observed, double error){
  if(truth == observed) return 1 - error;
  else return error;
}

//' The genotyping error probability for 1 diploid chromosome
//' @param truth The true genotype
//' @param observed The observed genotype
//' @param error The genotype error rate
// [[Rcpp::export]]
double genotypeErrorD(int truth, int observed, double error){
  if(truth == 0 && observed ==0) return pow(1 - error,2);
  if(truth == 0 && observed ==1) return 2 * (1 - error) * error;
  if(truth == 0 && observed ==2) return pow(error,2);
  if(truth == 1 && observed ==0) return (1 - error) * error;
  if(truth == 1 && observed ==1) return pow(1 - error,2) + pow(error,2);
  if(truth == 1 && observed ==2) return (1 - error) * error;
  if(truth == 2 && observed ==0) return pow(error,2);
  if(truth == 2 && observed ==1) return 2 * (1 - error) * error;
  if(truth == 2 && observed ==2) return pow(1 - error,2);

  return 0;
}


//----- true genotypes -----

//' Matrices of all possible genotype combinations between pairs, given MOI
//' @param gender_1 The MOI estimate of isolate 1
//' @param gender_2 The MOI estimate of isolate 2
// [[Rcpp::export]]
IntegerMatrix trueGenotypes(int gender_1, int gender_2){
  // 2 haploid chromosomes
  if(gender_1 == 1 && gender_2 == 1){
    IntegerMatrix true_genotype(4,2) ;
    true_genotype(0,0)=0; true_genotype(0,1)=0;
    true_genotype(1,0)=0; true_genotype(1,1)=2;
    true_genotype(2,0)=2; true_genotype(2,1)=0;
    true_genotype(3,0)=2; true_genotype(3,1)=2;
    return true_genotype ;
  }
  // 1 haploid chromosome and 1 diploid
  if((gender_1 == 1 && gender_2 == 2) || (gender_1 == 2 && gender_2 == 1)){
    IntegerMatrix true_genotype(6,2) ;
    true_genotype(0,0)=0; true_genotype(0,1)=0;
    true_genotype(1,0)=0; true_genotype(1,1)=2;
    true_genotype(2,0)=1; true_genotype(2,1)=0;
    true_genotype(3,0)=1; true_genotype(3,1)=2;
    true_genotype(4,0)=2; true_genotype(4,1)=0;
    true_genotype(5,0)=2; true_genotype(5,1)=2;
    return true_genotype ;
  }
  // 2 diploid chromosomes
  if(gender_1 == 2 && gender_2 == 2){
    IntegerMatrix true_genotype(9,2) ;
    true_genotype(0,0)=0; true_genotype(0,1)=0;
    true_genotype(1,0)=0; true_genotype(1,1)=1;
    true_genotype(2,0)=0; true_genotype(2,1)=2;
    true_genotype(3,0)=1; true_genotype(3,1)=0;
    true_genotype(4,0)=1; true_genotype(4,1)=1;
    true_genotype(5,0)=1; true_genotype(5,1)=2;
    true_genotype(6,0)=2; true_genotype(6,1)=0;
    true_genotype(7,0)=2; true_genotype(7,1)=1;
    true_genotype(8,0)=2; true_genotype(8,1)=2;
    return true_genotype ;
  }
  return 0;
}


//----- emission probability error summation -----

//' Calculating the emission probability sumation when missing genotype calls present
//' @param pop_allele_freq The population allele frequency of SNP t
//' @param genotype_1 The genotype of isolate 1 at SNP t
//' @param genotype_2 The genotype of isolate 2 at SNP t
//' @param error The genotype error rate
//' @param gender_1 The MOI estimate of isolate 1
//' @param gender_2 The MOI estimate of isolate 2
//' @param ibd_j The IBD state
// [[Rcpp::export]]
double emissionProbMissingGeno(double pop_allele_freq, int genotype_1, int genotype_2, double error, int gender_1, int gender_2, int ibd_j){
  IntegerMatrix trueGenotype ;
  trueGenotype = trueGenotypes(gender_1,gender_2) ; // All combinations of genotypes
  double emission_prob_error ;
  int noGenotypes, miss, miss1, miss2 ;
  noGenotypes = trueGenotype.nrow() ; // total number of combinations of genotpes
  emission_prob_error = 0 ;

  for(int g = 0; g<noGenotypes; ++g){

    if(genotype_1 == -1 && genotype_2 != -1){ // first genotype missing
      if(gender_1 == 2 && gender_2 == 2){
        for(int miss = 0; miss < 3; ++miss){
          emission_prob_error += emissionProbDD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),miss,error) * genotypeErrorD(trueGenotype(g,1),genotype_2,error) ;
        }
      }
      if(gender_1 == 1 && gender_2 == 2){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),miss,error) * genotypeErrorD(trueGenotype(g,0),genotype_2,error) ;
        }
      }
      if(gender_1 == 2 && gender_2 == 1){
        for(int miss = 0; miss < 3; ++miss){
          emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),miss,error) * genotypeErrorH(trueGenotype(g,1),genotype_2,error) ;
        }
      }
      if(gender_1 == 1 && gender_2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emission_prob_error += emissionProbHH(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),miss,error) * genotypeErrorH(trueGenotype(g,1),genotype_2,error) ;
        }
      }
    }

    if(genotype_1 != -1 && genotype_2 == -1){ // second genotype missing
      if(gender_1 == 2 && gender_2 == 2){
        for(int miss = 0; miss < 3; ++miss){
          emission_prob_error += emissionProbDD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),genotype_1,error) * genotypeErrorD(trueGenotype(g,1),miss,error) ;
        }
      }
      if(gender_1 == 1 && gender_2 == 2){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),genotype_1,error) * genotypeErrorD(trueGenotype(g,0),miss,error) ;
        }
      }
      if(gender_1 == 2 && gender_2 == 1){
        for(int miss = 0; miss < 3; ++miss){
          emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),genotype_1,error) * genotypeErrorH(trueGenotype(g,1),miss,error) ;
        }
      }
      if(gender_1 == 1 && gender_2 == 1){
        for(int miss = 0; miss < 3; miss = miss + 2){
          emission_prob_error += emissionProbHH(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),genotype_1,error) * genotypeErrorH(trueGenotype(g,1),miss,error) ;
        }
      }
    }

    if(genotype_1 == -1 && genotype_2 == -1){ // both genotypes missing
      if(gender_1 == 2 && gender_2 == 2){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emission_prob_error += emissionProbDD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),miss1,error) * genotypeErrorD(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(gender_1 == 1 && gender_2 == 2){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; ++miss2){
            emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),miss1,error) * genotypeErrorD(trueGenotype(g,0),miss2,error) ;
          }
        }
      }
      if(gender_1 == 2 && gender_2 == 1){
        for(int miss1 = 0; miss1 < 3; ++miss1){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),miss1,error) * genotypeErrorH(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
      if(gender_1 == 1 && gender_2 == 1){
        for(int miss1 = 0; miss1 < 3; miss1 = miss1 + 2){
          for(int miss2 = 0; miss2 < 3; miss2 = miss2 + 2){
            emission_prob_error += emissionProbHH(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),miss1,error) * genotypeErrorH(trueGenotype(g,1),miss2,error) ;
          }
        }
      }
    }

    if(genotype_1 != -1 && genotype_2 != -1){ // no genotypes missing
      if(gender_1 == 2 && gender_2 == 2){
        emission_prob_error += emissionProbDD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorD(trueGenotype(g,0),genotype_1,error) * genotypeErrorD(trueGenotype(g,1),genotype_2,error) ;
      }
      if(gender_1 == 1 && gender_2 == 2){
        emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,1),trueGenotype(g,0),ibd_j,1,2) * genotypeErrorH(trueGenotype(g,1),genotype_1,error) * genotypeErrorD(trueGenotype(g,0),genotype_2,error) ;
      }
      if(gender_1 == 2 && gender_2 == 1){
        emission_prob_error += emissionProbHD(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j,2,1) * genotypeErrorD(trueGenotype(g,0),genotype_1,error) * genotypeErrorH(trueGenotype(g,1),genotype_2,error) ;
      }
      if(gender_1 == 1 && gender_2 == 1){
        emission_prob_error += emissionProbHH(pop_allele_freq,trueGenotype(g,0),trueGenotype(g,1),ibd_j) * genotypeErrorH(trueGenotype(g,0),genotype_1,error) * genotypeErrorH(trueGenotype(g,1),genotype_2,error) ;
      }
    }
  }
  return emission_prob_error ;
}


//----- alpha -----

//' Calculate alpha
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
NumericMatrix calculateAlpha(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, double error, int gender_1, int gender_2){
  // Initialising all parameters:
  NumericVector scale(number_snps);
  NumericVector alpha_a(number_snps);
  NumericMatrix alpha(number_snps, number_states);
  NumericMatrix alpha_hat(number_snps, number_states);
  double alpha_a_sum, alpha_sum, emission_prob;

  // Initialisation
  alpha_sum = 0;
  for(int i = 0; i<number_states; ++i){
    emission_prob = emissionProbMissingGeno(pop_allele_freqs[0], genotypes(0,0), genotypes(0,1), error, gender_1, gender_2, i) ;
    alpha(0,i) = initial_prob[i] * emission_prob;
    alpha_sum += alpha(0,i);
  }
  scale[0] = roundDecimal( 1/alpha_sum, 4);
  for(int i = 0; i<number_states; ++i){
    alpha_hat(0,i) = roundDecimal( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(number_snps-1); ++t){
    alpha_sum = 0;
    for(int j = 0; j<number_states; ++j){
      alpha_a_sum = 0;
      for(int i = 0; i<number_states; ++i){
        if(gender_1 == 2 && gender_2 == 2){
          alpha_a[i] = alpha_hat(t,i) * transitionProbDD(initial_prob[0],initial_prob[1],initial_prob[2],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        if((gender_1 == 1 && gender_2 == 2) || (gender_1 == 2 && gender_2 == 1) ){
          alpha_a[i] = alpha_hat(t,i) * transitionProbHD(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        if(gender_1 == 1 && gender_2 == 1){
          alpha_a[i] = alpha_hat(t,i) * transitionProbHH(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        alpha_a_sum += alpha_a[i];
      }
      emission_prob = emissionProbMissingGeno(pop_allele_freqs[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender_1, gender_2, j) ;
      alpha(t+1,j) = alpha_a_sum * emission_prob;
      alpha_sum += alpha(t+1,j);
     }
     scale[t+1] = roundDecimal( 1/alpha_sum, 4);
     for(int i = 0; i<number_states; ++i){
       alpha_hat(t+1,i) = roundDecimal( alpha(t+1,i) * scale[t+1], 4);
     }
   }
   return alpha_hat;
}


//----- scale -----

//' Calculate scale
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
NumericVector calculateScale(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, double error, int gender_1, int gender_2){
  // Initialising all parameters:
  NumericVector scale(number_snps);
  NumericVector alpha_a(number_snps);
  NumericMatrix alpha(number_snps, number_states);
  NumericMatrix alpha_hat(number_snps, number_states);
  double alpha_a_sum, alpha_sum, emission_prob;

  // Initialisation
  alpha_sum = 0;
  for(int i = 0; i<number_states; ++i){
    emission_prob = emissionProbMissingGeno(pop_allele_freqs[0], genotypes(0,0), genotypes(0,1), error, gender_1, gender_2, i) ;
    alpha(0,i) = initial_prob[i] * emission_prob;
    alpha_sum += alpha(0,i);
  }
  scale[0] = roundDecimal( 1/alpha_sum, 4);
  for(int i = 0; i<number_states; ++i){
    alpha_hat(0,i) = roundDecimal( alpha(0,i) * scale[0], 4);
  }

  // Induction
  for(int t = 0; t<(number_snps-1); ++t){
    alpha_sum = 0;
    for(int j = 0; j<number_states; ++j){
      alpha_a_sum = 0;
      for(int i = 0; i<number_states; ++i){
        if(gender_1 == 2 && gender_2 == 2){
          alpha_a[i] = alpha_hat(t,i) * transitionProbDD(initial_prob[0],initial_prob[1],initial_prob[2],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        if((gender_1 == 1 && gender_2 == 2) || (gender_1 == 2 && gender_2 == 1) ){
          alpha_a[i] = alpha_hat(t,i) * transitionProbHD(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        if(gender_1 == 1 && gender_2 == 1){
          alpha_a[i] = alpha_hat(t,i) * transitionProbHH(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) ;
        }
        alpha_a_sum += alpha_a[i];
      }
      emission_prob = emissionProbMissingGeno(pop_allele_freqs[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender_1, gender_2, j) ;
      alpha(t+1,j) = alpha_a_sum * emission_prob;
      alpha_sum += alpha(t+1,j);
     }
     scale[t+1] = roundDecimal( 1/alpha_sum, 4);
     for(int i = 0; i<number_states; ++i){
       alpha_hat(t+1,i) = roundDecimal( alpha(t+1,i) * scale[t+1], 4);
     }
   }
   return scale;
}


//----- beta -----

//' Calculate beta
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param scale A numeric vector containing the scaling values used to scale alpha
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
NumericMatrix calculateBeta(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, NumericVector scale, double error, int gender_1, int gender_2){
  // Initialising all parameters:
  NumericMatrix beta(number_snps, number_states);
  NumericMatrix beta_hat(number_snps, number_states);
  double beta_sum, emission_prob;
  const int T = number_snps-1;

  // Initialisation
  for(int i = 0; i<number_states; ++i){
    beta(T,i)    = 1;
    beta_hat(T,i) = roundDecimal( beta(T,i) * scale[T], 4);
  }

  // Induction
  for(int t = (T-1); t>=0; t--){
    for(int i = 0; i<number_states; ++i){
      beta_sum = 0;
      for(int j = 0; j<number_states; ++j){
        emission_prob = emissionProbMissingGeno(pop_allele_freqs[t+1], genotypes(t+1,0), genotypes(t+1,1), error, gender_1, gender_2, j) ;
        if(gender_1 == 2 && gender_2 == 2){
          beta(t,j) = transitionProbDD(initial_prob[0],initial_prob[1],initial_prob[2],meiosis,positions_cM[t+1]-positions_cM[t],j,i) * emission_prob * beta_hat(t+1,j) ;
        }
        if((gender_1 == 1 && gender_2 == 2) || (gender_1 == 2 && gender_2 == 1)){
          beta(t,j) = transitionProbHD(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) * emission_prob * beta_hat(t+1,j) ;
        }
        if(gender_1 == 1 && gender_2 == 1){
          beta(t,j) = transitionProbHH(initial_prob[0],meiosis,positions_cM[t+1]-positions_cM[t],j,i) * emission_prob * beta_hat(t+1,j) ;
        }
        beta_sum += beta(t,j);
      }
      beta_hat(t,i) = roundDecimal( beta_sum * scale[t], 4);
    }
  }
  return beta_hat;
}


//----- Viterbi -----

//' Calculate the Viterbi sequence
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
IntegerVector calculateViterbi(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, double error, int gender_1, int gender_2){
  // Initialising all parameters:
  IntegerVector q_star(number_snps);
  NumericVector delta_a(number_snps);
  NumericMatrix delta(number_snps, number_states);
  NumericMatrix psi(number_snps, number_states);
  double emission_prob;
  const int T = number_snps-1;

  // Initialisation
  for(int i = 0; i<number_states; ++i){
    psi(0,i) = 0;
    emission_prob = emissionProbMissingGeno(pop_allele_freqs[0], genotypes(0,0), genotypes(0,1), error, gender_1, gender_2, i) ;
    delta(0,i) = log(initial_prob[i]) + log(emission_prob) ;
  }

  // Recursion
  for(int t = 1; t<number_snps; ++t){
    for(int j = 0; j<number_states; ++j){
      for(int i = 0; i<number_states; ++i){
        if(gender_1 == 2 && gender_2 == 2){
          log_trans = log(transitionProbDD(initial_prob[0],initial_prob[1],initial_prob[2],meiosis,positions_cM[t]-positions_cM[t-1],j,i)) ;
          if (log_trans < -100000) log_trans = log(0.000001) ;
          delta_a[i] = delta(t-1,i) +  log_trans ;
        }
        if((gender_1 == 1 && gender_2 == 2) || (gender_1 == 2 && gender_2 == 1)){
          log_trans = log(transitionProbHD(initial_prob[0],meiosis,positions_cM[t]-positions_cM[t-1],j,i)) ;
          if (log_trans < -100000) log_trans = log(0.000001) ;
          delta_a[i] = delta(t-1,i) +  log_trans ;
        }
        if(gender_1 == 1 && gender_2 == 1){
          log_trans = log(transitionProbHH(initial_prob[0],meiosis,positions_cM[t]-positions_cM[t-1],j,i)) ;
          if (log_trans < -100000) log_trans = log(0.000001) ;
          delta_a[i] = delta(t-1,i) +  log_trans ;
        }
      }
      emission_prob = emissionProbMissingGeno(pop_allele_freqs[t], genotypes(t,0), genotypes(t,1), error, gender_1, gender_2, j) ;
      double max_num = delta_a[0];
      int max_arg = 0;
      for(int k = 0; k<number_states; ++k){
        if(delta_a[k]>max_num){
          max_num = delta_a[k];
          max_arg = k; // could be problematic if there are more than 1 maximum
        }
      }
      delta(t,j) = max_num + log(emission_prob) ;
      psi(t,j) = max_arg ;
    }
  }

  // Termination
  q_star[T] = 0;
  double prob_star = delta(T,0);
  for(int k = 0; k<number_states; ++k){
    if(delta(T,k)>prob_star){
      q_star[T] = k; // could be problematic if there are more than 1 maximum
      prob_star = delta(T,k);
    }
  }

  // Path backtracking
  for(int t = (T-1); t>=0; t--){
    q_star[t] = psi(t+1,q_star[t+1]); // Not sure if this will work without q_star being a constant
  }
  return q_star;
}


//----- gamma -----

//' Calculate gamma
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
NumericMatrix calculateGamma(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, double error, int gender_1, int gender_2){
  NumericMatrix gamma(number_snps, number_states);
  NumericMatrix alpha(number_snps, number_states);
  NumericMatrix beta(number_snps, number_states);
  NumericVector scale(number_snps);
  NumericVector alpha_beta(number_snps);
  double alpha_beta_sum;

  alpha = calculateAlpha(number_states, initial_prob, meiosis, number_snps, genotypes, pop_allele_freqs, positions_cM, error, gender_1, gender_2);
  scale = calculateScale(number_states, initial_prob, meiosis, number_snps, genotypes, pop_allele_freqs, positions_cM, error, gender_1, gender_2);
  beta  = calculateBeta(number_states, initial_prob, meiosis, number_snps, genotypes, pop_allele_freqs, positions_cM, scale, error, gender_1, gender_2);

  for(int t = 0; t<number_snps; ++t){
    alpha_beta_sum = 0;
    for(int i = 0; i<number_states; ++i){
      alpha_beta[i] = alpha(t,i) * beta(t,i);
      alpha_beta_sum += alpha_beta[i];
    }
    for(int i = 0; i<number_states; ++i){
      gamma(t,i) = roundDecimal( (alpha(t,i) * beta(t,i))/alpha_beta_sum, 3) ;
    }
  }
  return gamma;
}


//----- log-likelihood -----

//' Calculate the log-likelihood of the data
//' @param number_states Integer. The number of IBD states in the model
//' @param initial_prob A numeric vector containing the initial state probabilities
//' @param meiosis Integer. The number of meiosis separating the two isolates
//' @param number_snps Integer. The number of SNPs
//' @param genotypes A integer martix containing the genotype calls for a pair of isolates
//' @param pop_allele_freqs A numeric vector of population allele frequencies
//' @param positions_cM A numeric vector of SNP genetic map positions in cM
//' @param error Numeric. The genotype error rate
//' @param gender_1 Integer. The MOI estimate of isolate 1
//' @param gender_2 Integer. The MOI estimate of isolate 2
// [[Rcpp::export]]
double calculateLogLikelihood(const int number_states, NumericVector initial_prob, int meiosis, const int number_snps, IntegerMatrix genotypes, NumericVector pop_allele_freqs, NumericVector positions_cM, double error, int gender_1, int gender_2){
  double sum_scale = 0;
  NumericVector scale(number_snps);

  scale = calculateScale(number_states, initial_prob, meiosis, number_snps, genotypes, pop_allele_freqs, positions_cM, error, gender_1, gender_2);

  for(int t = 0; t<number_snps; t++){
    sum_scale += log(scale[t]);
  }
  return roundDecimal( -sum_scale, 4);
}
