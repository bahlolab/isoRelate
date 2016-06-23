#include <Rcpp.h>
using namespace Rcpp;

//' Pair Combinations for Analysis
//' 
//' Creates a data frame containing family IDs and isolate IDs for each pair to be analysed. Each row
//' corresponds to a unique pair.
//' 
//' @param fid A character vector of all family IDs
//' @param iid A character vector of all individual ID
//' @export
// [[Rcpp::export]]
CharacterMatrix isolatePairs(CharacterVector fid, CharacterVector iid) {
  int number_isolates = fid.size();
  CharacterMatrix isolate_pairs(number_isolates*(number_isolates - 1)/2,4);
  int k = 0;
  
  for (int i = 0; i < (number_isolates - 1); i++) {
    for(int j = (i + 1); j < number_isolates; j++){
      isolate_pairs(k,0) = fid[i];
      isolate_pairs(k,1) = iid[i];
      isolate_pairs(k,2) = fid[j];
      isolate_pairs(k,3) = iid[j];
      k += 1;
    }
  }
  return isolate_pairs;
}
