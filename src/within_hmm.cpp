#include <Rcpp.h>
using namespace Rcpp;

// Below are all the functions used in the hidden Markov model
// to detect IBD within isoaltes.


//----- round decimal -----

//' Round digits to specified decimal places
//' @param number A number to round
//' @param number The number of digits to round to
// [[Rcpp::export]]
double roundDecimal_w(double number, int digits){
  double new_number;
  new_number = ceil( number * pow(10,digits) - 0.4999999 )/ pow(10,digits);
  if(new_number == -0) return 0;
  else return new_number;
}


//----- emission probabilities -----

//' The emission probability for isolate i, SNP j in state k.
//' @param BAFvalue A vector of length 1 with the BAF value for isolate i, SNP j
//' @param meanBAF0 A vector of possible major clone proportions
//' @param sdBAF0 A vector of standard deviations for; one for each major clone proportion
//' @param state0 The IBD state 
//' @param pfb The population BAF for SNP j
// [[Rcpp::export]]
NumericVector BAFemprobs(NumericVector BAFvalue, NumericVector meanBAF0, NumericVector sdBAF0, int state0, double pfb){ 
  NumericVector baf_emprob;
  double mean0 = meanBAF0[0];
  double mean50 = meanBAF0[state0-1];
  double sd0 = sdBAF0[0];
  double sd50 = sdBAF0[state0-1]; 
  
  if(state0 == 1) { // normal (A,B)
                    if(BAFvalue[0] == 0.0){ baf_emprob = (1-pfb); }
                    if(BAFvalue[0] == 1.0){ baf_emprob = pfb; }
                    if(BAFvalue[0] != 0.0 & BAFvalue[0] != 1.0){ 
                      baf_emprob = (1-pfb)*dnorm(BAFvalue,mean0,sd0) + pfb*dnorm(BAFvalue,1-mean0,sd0);
                    }
  }
  if(state0 != 1) { // multiple clones
                    if(BAFvalue[0] == 0.0){ baf_emprob = (1-pfb)*(1-pfb); }
                    if(BAFvalue[0] == 1.0){ baf_emprob = pfb*pfb; }
                    if(BAFvalue[0] != 0.0 & BAFvalue[0] != 1.0){ 
                      baf_emprob = (1-pfb)*(1-pfb)*dnorm(BAFvalue,mean0,sd0) + pfb*(1-pfb)*dnorm(BAFvalue,mean50,sd50) +
                        pfb*(1-pfb)*dnorm(BAFvalue,1-mean50,sd50) + pfb*pfb*dnorm(BAFvalue,1-mean0,sd0);
                    }
  }
  return baf_emprob;
}


//' The emission probabilities for isolate i calculated over all SNPs and states
//' @param BAFmatrix A vector of BAF values for isolate i
//' @param alleleFreq A vector of population BAF values
//' @param BAFmean A vector of possible major clone proportions
//' @param BAFsd A vector of standard deviations for; one for each major clone proportion
//' @param positionBP A vector of genetic map positions for all SNPs in base-pairs
//' @param bin A matrix of bins with parameter estimates (proportion estimates with corresponding SD)
// [[Rcpp::export]]
NumericMatrix calculate_em_probs(NumericVector BAFmatrix, NumericVector alleleFreq, NumericMatrix BAFmean, NumericMatrix BAFsd, IntegerVector positionBP, IntegerMatrix bin){
  int noSNPs = BAFmatrix.size();
  int noStates = BAFmean.ncol();
  NumericMatrix emProbs(noSNPs,noStates); 
  NumericVector BAF(noStates);
  double BAFsum;
  double noBins = bin.nrow();
  NumericVector BAFmean0(noStates);
  NumericVector BAFsd0(noStates);
  
  for(int t = 0; t < noSNPs; ++t){
    BAFsum = 0.0;
    for(int j = 0; j < noBins; ++j){
      if(positionBP[t] >= bin(j,0) && positionBP[t] <= bin(j,1)){
        BAFmean0 = BAFmean(j,_);
        BAFsd0 = BAFsd(j,_);
      }
    }
    for(int i = 1; i <= noStates; ++i){
      NumericVector BAFvalue(1,BAFmatrix[t]);
      BAF[i-1] = BAFemprobs(BAFvalue, BAFmean0, BAFsd0, i, alleleFreq[t])[0];
      BAFsum += BAF[i-1];
    }
    emProbs(t,_) = BAF/BAFsum;
  }
  return emProbs;
}



//----- transition probabilities -----

//' The transition probabilities for 2 haploid chromosomes
//' @param dist The distance between adjacent SNPs
//' @param transParam A matrix of transition rate probabilities
// [[Rcpp::export]]
NumericMatrix transMatrix(double dist, NumericMatrix transParam){
  double D = 100000;
  double offdiagonalSum;
  int noStates = transParam.ncol();
  NumericMatrix A(noStates,noStates);
  
  for(int i = 1; i <= noStates; ++i){
    offdiagonalSum = 0;
    for(int j = 1; j <= noStates; ++j){
      if(i != j){
        A(i-1,j-1) = transParam(i-1,j-1) * (1 - exp(-dist/D));
        if(A(i-1,j-1) > 1){ A(i-1,j-1) = 0.999; }
        offdiagonalSum += A(i-1,j-1);
      }
    }
    if(offdiagonalSum >= 1){
      for(int j = 1; j <= noStates; ++j){
        A(i-1,j-1) /= (offdiagonalSum/0.999);
      }
      offdiagonalSum = 0.999;
    }
    A(i-1,i-1) = 1 - offdiagonalSum;
  }
  return A;
}


//----- alpha -----

//' Calculate alpha 
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param initalProbs A numeric vector of initial state probabilities. One probability per state
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
// [[Rcpp::export]]
NumericMatrix calculate_alpha_w(NumericMatrix emissionProbs, NumericVector initalProbs, NumericVector positionBP, NumericMatrix transParam){
  int noSNPs = emissionProbs.nrow();
  int noStates = initalProbs.size();
  NumericMatrix transProbs(noStates,noStates);
  NumericMatrix alpha(noSNPs,noStates);
  NumericMatrix alphaHat(noSNPs,noStates);
  NumericVector scale(noSNPs);
  double alphaSum;
  
  // Initialization
  for(int i = 1; i <= noStates; ++i){
    alpha(0,i-1) = initalProbs[i-1] * emissionProbs(0,i-1);
    alphaSum += alpha(0,i-1);
  }
  scale[0] = roundDecimal_w( 1/alphaSum, 8);
  for(int i = 1; i <= noStates; ++i){
    alphaHat(0,i-1) = roundDecimal_w( alpha(0,i-1) * scale[0], 8);
  }
  
  // Induction
  for(int t = 1; t < noSNPs; ++t){
    transProbs = transMatrix(positionBP[t] - positionBP[t-1], transParam);
    alphaSum = 0.0;
    for(int j = 1; j <= noStates; ++j){
      double alphaAsum = 0.0;
      for(int i = 1; i <= noStates; ++i){
        alphaAsum += alphaHat(t-1,i-1) * transProbs(j-1, i-1);
      }
      alpha(t,j-1) = alphaAsum * emissionProbs(t,j-1);
      alphaSum += alpha(t,j-1);
    }
    scale[t] = roundDecimal_w( 1/alphaSum, 8);
    for(int i = 1; i <= noStates; ++i){
      alphaHat(t,i-1) = roundDecimal_w( alpha(t,i-1) * scale[t], 8);
    }
  }
  return alphaHat; 
}


//----- scale -----

//' Calculate scale
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param initalProbs A numeric vector of initial state probabilities. One probability per state
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
// [[Rcpp::export]]
NumericVector calculate_scale_w(NumericMatrix emissionProbs, NumericVector initalProbs, NumericVector positionBP, NumericMatrix transParam){
  int noSNPs = emissionProbs.nrow();
  int noStates = initalProbs.size();
  NumericMatrix transProbs(noStates,noStates);
  NumericMatrix alpha(noSNPs,noStates);
  NumericMatrix alphaHat(noSNPs,noStates);
  NumericVector scale(noSNPs);
  double alphaSum;
  
  // Initialization
  for(int i = 1; i <= noStates; ++i){
    alpha(0,i-1) = initalProbs[i-1] * emissionProbs(0,i-1);
    alphaSum += alpha(0,i-1);
  }
  scale[0] = roundDecimal_w( 1/alphaSum, 8);
  for(int i = 1; i <= noStates; ++i){
    alphaHat(0,i-1) = roundDecimal_w( alpha(0,i-1) * scale[0], 8);
  }
  
  // Induction
  for(int t = 1; t < noSNPs; ++t){
    transProbs = transMatrix(positionBP[t] - positionBP[t-1], transParam);
    alphaSum = 0.0;
    for(int j = 1; j <= noStates; ++j){
      double alphaAsum = 0.0;
      for(int i = 1; i <= noStates; ++i){
        alphaAsum += alphaHat(t-1,i-1) * transProbs(j-1, i-1);
      }
      alpha(t,j-1) = alphaAsum * emissionProbs(t,j-1);
      alphaSum += alpha(t,j-1);
    }
    scale[t] = roundDecimal_w( 1/alphaSum, 8);
    for(int i = 1; i <= noStates; ++i){
      alphaHat(t,i-1) = roundDecimal_w( alpha(t,i-1) * scale[t], 8);
    }
  }
  return scale; 
}


//----- beta -----

//' Calculate beta
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param initalProbs A numeric vector of initial state probabilities. One probability per state
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
//' @param scale A vector of scaling values used to scale alpha
// [[Rcpp::export]]
NumericMatrix calculate_beta_w(NumericMatrix emissionProbs, NumericVector initalProbs, NumericVector positionBP, NumericMatrix transParam, NumericVector scale){
  int noSNPs = emissionProbs.nrow();
  int noStates = initalProbs.size();
  NumericMatrix transProbs(noStates,noStates);
  NumericMatrix beta(noSNPs,noStates);
  NumericMatrix betaHat(noSNPs,noStates);
  double betaSum;
  const int T = noSNPs-1;
  
  // Initialization
  for(int i = 1; i <= noStates; ++i){
    beta(T,i-1) = 1;
    betaHat(T,i-1) = roundDecimal_w( beta(T,i-1) * scale[T], 8);
  }
  
  // Induction
  for(int t = (T-1); t >= 0; t--){
    transProbs = transMatrix(positionBP[t+1] - positionBP[t], transParam);
    for(int i = 1; i <= noStates; ++i){
      betaSum = 0.0;
      for(int j = 1; j <= noStates; ++j){
        betaSum += transProbs(i-1, j-1) * emissionProbs(t+1,j-1) * betaHat(t+1,j-1);
      }
      beta(t,i-1) = betaSum;
      betaHat(t,i-1) = roundDecimal_w( beta(t,i-1) * scale[t], 8);
    }
  }
  return betaHat; 
}


//----- gamma -----

//' Calculate gamma
//' @param alpha A matrix of alpha values calculated from the forward algorithm
//' @param beta A matrix of beta values calculated from the backward algorithm
// [[Rcpp::export]]
NumericMatrix calculate_gamma_w(NumericMatrix alpha, NumericMatrix beta){
  int noSNPs = alpha.nrow();
  int noStates = alpha.ncol();
  NumericMatrix gamma(noSNPs, noStates);
  double alphaBetaSum;
  
  for(int t = 0; t < noSNPs; ++t){
    alphaBetaSum = 0;
    for(int i = 1; i <= noStates; ++i){
      gamma(t,i-1) = alpha(t,i-1) * beta(t,i-1);
      alphaBetaSum += gamma(t,i-1);
    }
    for(int i = 1; i <= noStates; ++i){
      gamma(t,i-1) = roundDecimal_w( gamma(t,i-1)/alphaBetaSum, 3) ;
    }
  }
  return gamma;
}


//----- Viterbi -----

//' Calculate the Viterbi sequence
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param initalProbs A numeric vector of initial state probabilities. One probability per state
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
// [[Rcpp::export]]
IntegerVector calculate_viterbi_w(NumericMatrix emissionProbs, NumericVector initalProbs, NumericVector positionBP, NumericMatrix transParam){
  int noSNPs = emissionProbs.nrow();
  int noStates = initalProbs.size();
  NumericMatrix transProbs(noStates,noStates);
  IntegerVector qStar(noSNPs);
  NumericVector deltaA(noStates);
  NumericMatrix delta(noSNPs, noStates);
  NumericMatrix psi(noSNPs, noStates);
  const int T = noSNPs-1;
  
  // Initialization
  for(int i = 1; i <= noStates; ++i){
    delta(0,i-1) = log(initalProbs[i-1]) + log(emissionProbs(0,i-1));
    psi(0,i-1) = 0;
  }
  
  // Recursion
  for(int t = 1; t < noSNPs; ++t){
    transProbs = transMatrix(positionBP[t] - positionBP[t-1], transParam);
    for(int j = 1; j <= noStates; ++j){
      for(int i = 1; i <= noStates; ++i){
        deltaA[i-1] = delta(t-1,i-1) + log(transProbs(i-1, j-1));
      }
      double maxNum = deltaA[0];
      int maxArg = 0;
      for(int k = 0; k < noStates; ++k){
        if(deltaA[k] > maxNum){
          maxNum = deltaA[k];
          maxArg = k; // could be problematic if there are more than 1 maximum
        }
      }
      delta(t,j-1) = maxNum + log(emissionProbs(t,j-1)) ;
      psi(t,j-1) = maxArg ; 
    } 
  }
  
  // Termination
  qStar[T] = 0;
  double probStar = delta(T,0);
  for(int k = 0; k < noStates; ++k){
    if(delta(T,k)>probStar){
      qStar[T] = k; // could be problematic if there are more than 1 maximum
      probStar = delta(T,k);
    }
  }
  
  // Path backtracking
  for(int t = (T-1); t >= 0; t--){
    qStar[t] = psi(t+1,qStar[t+1]); // Not sure if this will work without qstar being a constant
  }
  return qStar;
}


//----- Xi -----

//' Calculate the Xi
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
//' @param alpha A matrix of alpha values calculated from the forward algorithm
//' @param beta A matrix of alpha values calculated from the backward algorithm
// [[Rcpp::export]]
NumericMatrix calculate_xi(NumericMatrix emissionProbs, NumericVector positionBP, NumericMatrix transParam, NumericMatrix alpha, NumericMatrix beta){
  int noSNPs = alpha.nrow();
  int noStates = alpha.ncol();
  NumericMatrix transProbs(noStates,noStates);
  NumericMatrix xi(noSNPs, noStates*noStates);
  int xiCol;
  double xiSum;
  
  for(int t = 0; t < (noSNPs-1); ++t){
    xiCol = 0;
    xiSum = 0.0;
    transProbs = transMatrix(positionBP[t+1] - positionBP[t], transParam);
    for(int i = 1; i <= noStates; ++i){
      for(int j = 1; j <= noStates; ++j){
        xi(t,xiCol) = alpha(t,i-1) * transProbs(i-1, j-1) * emissionProbs(t+1,j-1) * beta(t+1,j-1);
        xiSum += xi(t,xiCol);
        xiCol += 1;
      }
    }
    xi(t,_) = xi(t,_)/xiSum;
  }
  return xi;
}


//----- re-estimate transition probabilities -----

//' Re-estimate the transition probabilities
//' @param xi A numeric matrix of xi values
//' @param gamma A numeric matrix of gamma values (posterior probabilities)
// [[Rcpp::export]]
NumericMatrix reest_transition(NumericMatrix xi, NumericMatrix gamma){
  int noSNPs = gamma.nrow();
  int noStates = gamma.ncol();
  NumericMatrix transProbs(noStates,noStates);
  double gammaSum;
  double xiSum;
  double transSum;
  int xiCol;
  double paramChange = 0.00000001;
  
  xiCol = 0;
  for(int i = 1; i <= noStates; ++i){
    gammaSum = 0.0;
    for(int t = 0; t < noSNPs; ++t){
      gammaSum += gamma(t,i-1);
    }
    if(gammaSum == 0) gammaSum = paramChange;
    for(int j = 1; j <= noStates; ++j){
      xiSum = 0.0;
      for(int t = 0; t < noSNPs; ++t){
        xiSum += xi(t,xiCol);
      }
      transProbs(i-1,j-1) = paramChange + (1 - paramChange) * xiSum/gammaSum;
      xiCol += 1;
    }
  }
  for(int i = 1; i <= noStates; ++i){
    transSum = 0.0;
    for(int j = 1; j <= noStates; ++j){
      transSum += transProbs(i-1,j-1);
    }
    if(transSum != 1) transProbs(i-1,_) = transProbs(i-1,_)/transSum;
  }
  return transProbs;
}


//----- log-likelihood -----

//' Calculate the log-likelihood
//' @param scale A vector of scaling values used to scale alpha
// [[Rcpp::export]]
double calculate_logLikelihood_w(NumericVector scale){
  int noSNPs = scale.size();
  double sumScale = 0;
  for(int t = 0; t < noSNPs; t++){
    sumScale += log(scale[t]);
  }
  return roundDecimal_w( -sumScale, 4);
}


//----- Baum-Welch -----

//' Calculate the Baum-Welch algorithm
//' @param emissionProbs A numeric matrix of pre-calculated emission probabilities for each state and SNP
//' @param initalProbs A numeric vector of initial state probabilities. One probability per state
//' @param positionBP A vector of genetic map positions in base-pairs
//' @param transParam A matrix of transition probabilities
// [[Rcpp::export]]
NumericMatrix calculate_BaumWelch(NumericMatrix emissionProbs, NumericVector initalProbs, NumericVector positionBP, NumericMatrix transParam){
  int noSNPs = emissionProbs.nrow();
  int noStates = initalProbs.size();
  NumericMatrix alphaBW(noSNPs,noStates);
  NumericMatrix betaBW(noSNPs,noStates);
  NumericMatrix gammaBW(noSNPs,noStates);
  NumericMatrix xiBW(noSNPs,noStates*noStates);
  NumericMatrix reestTransProbs(noStates,noStates);
  NumericVector scaleBW(noSNPs);
  double delta = 10;
  double prevlogLike;
  double newlogLike;
  double DELTA = 1.0;
  int iteration = 0;
  
  
  // Initial run
  alphaBW   = calculate_alpha_w(emissionProbs, initalProbs, positionBP, transParam);
  scaleBW   = calculate_scale_w(emissionProbs, initalProbs, positionBP, transParam);
  betaBW    = calculate_beta_w(emissionProbs, initalProbs, positionBP, transParam, scaleBW);
  gammaBW   = calculate_gamma_w(alphaBW, betaBW);
  xiBW      = calculate_xi(emissionProbs, positionBP, transParam, alphaBW, betaBW);
  prevlogLike = calculate_logLikelihood_w(scaleBW);
  
  while(delta > DELTA){
    reestTransProbs = reest_transition(xiBW, gammaBW);
    
    alphaBW  = calculate_alpha_w(emissionProbs, initalProbs, positionBP, reestTransProbs);
    scaleBW  = calculate_scale_w(emissionProbs, initalProbs, positionBP, reestTransProbs);
    betaBW   = calculate_beta_w(emissionProbs, initalProbs, positionBP, reestTransProbs, scaleBW);
    gammaBW  = calculate_gamma_w(alphaBW, betaBW);
    xiBW     = calculate_xi(emissionProbs, positionBP, reestTransProbs, alphaBW, betaBW);
    newlogLike = calculate_logLikelihood_w(scaleBW);
    
    delta = abs( newlogLike - prevlogLike) ;
    prevlogLike = newlogLike;
    iteration = iteration + 1;
    if(iteration > 10) delta = 0;
  }
  
  return reestTransProbs;
}

