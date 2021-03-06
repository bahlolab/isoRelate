% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{transitionProbDD}
\alias{transitionProbDD}
\title{The transition probabilities for 2 diploid chromosomes}
\usage{
transitionProbDD(omega_0, omega_1, omega_2, meiosis, dist_cM, ibd_current,
  ibd_previous)
}
\arguments{
\item{omega_0}{The probability of sharing 0 alleles IBD}

\item{omega_1}{The probability of sharing 1 allele IBD}

\item{omega_2}{The probability of sharing 2 alleles IBD}

\item{meiosis}{The number of meiosis separating the two isoaltes}

\item{dist_cM}{The genetic map distance (cM) between SNP i and SNP j}

\item{ibd_current}{The IBD state of SNP j}

\item{ibd_previous}{The IBD state of SNP i}
}
\description{
The transition probabilities for 2 diploid chromosomes
}
