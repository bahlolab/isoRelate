[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.197857.svg)](https://doi.org/10.5281/zenodo.197857)


# isoRelate: software for inferring pairwise identity by descent in haploid species


## Features

* Estimates the proportion of genome shared IBD between pairs of isolates
* Detects genomic regions that are identical by descent (IBD) between pairs of isolates
* Allows IBD detection in isolates with multiple infections
* Identifies genomic loci under positive selection
* Creates networks of related isolates
* Includes multiple graphical functions to explore the results


## How to install isoRelate

isoRelate is currently available to install as a development version from Github:

```{r}
# first install the devtools package
install.packages("devtools")

# install isoRelate using devtools packages
devtools::install_github("bahlolab/isoRelate")
```


## How to use isoRelate

See the introduction vignette for details.


## Reference

Henden L, Lee S, Mueller I, Barry A, Bahlo M. 2016. Detecting selection signals in *Plasmodium falciparum* using identity-by-descent analysis. *bioRxiv*. doi: http://dx.doi.org/10.1101/088039
