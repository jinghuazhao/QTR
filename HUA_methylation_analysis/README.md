## The BiSeq problem

It becomes apparent that beyond certain number of samples, BiSeq suddenly has integer overflow/requests enormous amount of memory.

## A new implementation with BiSeq

1. Work in batches and then combine them.

[BiSeq.R](BiSeq.R) replaces 1_BiSeq.R in both experiments on *lmer* and *matie* below. It handles the triplet with suffexes as 1a2a,1b3b,2c3c rather than a sequence of numbers, as readily seen from *SUA_sample_information.csv*.

Working directories for programs in **scripts_lmer/** and **scripts_matie/** are now **rrbs_clean_data_lmer/** and **rrbs_clean_data_matie/**, respectively.

Unfortunately, the `combine` method of BSraw remains crashed. At least, it illustrates three batches here which could be saved, loaded again from a new session of R and retry, which turned to be successful.

2. Amendment to R/BiSeq package.

* [BiSeq 1.28.0](https://www.bioconductor.org/packages/release/bioc/src/contrib/BiSeq_1.28.0.tar.gz) is hosted at Bioconductor.
* [BiSeq 1.28.1](BiSeq_1.28.1.tar.gz) has the following changes in *R/readBismark.R* from the package above,
```r
# from
  tReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
  mReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
# to
  tReads <- matrix(0L, nrow=length(fData), ncol=length(methData))
  mReads <- matrix(0L, nrow=length(fData), ncol=length(methData))
```
then we don't have integer overflow on dimensions of the integer matrices but may have memory problem for a huge request (>250G).

*R/methods-BSraw.R* involves `length = nr*nc` and can be dealt with similarly.

## Possible way forward

Can we consider analysis by chromosome instead of sample? A data partition script is implemented as [split.sh]()(split.sh).

## Related work

Feng, H. and H. Wu (2019). "Differential methylation analysis for bisulfite sequencing using DSS." Quantitative Biology 7(4): 327-334, 
https://doi.org/10.1007/s40484-019-0183-8, [R code](https://static-content.springer.com/esm/art%3A10.1007%2Fs40484-019-0183-8/MediaObjects/40484_2019_183_MOESM2_ESM.zip).

## Acknowledgements

The use of `combine` method from BSraw class concurs with suggestion from the BiSeq package author, Katja Hebestreit <katja.hebestreit@gmail.com> as well as hint over analysis by chromosome.
