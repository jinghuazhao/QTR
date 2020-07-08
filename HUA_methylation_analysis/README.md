## The BiSeq problem

It becomes apparent that beyond certain number of samples, BiSeq suddenly has integer overflow/request enormous amount of memory.

## A new implementation with BiSeq

1. Work in batches and then combine them.

[BiSeq.R](BiSeq.R) replaces 1_BiSeq.R in both experiments below. It handles the triplet with suffexes as 1a2a,1b3b,2c3c rather than a sequence of numbers, as readily visible from *SUA_sample_information.csv*.

Working directories for programs in **scripts_lmer/** and **scripts_matie/** now use **rrbs_clean_data_lmer/** and **rrbs_clean_data_matie/** are the working directories, respectively.

2. Amendment to R/BiSeq package.

* [BiSeq 1.28.0](https://www.bioconductor.org/packages/release/bioc/src/contrib/BiSeq_1.28.0.tar.gz) is hosted at Bioconductor.
* [BiSeq 1.28.1](BiSeq_1.28.1.tar.gz) contains minor changes to the following two statements in *R/readBismark.R* of the original package,
```r
# from
  tReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
  mReads <- matrix(integer(length = length(fData) * length(methData)), nrow=length(fData))
# to
  tReads <- as.integer(matrix(numeric(length = as.numeric(length(fData)) * as.numeric(length(methData))), nrow=length(fData)))
  mReads <- as.integer(matrix(numeric(length = as.numeric(length(fData)) * as.numeric(length(methData))), nrow=length(fData)))
```
where `integer` is changed to `numeric` -- we don't have integer overflow but may have memory problem for a huge request (>250G); `as.integer()` is required since

## Legacy experiments

These efforts were incremental in nature and did not resolve the issue -- they are kept only historical reasons.

1. An attempt was to drop CpG sites from BiSeq/ into src/ through [zgrep.sh](zgrep.sh).

> ... class “BSraw” object: The totalReads and methReads matrices of an BSraw object must contain integer data.

3. DSS as a possible alternative.

* [Paper](https://doi.org/10.1007/s40484-019-0183-8)
* [R code](https://static-content.springer.com/esm/art%3A10.1007%2Fs40484-019-0183-8/MediaObjects/40484_2019_183_MOESM2_ESM.zip)
