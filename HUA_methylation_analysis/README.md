# scripts_lmer/ and scripts_matie/

* 1_BiSeq.R in both directories read a .csv phenotype file instead of .xlsx as well as src/ directory instead of BiSeq/. The former makes it easier to compromise the triplet with suffexes as 1a2a,1b3b,2c3c rather than a sequence of numbers; the latter contains those without "_" in the chromosome name through [zgrep.sh](zgrep.sh).
* rrbs_clean_data_lmer/ and rrbs_clean_data_matie/ are the working directories, respectively.

# BiSeq

* [Source](https://www.bioconductor.org/packages/release/bioc/src/contrib/BiSeq_1.28.0.tar.gz)
* [Revised version](BiSeq_1.28.1.tar.gz): this can be seen with `tar xvfz BiSeq_1.28.1.tar.gz` through `R/readBismark.R` mainly the following two statements,
```r
  tReads <- matrix(numeric(length = as.numeric(length(fData)) * as.numeric(length(methData))), nrow=length(fData))
  mReads <- matrix(numeric(length = as.numeric(length(fData)) * as.numeric(length(methData))), nrow=length(fData))
```
where `integer` is changed to `numeric` -- we don't have integer overflow but may have memory problem for a huge request.

# Tutorial on DSS

* [Paper](https://doi.org/10.1007/s40484-019-0183-8)
* [R code](https://static-content.springer.com/esm/art%3A10.1007%2Fs40484-019-0183-8/MediaObjects/40484_2019_183_MOESM2_ESM.zip)
