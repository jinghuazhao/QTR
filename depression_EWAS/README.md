## NOTE

Script [init.sh](init.sh) is used to generate files for chr1-22,M,X,Y only at BiSeq_cleaned.

### Whole-genome version

It is then necessary to change the directory name BiSeq to BiSeq_cleaned in `1_BiSeq.R`.

### Partition by chromosome

This is achieved with parition function inside [BiSeq.sh](BiSeq.sh) and necessary change in `depression_scripts_chromosomes_lmer/` for
* `1_BiSeq.R`
* `3_quality_control_and_refactor.R`.
* `8_sigplot.R` as there was a bug with directory name.
