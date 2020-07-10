#!/usr/bin/bash

function setup()
{
  export list=$(ls BiSeq/* | xargs -l basename -s .gz)
  export chrs=$(echo $(seq 22) X Y M)
  export out=partitioned

  mkdir src
  parallel --env out -j10 -C' ' 'zcat BiSeq/{1}.gz | awk -v chr=chr{2} "\$1==chr" | gzip -f > ${out}/{1}_chr{2}.gz' ::: ${list} ::: ${chrs}
  parallel --env out -j10 -C' ' 'zgrep "_" BiSeq/{1}.gz | gzip -f > ${out}/{}_Z.gz' ::: ${list}
}

export suffix=lmer

export chrs=$(echo $(seq 22) X Y M Z)
for chromosome in ${chrs}
do
  export chr=chr${chromosome}
  export srcdir=scripts_${suffix}
  echo --- running ${suffix} chromosome ${chromosomes} ---
  R --no-save -q < ${srcdir}/1_BiSeq.R > 1.log
  R --no-save -q < ${srcdir}/2_format_methylation_data.R > 2.log
  R --no-save -q < ${srcdir}/3_quality_control_and_refactor.R > 3.log
  R --no-save -q < ${srcdir}/4_analysis_lmer.R > 4.log
done
R --no-save -q < ${srcdir}/5_read_result.R
R --no-save -q < ${srcdir}/6_get_annotation.R
R --no-save -q < ${srcdir}/7_combined-pvalues.sh
R --no-save -q < ${srcdir}/8_sigplot.R
                 

