#!/usr/bin/bash

function partition()
{
  export list=$(ls BiSeq_cleaned/* | xargs -l basename -s .gz)
  export chrs=$(echo $(seq 22) X Y M)
  export out=partitioned

  if [ ! -d ${out} ]; then mkdir ${out}; fi
  parallel --env out -j10 -C' ' 'zcat BiSeq_cleaned/{1}.gz | awk -v chr=chr{2} "\$1==chr" | gzip -f > ${out}/{1}_chr{2}.gz' ::: ${list} ::: ${chrs}

  export suffix=lmer
  export chrs=$(echo $(seq 22) X Y M)
  if [ ! -d rrbs_clean_data_${suffix} ]; then mkdir rrbs_clean_data_${suffix}; fi
  for chromosome in ${chrs}
  do
    export chr=chr${chromosome}
    export scripts_dir=depression_scripts_chromosomes_${suffix}
    echo --- running ${suffix} chromosome ${chromosome} ---
    R --no-save -q < ${scripts_dir}/1_BiSeq.R > 1.log
    R --no-save -q < ${scripts_dir}/2_format_methylation_data.R > 2.log
    R --no-save -q < ${scripts_dir}/3_quality_control_and_refactor.R > 3.log
    R --no-save -q < ${scripts_dir}/4_analysis_lmer.R > 4.log
  done
}

function all()
{
  R --no-save -q < 1_BiSeq.R > 1.log
  R --no-save -q < 2_format_methylation_data.R > 2.log
  R --no-save -q < 3_quality_control_and_refactor.R > 3.log
  R --no-save -q < 4_analysis_lmer.R > 4.log
}

# below for both partition and all
R --no-save -q < ${scripts_dir}/5_read_result.R
R --no-save -q < ${scripts_dir}/6_get_annotation.R
7_combined-pvalues.sh
R --no-save -q < ${scripts_dir}/8_sigplot.R
