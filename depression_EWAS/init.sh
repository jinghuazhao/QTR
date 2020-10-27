#!/usr/bin/bash

export src=BiSeq
export prefix=CpG_context_Human_
export suffix=_good_1_val_1_bismark_bt2_pe.bismark.cov.gz
export dest=BiSeq_cleaned

if [ ! -d ${dest} ]; then mkdir ${dest}; fi

for n in $(ls ${src} | sed 's/'"${prefix}"'//;s/'"${suffix}"'//')
do
  echo ${src}/${prefix}${n}${suffix}
  zcat ${src}/${prefix}${n}${suffix} | awk '$1~/^chr[0-9]$|^chr[0-9][0-9]$|chrX|chrY|chrM/' | gzip -f > ${dest}/${n}.gz
done

# more sophisticated
# ls ${src} | xargs -I {} sh -c "basename {} ${suffix} | sed 's/"${prefix}"//'"
