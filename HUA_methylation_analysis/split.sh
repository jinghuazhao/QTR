#!/usr/bin/bash

export list=$(ls BiSeq/* | xargs -l basename -s .gz)
export chrs=$(echo $(seq 22) X Y M)

mkdir src
parallel -j10 -C' ' 'zcat BiSeq/{1}.gz | awk -v chr=chr{2} '\$1==chr' | gzip -f > src/{1}_chr{2}.gz' ::: ${list} ::: ${chrs}
parallel -j10 -C' ' 'zgrep "_" BiSeq/{1}.gz | gzip -f > src/{}_.gz' ::: ${list}
