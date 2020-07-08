#!/usr/bin/bash

mkdir src
ls BiSeq/* | xargs -l basename -s .gz | parallel -j10 -C' ' 'zgrep -v "_" BiSeq/{}.gz | gzip -f > src/{}.gz'
