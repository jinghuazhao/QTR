ls src/* | xargs -l basename -s .gz | parallel -j10 -C' ' 'zgrep -v "_" src/{}.gz | gzip -f > BiSeq/{}.gz'
