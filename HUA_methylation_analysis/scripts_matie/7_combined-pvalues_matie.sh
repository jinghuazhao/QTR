cd combined-pvalues
comb-p acf -d 1:500:50 -c 4 data/pvals.bed > data/acf.txt
comb-p slk --acf data/acf.txt -c 4 data/pvals.bed > data/pvals.acf.bed
comb-p peaks --dist 500 --seed 0.1 data/pvals.acf.bed > data/pvals.regions.bed
comb-p region_p -p data/pvals.bed -r data/pvals.regions.bed -s 50 -c 4 > data/regions.sig.bed

