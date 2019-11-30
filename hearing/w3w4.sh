# 12-8-2019 JHZ


# list of SNPs

head -1 w3w4.csv | \
sed 's/,/\n/g' | \
awk 'NR > 1{print $1}' > w3w4.rsid

# SNP annotation through PhenoScanner

R --no-save -q <<END
# devtools::install("phenoscanner/phenoscanner")
  library(phenoscanner)
  w3w4 <- scan("w3w4.rsid", what="")
  length(w3w4)
  w3w4_ps <- phenoscanner(w3w4)
  names(w3w4_ps)
  attach(w3w4_ps)
  map <- within(snps,{z <- 0; chr <- as.numeric(chr)})[c("chr","rsid","z","pos_hg19","a1","a2")]
  detach(w3w4_ps)
  write.table(map,file="w3w4.snps", col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
END

## Build bed/bim/fam files from map/ped files

echo "rs258751 needs to be added"
join -v2 <(cut -f2 w3w4.snps | sort -k1,1) <(sort w3w4.rsid)
(
  head -10 w3w4.snps
  echo -e "5\trs258751\t0\t142662280\tA\tG"
  tail -n 49 w3w4.snps
) > w3w4.map

awk 'NR>1' w3w4.csv | \
awk -vFS="," -vOFS=" " '
{
  id = sprintf("%d %d %d %d %d %d", NR, $1, 0, 0, 1, 1)
  printf id
  for(i=2; i<=NF; i++)
  {
    if ($i == "") g="0 0";
    else {
      a1=substr($i,1,1)
      a2=substr($i,2,1)
      if (a2=="") a2=a1
      g=a1 " " a2
    }
    printf OFS g
  }
  printf "\n"
}' > w3w4.ped

echo "entries per line"
awk '{print NF}' w3w4.ped | sort | uniq

./plink --file w3w4 --make-bed --out w3w4

# Handling of phenotypes/covariates via SPSS (example)

R --no-save -q <<END
library(foreign)
d <-read.spss("zhugu_GWAS_20181119/Hearing/mxjobs/pta5/duan2_LR_u1_u.sav")
d <- as.data.frame(d)
head(d)
# simulation
pheno <- matrix(rnorm(868*2),ncol=2)
covar <- rnorm(868)
# write them as text files
write.table(pheno,file="w3w4.pheno",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(covar,file="w3w4.covar",col.names=FALSE,row.names=FALSE,quote=FALSE)
END

# lm/lmm by GEMMA

export gemma=GEMMA-0.98.1/bin/gemma

## 

GEMMA-0.98.1/bin/gemma -bfile w3w4 -gk 1 -o w3w4
$gemma -bfile w3w4 -lmm 1 -o w3w4 -p w3w4.pheno -n 1 -c w3w4.covar -k output/w3w4.cXX.txt -o w3w4-1
