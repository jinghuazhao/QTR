# 16-8-2019 JHZ

export imp2=zhugu_GWAS_20181119/Hearing/imp2
seq 22|awk -vimp2=$imp2 -vp=merge. '{print imp2 "/" p $1}' > merge-list
plink --merge-list merge-list --extract w3w4.rsid --make-bed --out merge

