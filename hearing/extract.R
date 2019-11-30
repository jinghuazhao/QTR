# 16-8-2019 JHZ

dz <- read.csv("twin hearing  DZ.csv", as.is=TRUE)
fam <- read.table("imp2/merge.1.fam", as.is=TRUE, col.names=c("pid","iid","fid","mid","gender","aff"))

dz_fam <- merge(dz,fam, by.x="no",by.y="iid")
ids <- with(dz_fam, cbind(pid,no))
write.table(ids, file="w3w4.ids", col.names=FALSE, row.names=FALSE, quote=FALSE)
head(dz_fam)

d <- within(dz_fam, {
   bh1 <- left1
   if (right1 < left1) bh1=right1
   bh1 <- gap::invnormal(sqrt(bh1))
})
pheno <- d[c("pid","no","bh1")]
write.table(pheno, file="w3w4.bh1", col.names=FALSE, row.names=FALSE, quote=FALSE)
covar <- d[c("pid","no","sex","age","education")]
write.table(covar, file="w3w4.cov", col.names=FALSE, row.names=FALSE, quote=FALSE)

pheno <- d[c("bh1")]
write.table(pheno, file="gemma.bh1", col.names=FALSE, row.names=FALSE, quote=FALSE)
covar <- d[c("sex","age","education")]
write.table(covar, file="gemma.cov", col.names=FALSE, row.names=FALSE, quote=FALSE)
