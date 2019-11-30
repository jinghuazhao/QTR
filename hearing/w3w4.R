# 6-8-2019 JHZ

Sys.setlocale("LC_ALL", "C") 
d3 <-"SNP最终结果/SNP检测报告-QDOE2018H10038-W3/SNP检测结果/SNP检测结果/"
f3 <- "PlateData-XJJ0043-W3-868_tr.csv"
w3 <- read.table(paste0(d3,f3),sep=",",header=TRUE)
names(w3)[1:6] <- c("plat11","plat12","plat13","Sample.Id","name1","PE1")
dim(w3)
head(w3)
d4 <-"SNP最终结果/SNP检测报告-QDOE2018H10038-W4/SNP检测结果/SNP检测结果/"
f4 <- "XJJ0043-W4-868_tr.csv"
w4 <- read.table(paste0(d4,f4),sep=",",header=TRUE)
names(w4)[1:6] <- c("plat21","plat22","plat23","Sample.Id","name2","PE2")
dim(w4)
head(w4)
m1 <- c("plat11","plat12","plat13","name1","PE1")
m2 <- c("plat21","plat22","plat23","name2","PE2")
w3w4 <- merge(w3[setdiff(names(w3),m1)],w4[setdiff(names(w4),m2)],by="Sample.Id")
write.table(w3w4,file="w3w4.csv",sep=",",row.names=FALSE,quote=FALSE)
