# 6-8-2019 JHZ

r1 <- read.csv("XJJ0043-W3-1.csv",as.is=TRUE,skip=2)
dim(r1)
r2 <- read.csv("XJJ0043-W3-2.csv",as.is=TRUE,skip=2)
dim(r2)
r3 <- read.csv("XJJ0043-W3-96.csv",as.is=TRUE,skip=2)
dim(r3)
xin <- read.csv("XJJ0043-W3-96-XIN.csv",as.is=TRUE,skip=2)
dim(xin)
r <- rbind(r1[,c(1,3,4)],r2[,c(1,3,4)],r3[,c(1,3,4)],xin[,c(1,3,4)])
dim(r)
nrow(r)/30
write.table("rawdata.txt",row.names=FALSE,quote=FALSE)
