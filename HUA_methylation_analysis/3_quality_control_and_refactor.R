# setwd("/media/data/project/HUA_methylation_analysis/")
getwd()

options(stringsAsFactors = FALSE)

datafile <- 'rrbs_clean_data/predictedMeth_beta.RDS'
rangefile <- 'rrbs_clean_data/predictedMeth_range.RDS'

betav <- readRDS(datafile)
rangev <- readRDS(rangefile)
# I am not sure whether the following two scripts should be run or not
qualityb <- (rowMeans(betav, na.rm = TRUE) < 0.05) | (rowSums(is.na(betav)) > 4)
betav <- betav[!qualityb, ]
offset <- 1e-5
betav[betav == 0] <- offset
betav[betav == 1] <- 1 - offset
mv <- log2(betav / (1 - betav))
saveRDS(mv, gsub("_beta.RDS", "_m_quality.RDS", datafile))

# I am not sure whether the following script should be run or not
rangev <- rangev[!qualityb, ]
saveRDS(rangev, gsub("_range.RDS", "_range_quality.RDS", rangefile))

library(plyr)
library(doMC)
doMC::registerDoMC(cores = 14)
colnamesmv <- colnames(mv)
rownamesmv <- rownames(mv)
mv <- alply(.data = mv, .margins = 1, .fun = function(x) {x[which(is.na(x))] <- median(x, na.rm = TRUE)	return(x)}, .parallel = TRUE)
mv <- do.call(rbind, mv)
mv <- as.data.frame(mv)
rownames(mv) <- rownamesmv
colnames(mv) <- colnamesmv

saveRDS(mv, 'rrbs_clean_data/predictedMeth_m.RDS')

options(stringsAsFactors = FALSE)
mvtable <- mv
mvtable <- data.frame(ID = rownames(mvtable), mvtable)

refactor_datafile <- "rrbs_clean_data/refactor_mv.txt"
write.table(mvtable, refactor_datafile, quote = FALSE, sep = '\t')

source("refactor_modify.R")
k = 5
refactor_obj <- refactor(refactor_datafile,k, out='rrbs_clean_data/refactor')
PCs <- as.data.frame(refactor_obj$standard_pca)
rownames(PCs) <- colnamesmv

saveRDS(PCs,  "rrbs_clean_data/refactor_obj.RDS")
