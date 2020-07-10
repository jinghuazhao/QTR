setwd("/media/data/project/HUA_methylation_analysis_matie/")
getwd()

datafile <- 'rrbs_clean_data_matie/predictedMeth_beta.RDS'
rangefile <- 'rrbs_clean_data_matie/predictedMeth_range.RDS'
betav <- readRDS(datafile)
rangev <- readRDS(rangefile)
#qualityb <- (rowMeans(betav, na.rm = TRUE) < 0.05) | (rowSums(is.na(betav)) > 4)
#betav <- betav[!qualityb, ]
offset <- 1e-5
betav[betav == 0] <- offset
betav[betav == 1] <- 1 - offset
mv <- log2(betav / (1 - betav))

saveRDS(mv, gsub("_beta.RDS", "_m_quality.RDS", datafile))
rangev <- as.data.frame(rangev)
#rangev <- rangev[!qualityb, ]
saveRDS(rangev, gsub("_range.RDS", "_range_quality.RDS", rangefile))

library(plyr)
library(doMC)
doMC::registerDoMC(cores = 14)
mv <- alply(.data = mv, .margins = 1, .fun = function(x) {
				x[which(is.na(x))] <- median(x, na.rm = TRUE)
				return(x)
}, .parallel = TRUE)
mv <- do.call(rbind, mv)

options(stringsAsFactors = FALSE)
mvtable <- as.data.frame(mv)
mvtable <- data.frame(ID = rownames(mvtable), mvtable)

datafile <- "rrbs_clean_data_matie/refactor_mv.txt"
write.table(mvtable, datafile, quote = FALSE, sep = '\t')

source("scripts_matie/refactor_modify.R")
k = 5
refactor_obj <- refactor(datafile,k)


datafile <- 'rrbs_clean_data_matie/predictedMeth_m.RDS'
saveRDS(mv, datafile)
rm(betav, mv, qualityb, offset)

refactor_obj_name <- "rrbs_clean_data_matie/refactor_obj.RDS"
saveRDS(refactor_obj, refactor_obj_name)
