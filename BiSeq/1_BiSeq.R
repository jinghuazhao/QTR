#the "setwd" was added by myself, it is important to set the working directory firstly
setwd("/home/w/working/20190326_认知功能甲基化分析")
getwd()
options(stringsAsFactors = FALSE)
#parallel can't be loaded
#install.packages("parallel")
require(parallel)

library(BiSeq)

bm_files <- list.files("BiSeq", full.names = TRUE)
bm_ids <- sub(".*([0-9]{4}).*", "\\1", bm_files, perl = TRUE)
filtered_ids <- which(! bm_ids %in% c("3003", "3004"))
bm_files <- bm_files[filtered_ids]
bm_ids <- bm_ids[filtered_ids]


library(openxlsx)
bm_phen <- read.xlsx("20190326_sample_information(copy).xlsx")
bm_group <- bm_phen[match(bm_ids, bm_phen$no), "cognitive_function_score_group"]
bm_group <- as.factor(bm_group)

#problem: "could not find function of "readBismark"
rrbs <- readBismark(bm_files, colData = DataFrame(row.names = bm_ids, group = bm_group))
#the directory of "rrbs_clean_data1" was created by myself
dir.create("rrbs_clean_data1")
saveRDS(rrbs, "rrbs_clean_data1/rrbs.RDS")

# pdf('rrbs_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png("rrbs_clean_data1/rrbs_covBoxplots.png", width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs, col = "cornflowerblue", las = 2)
dev.off()

  rrbs.clust.unlim <- clusterSites(object = rrbs,
                                 groups = colData(rrbs)$group,
                                 perc.samples = 4/5,
                                 min.sites = 20,
                                 max.dist = 100,
                                 mc.cores = 1)

ind.cov <- totalReads(rrbs.clust.unlim) > 0
quant <- quantile(totalReads(rrbs.clust.unlim)[ind.cov], 0.9)
quant
rrbs.clust.lim <- limitCov(rrbs.clust.unlim, maxCov = quant)
saveRDS(rrbs.clust.lim, 'rrbs_clean_data1/rrbs.clust.lim.RDS')

# pdf('rrbs.clust.lim_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png('rrbs_clean_data1/rrbs.clust.lim_covBoxplots.png', width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs.clust.lim, col = "cornflowerblue", las = 2)
dev.off()

predictedMeth <- predictMeth(object = rrbs.clust.lim, mc.cores = 1)
saveRDS(predictedMeth, 'rrbs_clean_data1/predictedMeth.RDS')

# bs_range <- as.data.frame(rowRanges(predictedMeth))
# bs_methylation <- methLevel(predictedMeth)
