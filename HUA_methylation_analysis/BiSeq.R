# setwd("/media/data/project/HUA_methylation_analysis/")
getwd()

options(stringsAsFactors = FALSE)
require(parallel)
require(BiSeq)
require(openxlsx)
# bm_phen <- read.xlsx("SUA_sample_information.xlsx")
bm_phen <- read.csv("SUA_sample_information.csv")

list_files <- list.files("src", full.names = TRUE)
N <- length(list_files)
y <- setNames(list_files,paste0("sample",1:N))

readBatch <- function(set)
# sets of samples
{
  bm_files <- as.character(y[paste0("sample",set)])
  bm_ids <- sub(".*([0-9]{4}[a-c]?).*", "\\1", bm_files, perl = TRUE)
  filtered_ids <- which(! bm_ids %in% c("3111","3112","4149", "4150"))
  bm_files <- bm_files[filtered_ids]
  bm_ids <- bm_ids[filtered_ids]
  bm_group <- bm_phen[match(bm_ids, bm_phen$no), "SUA_level"]
  bm_group <- as.factor(bm_group)
  rrbs <- readBismark(bm_files, colData = DataFrame(row.names = bm_ids, group = bm_group))
}

dir.create("rrbs_clean_data_lmer")
rrbs1 <- readBatch(1:45)
rrbs2 <- readBatch(46:90)
rrbs3 <- readBatch(91:134)
save(rrbs1,rrbs2,rrbs3,file="rrbs.rda")
rrbs12 <- combine(rrbs1,rrbs2)
rrbs <- combine(rrbs12,rrbs3)
saveRDS(rrbs, "rrbs_clean_data_lmer/rrbs.RDS")

# pdf('rrbs_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png("rrbs_clean_data_lmer/rrbs_covBoxplots.png", width = 900, height = 480, units = 'px', pointsize = 12)
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
saveRDS(rrbs.clust.lim, 'rrbs_clean_data_lmer/rrbs.clust.lim.RDS')

# pdf('rrbs.clust.lim_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png('rrbs_clean_data_lmer/rrbs.clust.lim_covBoxplots.png', width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs.clust.lim, col = "cornflowerblue", las = 2)
dev.off()

predictedMeth <- predictMeth(object = rrbs.clust.lim, mc.cores = 1)
saveRDS(predictedMeth, 'rrbs_clean_data_lmer/predictedMeth.RDS')

# bs_range <- as.data.frame(rowRanges(predictedMeth))
# bs_methylation <- methLevel(predictedMeth)