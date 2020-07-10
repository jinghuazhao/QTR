# setwd("/media/data/project/HUA_methylation_analysis/")
getwd()

options(stringsAsFactors = FALSE)
require(parallel)

library(BiSeq)

chr <- Sys.getenv("chr")
chr <- Sys.getenv("suffix")
outdir <- paste0("rrbs_clean_data_",suffix,"/")
rrbs_file <- paste0(outdir,"rrbs-",chr,".RDS")
rrbs.clust.lim_file <- paste0(outdir,"rrbs.clust.lim-",chr,".RDS")
predictedMeth_file <- paste0(outdir,"predictedMeth-",chr,".RDS")
rrbs_covBoxPlots_file <- paste0(outdir,rrbs_covBoxPlots-",chr.".png")
rrbs.clust.lim_covBoxplots_file <- paste0(outdir,"rrbs.clust.lim_covBoxplots-",chr,"png")

list_files <- list.files("BiSeq", full.names = TRUE)
bm_ids <- sub(".*([0-9]{4}[a-c]?).*", "\\1", list_files, perl = TRUE)
filtered_ids <- which(! bm_ids %in% c("3003", "3004", "3111", "3112"))
bm_files <- paste0("partitioned/",sub("cov",paste0("cov_",chr),basename(list_files[filtered_ids])))
bm_ids <- bm_ids[filtered_ids]

library(openxlsx)
bm_phen <- read.xlsx("SUA_sample_information.xlsx")
bm_group <- bm_phen[match(bm_ids, bm_phen$no), "SUA_level"]
bm_group <- as.factor(bm_group)

rrbs <- readBismark(bm_files, colData = DataFrame(row.names = bm_ids, group = bm_group))
dir.create(outdir)
saveRDS(rrbs, rrbs_file)

# pdf('rrbs_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png(paste0(rrbs_covBoxplots_file, width = 900, height = 480, units = 'px', pointsize = 12)
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
saveRDS(rrbs.clust.lim, rrbs.clust.lim_file)

# pdf('rrbs.clust.lim_covBoxplots.pdf', paper = 'a4r', width = 0, height = 0)
png(rrbs.clust.lim_covBoxplots_file, width = 900, height = 480, units = 'px', pointsize = 12)
covBoxplots(rrbs.clust.lim, col = "cornflowerblue", las = 2)
dev.off()

predictedMeth <- predictMeth(object = rrbs.clust.lim, mc.cores = 1)
saveRDS(predictedMeth, predictedMeth_file)

# bs_range <- as.data.frame(rowRanges(predictedMeth))
# bs_methylation <- methLevel(predictedMeth)
