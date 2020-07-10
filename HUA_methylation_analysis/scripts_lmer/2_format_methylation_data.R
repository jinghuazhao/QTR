# setwd("/media/data/project/HUA_methylation_analysis/")
getwd()

options(stringsAsFactors = FALSE)

chr <- Sys.getenv("chr")
suffix <- Sys.getenv("suffix")
outdir <- paste0("rrbs_clean_data_",suffix,"/")
rrbs.clust.lim_file <- paste0(outdir,"rrbs.clust.lim-",chr,".RDS")
predictedMeth_file <- paste0(outdir,"predictedMeth-",chr,".RDS")

bsraw_matrix <- function(filename) {

	require(BiSeq)

	bsraw <- readRDS(filename)

	mgranges <- rowRanges(bsraw)
        mgranges <- as.data.frame(mgranges)
	rownames(mgranges) <- paste0('r', 1:nrow(mgranges))

	## move first row to the last row
	mgranges_new <- rbind.data.frame(mgranges[2:nrow(mgranges), ], mgranges[1, ])
	idx_keep <- rownames(mgranges)[!((mgranges$seqnames == mgranges_new$seqnames) & (mgranges$start == (mgranges_new$start - 1)))]
	mgranges <- mgranges[idx_keep, ]

	mtotal <- totalReads(bsraw)
	mmeth <- methReads(bsraw)

	mbeta <- mmeth/mtotal
	mbeta[which(is.nan(mbeta), arr.ind = TRUE)] <- NA
	rownames(mbeta) <- paste0('r', 1:nrow(mbeta))

	## filter
	mbeta <- mbeta[idx_keep, ]

	file_range <- gsub('.RDS', '_range.RDS', filename)
	file_beta <- gsub('.RDS', '_beta.RDS', filename)

	saveRDS(mgranges, file_range)
	saveRDS(mbeta, file_beta)
}

bsraw_matrix(rrbs.clust.lim_file)

bsrel_matrix <- function(filename) {

	require(BiSeq)

	bsrel <- readRDS(filename)

	mgranges <- rowRanges(bsrel)
	mgranges <- as.data.frame(mgranges)
	rownames(mgranges) <- paste0('r', 1:nrow(mgranges))

	## move first row to the last row
	mgranges_new <- rbind.data.frame(mgranges[2:nrow(mgranges), ], mgranges[1, ])
	idx_keep <- rownames(mgranges)[!((mgranges$seqnames == mgranges_new$seqnames) & (mgranges$start == (mgranges_new$start - 1)))]
	mgranges <- mgranges[idx_keep, ]

	mbeta <- methLevel(bsrel)
	rownames(mbeta) <- paste0('r', 1:nrow(mbeta))

	## filter
	mbeta <- mbeta[idx_keep, ]

	file_range <- gsub('.RDS', '_range.RDS', filename)
	file_beta <- gsub('.RDS', '_beta.RDS', filename)

	saveRDS(mgranges, file_range)
	saveRDS(mbeta, file_beta)
}

bsrel_matrix(predictedMeth_file)
