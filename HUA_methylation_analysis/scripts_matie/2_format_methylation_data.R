# setwd("/media/data/project/HUA_methylation_analysis_matie/")
getwd()
options(stringsAsFactors = FALSE)

bsraw_matrix <- function(filename) {

	require(BiSeq)

	bsraw <- readRDS(filename)

	mgranges <- rowRanges(bsraw)

	mtotal <- totalReads(bsraw)
	mmeth <- methReads(bsraw)

	mbeta <- mmeth/mtotal
	mbeta[which(is.nan(mbeta), arr.ind = TRUE)] <- NA
	rownames(mbeta) <- paste0('r', 1:nrow(mbeta))

	file_range <- gsub('.RDS', '_range.RDS', filename)
	file_beta <- gsub('.RDS', '_beta.RDS', filename)

	saveRDS(mgranges, file_range)
	saveRDS(mbeta, file_beta)
}

bsraw_matrix('rrbs_clean_data_matie/rrbs.clust.lim.RDS')


bsrel_matrix <- function(filename) {

	require(BiSeq)

	bsrel <- readRDS(filename)

	mgranges <- rowRanges(bsrel)

	mbeta <- methLevel(bsrel)
	rownames(mbeta) <- paste0('r', 1:nrow(mbeta))

	file_range <- gsub('.RDS', '_range.RDS', filename)
	file_beta <- gsub('.RDS', '_beta.RDS', filename)

	saveRDS(mgranges, file_range)
	saveRDS(mbeta, file_beta)
}

bsrel_matrix('rrbs_clean_data_matie/predictedMeth.RDS')
