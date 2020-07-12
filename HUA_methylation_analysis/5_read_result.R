# setwd("/media/data/project/HUA_methylation_analysis_lmer/")
getwd()

options(stringsAsFactors = FALSE)

resultname <- 'rrbs_clean_data/lmer.RDS'

result <- readRDS(resultname)

range_df <- readRDS('rrbs_clean_data/predictedMeth_range_quality.RDS')
range_df$unit <- rownames(range_df)

library(plyr)
result <- join(result, range_df, by = 'unit', type = "left")
result$padj <- p.adjust(result$`SUA:Pr(>|t|)`, method = 'fdr')
result$seqnames <- as.character(result$seqnames)
result <- result[(result$seqnames %in% paste0('chr', 1:22)), ]
saveRDS(result, "rrbs_clean_data/result_range.RDS")


library(openxlsx)
sigresult <- result[result$seqnames %in% paste0('chr', 1:22), ]
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'end', 'SUA:Pr(>|t|)', 'SUA:Estimate', 'padj')]
sigresult <- sigresult[sigresult$'SUA:Pr(>|t|)' < 0.05, ]
sigresult <- sigresult[order(sigresult$'SUA:Pr(>|t|)'), ]
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'SUA:Pr(>|t|)', 'SUA:Estimate', 'padj')]
colnames(sigresult) <- c('id_row_number', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
saveRDS(sigresult, file = 'rrbs_clean_data/sigresult.RDS')
write.xlsx(sigresult, file = 'rrbs_clean_data/sigresult.xlsx')

p_cut_v <- c(0.001, 0.01, 0.05, 1)

library(openxlsx)
for(p_cut in p_cut_v) {
  sigfilename <- gsub('.RDS', paste0('_', p_cut, '.xlsx'), resultname)
  sigfile <- result[result$`SUA:Pr(>|t|)` < p_cut, ]
  write.xlsx(sigfile, sigfilename)

  bedfilename <- gsub('.RDS', paste0('_', p_cut, '.BED'), resultname)
  bedfile <- result[result$`SUA:Pr(>|t|)` < p_cut, c('seqnames', 'start', 'end', 'unit')]
  write.table(bedfile, bedfilename, sep = ' ', row.names = FALSE, col.names = FALSE, quote = FALSE)
}

dir.create('combined-pvalues/data')
bedfilename <- 'combined-pvalues/data/pvals.bed'
bedfile <- result[(result$seqnames %in% paste0('chr', 1:22)), c('seqnames', 'start', 'end', 'SUA:Pr(>|t|)')]
bedfile$end <- bedfile$start + 1
bedfile$seqnames <- as.character(bedfile$seqnames)
# bedfile$seqnames <- factor(bedfile$seqnames, paste0('chr', 1:22))
bedfile <- bedfile[order(bedfile$seqnames, bedfile$start), ]
colnames(bedfile) <- c('chrom', 'start', 'end', 'p')
bedfile$p <- format(bedfile$p, digits=16, scientific=F)
write.table(bedfile, bedfilename, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

library(CMplot)
cmplotdata <- result[result$seqnames %in% paste0('chr', 1:22), ]
cmplotdata <- cmplotdata[, c('unit', 'seqnames', 'start', 'SUA:Pr(>|t|)')]
colnames(cmplotdata) <- c('SNP', 'Chromosome', 'Position', 'p-value')
cmplotdata$SNP <- as.factor(cmplotdata$SNP)
cmplotdata$Chromosome <- as.character(cmplotdata$Chromosome)
cmplotdata$Chromosome <- gsub('chr', '', cmplotdata$Chromosome)
cmplotdata$Chromosome <- as.factor(as.numeric(cmplotdata$Chromosome))

CMplot(cmplotdata,plot.type="c",chr.labels=paste("chr",c(1:22),sep=""),r=0.4,cir.legend=TRUE,
       outward=TRUE,cir.legend.col="black",cir.chr.h=0.3,chr.den.col="black",file="pdf", file.output = TRUE,
       memo="",dpi=600,verbose=TRUE)
file.rename('Circular-Manhattan.p-value.pdf', 'rrbs_clean_data/Circular-Manhattan_p-value.pdf')
CMplot(cmplotdata,plot.type="q",threshold=1e-6,
       signal.pch=19,signal.cex=1,box=FALSE,multracks=
       FALSE,memo="",dpi=600,file = "pdf",file.output=TRUE,verbose=TRUE)
file.rename('QQplot.p-value.pdf', 'rrbs_clean_data/QQplot_p-value.pdf')
