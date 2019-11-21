setwd("/home/w/working/20190326_认知功能甲基化分析")
getwd()
options(stringsAsFactors = FALSE)

resultname <- 'rrbs_clean_data1/lmer.RDS'

result <- readRDS(resultname)
result$unit <- as.character(result$unit)
range <- readRDS('rrbs_clean_data1/predictedMeth_range_quality.RDS')

range_df <- as.data.frame(range)
range_df$unit <- paste0('r', 1:nrow(range_df))
range_df$unit <- as.character(range_df$unit)

library(plyr)
result <- join(result, range_df, by = 'unit')
result$padj <- p.adjust(result$`cognitive_function_score:Pr(>|t|)`, method = 'fdr')
result$seqnames <- as.character(result$seqnames)
saveRDS(result, "rrbs_clean_data1/result_range.RDS")

#there was an error while installing the package of CMplot, which may be due to the network. And this package was installed by weilong actually
#install.packages("CMplot")
library(CMplot)
cmplotdata <- result[result$seqnames %in% paste0('chr', 1:22), ]
cmplotdata <- cmplotdata[, c('unit', 'seqnames', 'start', 'cognitive_function_score:Pr(>|t|)')]
colnames(cmplotdata) <- c('SNP', 'Chromosome', 'Position', 'p-value')
cmplotdata$SNP <- as.factor(cmplotdata$SNP)
cmplotdata$Chromosome <- as.character(cmplotdata$Chromosome)
cmplotdata$Chromosome <- gsub('chr', '', cmplotdata$Chromosome)
cmplotdata$Chromosome <- as.factor(as.numeric(cmplotdata$Chromosome))
# CMplot(cmplotdata,plot.type="c",chr.labels=paste("chr",c(1:22,"X"),sep=""),r=0.4,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=0.3,chr.den.col="black",file="jpg",
#        memo="",dpi=600,file.output=TRUE,verbose=TRUE)
# CMplot(cmplotdata,plot.type="q",threshold=1e-6,
#        signal.pch=19,signal.cex=1,box=FALSE,multracks=
#          FALSE,file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE)


p_cut_v <- c(0.001, 0.01, 0.05)

for(p_cut in p_cut_v) {
  sigfilename <- gsub('.RDS', paste0('_', p_cut, '.csv'), resultname)
  sigfile <- result[result$`cognitive_function_score:Pr(>|t|)` < p_cut, ]

  bedfilename <- gsub('.RDS', paste0('_', p_cut, '.BED'), resultname)
  bedfile <- result[result$`cognitive_function_score:Pr(>|t|)` < p_cut, c('seqnames', 'start', 'end', 'unit')]
  write.table(bedfile, bedfilename, sep = ' ', row.names = FALSE, col.names = FALSE, quote = FALSE)
}


dir.create('combined-pvalues/data11')
bedfilename <- 'combined-pvalues/data11/pvals.bed'
bedfile <- result[(result$seqnames %in% paste0('chr', 1:22)), c('seqnames', 'start', 'end', 'cognitive_function_score:Pr(>|t|)')]
bedfile$end <- bedfile$start + 1
bedfile$seqnames <- as.character(bedfile$seqnames)
# bedfile$seqnames <- factor(bedfile$seqnames, paste0('chr', 1:22))
bedfile <- bedfile[order(bedfile$seqnames, bedfile$start), ]
colnames(bedfile) <- c('chrom', 'start', 'end', 'p')
bedfile$p <- format(bedfile$p, digits=16, scientific=F)
write.table(bedfile, bedfilename, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


result_range_filename <- gsub('.RDS', '_range.RDS', resultname)


library(CMplot)
library(ggplot2)

CMplot(cmplotdata,plot.type="c",chr.labels=paste("chr",c(1:22,"X"),sep=""),r=0.4,cir.legend=TRUE,
       outward=TRUE,cir.legend.col="black",cir.chr.h=0.3,chr.den.col="black",file="pdf", file.output = TRUE,
       memo="",dpi=600,verbose=TRUE)

file.rename('Circular-Manhattan.p-value.pdf', 'rrbs_clean_data1/Circular-Manhattan_p-value.pdf')
CMplot(cmplotdata,plot.type="q",threshold=1e-6,
       signal.pch=19,signal.cex=1,box=FALSE,multracks=
         FALSE,memo="",dpi=600,file = "pdf",file.output=TRUE,verbose=TRUE)
file.rename('QQplot.p-value.pdf', 'rrbs_clean_data1/QQplot_p-value.pdf')
# ggplotdata <- result[result$seqnames %in% paste0('chr', 1:22), ]
# ggplotdata <- ggplotdata[, c('unit', 'seqnames', 'start', 'cognitive_function_score:Pr(>|t|)', 'cognitive_function_score:Estimate')]
# colnames(ggplotdata) <- c('SNP', 'Chromosome', 'Position', 'pvalue', 'effectsize')
# pdf(file = "volcano-plot.pdf", width = 8, height = 8)
# # png(filename = "volcano-plot.png", width = 2500, height = 2500, res = 600)
# ggplot(data=ggplotdata, mapping=aes(x=effectsize, y=-log10(pvalue))) +
#   geom_point(size = 0.5) +
#   theme(legend.position = "none") +
#   theme_bw() +
#   xlab("Coefficient") + ylab("-log10 p-value")
# dev.off()



sigresult <- result[result$seqnames %in% paste0('chr', 1:22), ]
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'end', 'cognitive_function_score:Pr(>|t|)', 'cognitive_function_score:Estimate', 'padj')]
sigresult <- sigresult[sigresult$'cognitive_function_score:Pr(>|t|)' < 0.05, ]
sigresult <- sigresult[order(sigresult$'cognitive_function_score:Pr(>|t|)'), ]
saveRDS(sigresult, file = 'rrbs_clean_data1/sigresult.RDS')
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'cognitive_function_score:Pr(>|t|)', 'cognitive_function_score:Estimate', 'padj')]
colnames(sigresult) <- c('id_row_number', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
write.csv(sigresult, file = "rrbs_clean_data1/sigresult.csv", row.names = FALSE)

