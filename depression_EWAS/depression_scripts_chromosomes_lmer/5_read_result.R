# setwd("/media/data/project/HUA_methylation_analysis_lmer/")
getwd()

options(stringsAsFactors = FALSE)

suffix <- Sys.getenv("suffix")
outdir <- paste0("rrbs_clean_data_",suffix,"/")

list_files <- list.files(outdir,"lmer")
chrs <- gsub("lmer-","",gsub(".RDS","",list_files))
results <- range_df <- list()
i <- 1
for (chr in chrs)
{
  resultname <- paste0(outdir,"lmer-",chr,".RDS")
  rangename <- paste0(outdir,"predictedMeth-",chr,"_range_quality.RDS")
  results[[i]] <- readRDS(resultname)
  range_df[[i]] <- readRDS(rangename)
  range_df[[i]]$unit <- rownames(range_df[[i]])
  i <- i+1
}
library(plyr)
result <- join(results[[1]], range_df[[1]], by = 'unit', type = "left")
i <- 1
for (chr in chrs[-1])
{
  result_range_df <- join(results[[i]], range_df[[i]], by = 'unit', type = "left")
  result <- rbind(result,result_range_df)
  i <- i+1
}
result$padj <- p.adjust(result$`CpG:Pr(>|t|)`, method = 'fdr')
result$seqnames <- as.character(result$seqnames)
result <- result[(result$seqnames %in% paste0('chr', 1:22)), ]
saveRDS(result, paste0(outdir,"result_range.RDS"))

library(openxlsx)
sigresult <- result[result$seqnames %in% paste0('chr', 1:22), ]
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'end', 'CpG:Pr(>|t|)', 'CpG:Estimate', 'padj')]
sigresult <- sigresult[sigresult$'CpG:Pr(>|t|)' < 0.05, ]
sigresult <- sigresult[order(sigresult$'CpG:Pr(>|t|)'), ]
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'CpG:Pr(>|t|)', 'CpG:Estimate', 'padj')]
colnames(sigresult) <- c('id_row_number', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
saveRDS(sigresult, file = paste0(outdir,"sigresult.RDS"))
write.xlsx(sigresult, file = paste0(outdir,"sigresult.xlsx"))

p_cut_v <- c(0.001, 0.01, 0.05, 1)

library(openxlsx)
for(p_cut in p_cut_v) {
  sigfilename <- paste0("lmer-",chr,"_", p_cut, '.xlsx')
  sigfile <- result[result$`CpG:Pr(>|t|)` < p_cut, ]
  write.xlsx(sigfile, sigfilename)
  bedfilename <- paste0(outdir,"lmer-",chr, "_", p_cut, '.BED')
  bedfile <- result[result$`CpG:Pr(>|t|)` < p_cut, c('seqnames', 'start', 'end', 'unit')]
  write.table(bedfile, bedfilename, sep = ' ', row.names = FALSE, col.names = FALSE, quote = FALSE)
}

dir.create('combined-pvalues/data')
bedfilename <- 'combined-pvalues/data/pvals.bed'
bedfile <- result[(result$seqnames %in% paste0('chr', 1:22)), c('seqnames', 'start', 'end', 'CpG:Pr(>|t|)')]
bedfile$end <- bedfile$start + 1
bedfile$seqnames <- as.character(bedfile$seqnames)
# bedfile$seqnames <- factor(bedfile$seqnames, paste0('chr', 1:22))
bedfile <- bedfile[order(bedfile$seqnames, bedfile$start), ]
colnames(bedfile) <- c('chrom', 'start', 'end', 'p')
bedfile$p <- format(bedfile$p, digits=16, scientific=F)
write.table(bedfile, bedfilename, sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

library(CMplot)
cmplotdata <- result[result$seqnames %in% paste0('chr', 1:22), ]
cmplotdata <- cmplotdata[, c('unit', 'seqnames', 'start', 'CpG:Pr(>|t|)')]
colnames(cmplotdata) <- c('SNP', 'Chromosome', 'Position', 'p-value')
cmplotdata$SNP <- as.factor(cmplotdata$SNP)
cmplotdata$Chromosome <- as.character(cmplotdata$Chromosome)
cmplotdata$Chromosome <- gsub('chr', '', cmplotdata$Chromosome)
cmplotdata$Chromosome <- as.factor(as.numeric(cmplotdata$Chromosome))

CMplot(cmplotdata,plot.type="c",chr.labels=paste("chr",c(1:22),sep=""),r=0.4,cir.legend=TRUE,
       outward=TRUE,cir.legend.col="black",cir.chr.h=0.3,chr.den.col="black",file="pdf", file.output = TRUE,
       memo="",dpi=600,verbose=TRUE)
file.rename('Circular-Manhattan.p-value.pdf', 'rrbs_clean_data_lmer/Circular-Manhattan_p-value.pdf')
CMplot(cmplotdata,plot.type="q",threshold=1e-6,
       signal.pch=19,signal.cex=1,box=FALSE,multracks=
         FALSE,memo="",dpi=600,file = "pdf",file.output=TRUE,verbose=TRUE)
file.rename('QQplot.p-value.pdf', paste0(outdir,"QQplot_p-value.pdf"))
