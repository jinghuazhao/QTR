setwd("/media/data/project/HUA_methylation_analysis_matie/")
getwd()
range_result_name <- 'rrbs_clean_data_matie/result_range.RDS'

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt", version = "3.8")



result <- readRDS(range_result_name)
result$seqnames <- gsub('chr', '', result$seqnames)
result$seqnames <- as.numeric(result$seqnames)
result <- result[!is.na(result$seqnames), ]


library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
# attr: entrezgene, ensembl_gene_id
attr <- listAttributes(ensembl)
filter <- listFilters(ensembl)


id_type <- 'entrezgene_id'
anno <- getBM(attributes = c(id_type, 'chromosome_name', 'start_position', 'end_position'), mart = ensembl)


library(plyr)
result_anno <- adply(.data = result, .margins = 1, .fun = function(x) {
  anno[(x$seqnames == anno$chromosome_name) & (x$start >= anno$start_position) & (x$end <= anno$end_position), ]
}, .progress = "text")

result_anno_filename <- gsub('.RDS', paste0('_', id_type, '.RDS'), range_result_name)
saveRDS(result_anno, result_anno_filename)
result_anno_csv <- gsub('.RDS', paste0('_', id_type, '.csv'), range_result_name)
write.csv(result_anno, file = result_anno_csv, row.names = FALSE)




id_type <- 'ensembl_gene_id'
anno <- getBM(attributes = c(id_type, 'chromosome_name', 'start_position', 'end_position', 'description', 'hgnc_symbol'), mart = ensembl)


library(plyr)
result_anno <- adply(.data = result, .margins = 1, .fun = function(x) {
  rst <- anno[(x$seqnames == anno$chromosome_name) & (x$start >= anno$start_position) & (x$end <= anno$end_position), ]
  extra <- data.frame(up_bp = NA, up_ensembl_gene_id = NA, up_hgnc_symbol = NA, down_bp = NA, down_ensembl_gene_id = NA, down_hgnc_symbol = NA)
  if(nrow(rst) == 0) {
    rst[1, ] <- NA

    rst_up <- anno[(x$seqnames == anno$chromosome_name) & (x$start < anno$start_position) & (x$end <= anno$end_position), ]
    rst_up <- rst_up[order(rst_up$start_position), ]
    if(nrow(rst_up) > 0) {
      rst_up <- rst_up[1, ]
      extra$up_bp <- x$start - rst_up$start_position
      extra$up_ensembl_gene_id <- rst_up$ensembl_gene_id
      extra$up_hgnc_symbol <- rst_up$hgnc_symbol
    }

    rst_down <- anno[(x$seqnames == anno$chromosome_name) & (x$start >= anno$start_position) & (x$end > anno$end_position), ]
    rst_down <- rst_down[order(rst_down$end_position, decreasing = TRUE), ]
    if(nrow(rst_down) > 0) {
      rst_down <- rst_down[1, ]
      extra$down_bp <- x$end - rst_down$end_position
      extra$down_ensembl_gene_id <- rst_down$ensembl_gene_id
      extra$down_hgnc_symbol <- rst_down$hgnc_symbol
    }
  }
  rst <- cbind(rst, extra)
  return(rst)
}, .progress = "text")

result_anno_filename <- gsub('.RDS', paste0('_', id_type, '.RDS'), range_result_name)
saveRDS(result_anno, result_anno_filename)
result_anno_csv <- gsub('.RDS', paste0('_', id_type, '.csv'), range_result_name)
write.csv(result_anno, file = result_anno_csv, row.names = FALSE)


sigresult <- readRDS('rrbs_clean_data_matie/sigresult.RDS')
anno <- readRDS('rrbs_clean_data_matie/result_range_ensembl_gene_id.RDS')
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'SUA:Pr(>|t|)', 'SUA:Estimate', 'padj')]
colnames(sigresult) <- c('unit', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
anno <- anno[, c('unit', "ensembl_gene_id", "hgnc_symbol", "up_bp", "up_ensembl_gene_id", "up_hgnc_symbol","down_bp", "down_ensembl_gene_id", "down_hgnc_symbol" )]
library(plyr)
sigresult <- join(sigresult, anno)
write.csv(sigresult, file = "rrbs_clean_data_matie/sigresult_anno.csv", row.names = FALSE)
write.xlsx(sigresult, file = "rrbs_clean_data_matie/sigresult_anno.xlsx", row.names = FALSE)


sigresult <- readRDS('rrbs_clean_data_matie/sigresult.RDS')
#Note:change the "entrezgene" to "entrezgene_id"
anno <- readRDS('rrbs_clean_data_matie/result_range_entrezgene_id.RDS')
sigresult <- sigresult[, c('unit', 'seqnames', 'start', 'SUA:Pr(>|t|)', 'SUA:Estimate', 'padj')]
colnames(sigresult) <- c('unit', 'Chromosome', 'Position', 'pvalue', 'effectsize', 'padj')
anno <- anno[, c('unit', "entrezgene_id")]
library(plyr)
sigresult <- join(sigresult, anno)
write.csv(sigresult, file = "rrbs_clean_data_matie/sigresult_anno_entrezgene.csv", row.names = FALSE)
write.xlsx(sigresult, file = "rrbs_clean_data_matie/sigresult_anno_entrezgene.xlsx", row.names = FALSE)
