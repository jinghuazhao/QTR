setwd("/home/w/working/20190326_认知功能甲基化分析")
getwd()


options(stringsAsFactors = FALSE)

#the file "combined-pvalues/data/regions.sig.bed" was not generated while I practiced
result <- read.table('combined-pvalues/data/regions.sig.bed')

sigresult <- result[result$V7 < 0.1, ]
sigresult <- sigresult[order(sigresult$V7), ]
numberOfSignificance <- nrow(sigresult)


library(xtable)
sigresult_tex <- xtable(sigresult)
print(sigresult_tex, file = "rrbs_clean_data1/sigDMR.tex")
write.csv(sigresult, file = "rrbs_clean_data1/sigDMR.csv", row.names = FALSE)

resultname <- 'rrbs_clean_data1/lmer.RDS'

result <- readRDS(resultname)
result$unit <- as.character(result$unit)
range <- readRDS('rrbs_clean_data1/predictedMeth_range_quality.RDS')

range_df <- as.data.frame(range)
range_df$unit <- paste0('r', 1:nrow(range_df))

library(plyr)
result <- join(result, range_df, by = 'unit')

cpgresult <- result


library(ggplot2)

for(i in 1:nrow(sigresult)) {
  pdf(file = paste0('rrbs_clean_data1/DMR_', i, '.pdf'), width = 8.3, height = 4)
  plotdata1 <- sigresult[i, ]
  plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
  print(ggplot(data = plotdata2, mapping = aes(y = plotdata2$`cognitive_function_score:Estimate`, x = plotdata2$start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(paste0('DMR ', plotdata1$V1, ':', plotdata1$V2, '-', plotdata1$V3)) + theme_bw())
  dev.off()
}

library(grid)
library(gridExtra)
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) sprintf(paste0("%-", decimals, "s"), x)
}

pdf(file = paste0("rrbs_clean_data1/DMRs.pdf"), width = 10, height = 16.5)
num_per_page <- 8
col_per_page <- 2
row_per_page <- 4
num_page <- ceiling(nrow(sigresult)/num_per_page)
for(npage in 1:num_page) {
	startplot <- (npage -1) * num_per_page + 1
	if((startplot + num_per_page) > nrow(sigresult)) {
		endplot <- nrow(sigresult)
	} else {
		endplot <- startplot + num_per_page - 1
	}
	cat(paste0("start: ", startplot, ", end: ", endplot, "\n"))
	plot_list <- list()
	for(i in startplot:endplot) {
		plotdata1 <- sigresult[i, ]
		plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
		# plot_list[[i]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(paste0(LETTERS[i], ' DMR ', plotdata1$V1, ':', plotdata1$V2, '-', plotdata1$V3)) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(2)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
		plot_list[[i-startplot+1]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(bquote(bold(.(LETTERS[i-startplot+1])) * ' DMR ' * .(plotdata1$V1) * ':' * .(plotdata1$V2) * '-' * .(plotdata1$V3))) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(4)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
	}
	do.call('grid.arrange', c(plot_list, ncol = col_per_page, nrow = row_per_page))
}
dev.off()


library(grid)
library(gridExtra)
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) sprintf(paste0("%-", decimals, "s"), x)
}
pdf(file = paste0('rrbs_clean_data/DMRsr.pdf'), width = 15, height = 11.7)
num_per_page <- 9
col_per_page <- 3
row_per_page <- 3
num_page <- ceiling(nrow(sigresult)/num_per_page)
for(npage in 1:num_page) {
	startplot <- (npage -1) * num_per_page + 1
	if((startplot + num_per_page) > nrow(sigresult)) {
		endplot <- nrow(sigresult)
	} else {
		endplot <- startplot + num_per_page - 1
	}
	cat(paste0("start: ", startplot, ", end: ", endplot, "\n"))
	plot_list <- list()
	for(i in startplot:endplot) {
		plotdata1 <- sigresult[i, ]
		plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
		# plot_list[[i]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(paste0(LETTERS[i], ' DMR ', plotdata1$V1, ':', plotdata1$V2, '-', plotdata1$V3)) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(2)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
		plot_list[[i-startplot+1]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(bquote(bold(.(i)) * ' DMR ' * .(plotdata1$V1) * ':' * .(plotdata1$V2) * '-' * .(plotdata1$V3))) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(4)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
	}
	do.call('grid.arrange', c(plot_list, ncol = col_per_page, nrow = row_per_page))
}
dev.off()




library(ggpubr)
plot_list <- list()
fmt_dcimals <- function(decimals=0){
    # return a function responpsible for formatting the 
    # axis labels with a given number of decimals 
    function(x) sprintf(paste0("%-", decimals, "s"), x)
}
for(i in 1:12) {
  plotdata1 <- sigresult[i, ]
  plotdata2 <- cpgresult[(cpgresult$seqnames == plotdata1$V1) & (cpgresult$start >= plotdata1$V2) & (cpgresult$end <= plotdata1$V3), ]
  # plot_list[[i]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(paste0(LETTERS[i], ' DMR ', plotdata1$V1, ':', plotdata1$V2, '-', plotdata1$V3)) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(2)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  plot_list[[i]] <- ggplot(data = plotdata2, mapping = aes(y = `cognitive_function_score:Estimate`, x = start)) + geom_smooth() + geom_point() + xlab('BP') + ylab('Coefficient') + ggtitle(bquote('DMR ' * .(plotdata1$V1) * ':' * .(plotdata1$V2) * '-' * .(plotdata1$V3))) + theme_bw() + scale_y_continuous(labels = fmt_dcimals(4)) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
pdf(file = paste0('rrbs_clean_data1/DMRsr_ar.pdf'), width = 15, height = 11.7)
figure <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], ncol = 3, labels = LETTERS[1:12])
# annotate_figure(figure, top = text_grob("", vjust = 0.2, x = 0.08), fig.lab = "Figure 3", fig.lab.face = "bold")
figure
dev.off()

