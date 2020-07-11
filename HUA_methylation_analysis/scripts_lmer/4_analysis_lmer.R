# setwd("/media/data/project/HUA_methylation_analysis_lmer/")
getwd()

options(stringsAsFactors = FALSE)

chr <- Sys.getenv("chr")
suffix <- Sys.getenv("suffix")
outdir <- paste0("rrbs_clean_data_",suffix,"/")

library(plyr)
library(doMC)
doMC::registerDoMC(cores = 14)

phenfile <- 'SUA_sample_information.xlsx'
datafile <- paste0(outdir,"predictedMeth-",chr,"_m.RDS")
refactor_obj_name <- paste0(outdir,"refactor-",chr,"_obj.RDS")

library(openxlsx)
phen <- read.xlsx(phenfile)
data_d <- readRDS(datafile)
data_d <- as.data.frame(data_d)

rownames(phen) <- phen$no
phen <- phen[colnames(data_d), ]

SUA <- phen$SUA
age <- phen$age
gender <- as.factor(phen$gender)
fid <- as.factor(phen$family)
PCs <- readRDS(refactor_obj_name)

require(gdata)
require(lmerTest)

dummy <- as.numeric(data_d[1, ])
dummy <- summary(lmer(dummy ~ SUA + age + gender + systolic + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + (1|fid), data=phen))$coefficient
dummy <- unmatrix(dummy, byrow = TRUE)
dummy <- cbind.data.frame(unit = rownames(data_d[1, ]), t(dummy))
dummy[1, ] <- NA

require(plyr)
require(doMC)
doMC::registerDoMC(cores = 14)
require(gdata)
require(data.table)
result <- rbindlist(alply(data_d, 1, function(obs) {
	tryCatch({
		sumt <- summary(lmer(as.numeric(obs) ~ SUA + age + gender + systolic + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + (1|fid), data=phen))$coefficient
		sumt <- unmatrix(sumt, byrow = TRUE)
		sumt <- cbind.data.frame(unit = rownames(obs), t(sumt))
		sumt$unit <- rownames(obs)
		return(sumt)
	}, error = function(e) {
		dummyresult <- dummy
		dummyresult$unit <- rownames(obs)
		return(dummyresult)
	})}, .progress = "none", .parallel = TRUE))

result <- as.data.frame(result)

saveRDS(result, paste0(outdir,"lmer-",chr,".RDS')
