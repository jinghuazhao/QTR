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

d <- data_d[1, ]
CpG <- as.numeric(d)
lmer.summary <- summary(lmer(SUA ~ CpG + age + gender + systolic + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + (1|fid),
                             data=phen))$coefficient
lmer.vector <- unmatrix(lmer.summary, byrow = TRUE)
dummy <- cbind.data.frame(unit = rownames(d), t(lmer.vector))
dummy[1, ] <- NA

require(plyr)
require(doMC)
doMC::registerDoMC(cores = 14)
require(gdata)
require(data.table)
result <- rbindlist(alply(data_d, 1, function(obs) {
	tryCatch({
                CpG <- as.numeric(obs)
		s <- summary(lmer(SUA~CpG + age + gender + systolic + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + (1|fid),
                                  data=phen))$coefficient
		us <- unmatrix(s, byrow = TRUE)
		sumt <- cbind.data.frame(unit = rownames(obs), t(us))
		sumt$unit <- rownames(obs)
		return(sumt)
	}, error = function(e) {
		dummyresult <- dummy
		dummyresult$unit <- rownames(obs)
		return(dummyresult)
	})}, .progress = "none", .parallel = TRUE))

result <- as.data.frame(result)

saveRDS(result, paste0(outdir,"lmer-",chr,".RDS"))
