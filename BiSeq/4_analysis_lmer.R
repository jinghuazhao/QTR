setwd("/home/w/working/20190326_认知功能甲基化分析")
getwd()

library(plyr)
library(doMC)
doMC::registerDoMC(cores = 14)

phenfile <- '20190326_sample_information(copy).xlsx'
datafile <- 'rrbs_clean_data1/predictedMeth_m.RDS'
refactor_obj_name <- 'rrbs_clean_data1/refactor_obj.RDS'

library(openxlsx)
phen <- read.xlsx(phenfile)
data_d <- readRDS(datafile)
refactor_obj <- readRDS(refactor_obj_name)

data_d <- as.data.frame(data_d)
rownames(data_d) <- paste0("r", rownames(data_d))

rownames(phen) <- phen$no
phen <- phen[colnames(data_d), ]

cognitive_function_score <- phen$cognitive_function_score
age <- phen$age
gender <- as.factor(phen$gender)
fid <- as.factor(phen$family)
PCs <- as.data.frame(refactor_obj$standard_pca)
eversmoking <- as.factor(phen$eversmoking)
everdrinking <- as.factor(phen$everdrinking)
diastolic <- phen$diastolic
SUA <- phen$SUA


require(gdata)
require(lmerTest)

dummy <- as.numeric(data_d[1, ])
dummy <- summary(lmer(dummy ~ cognitive_function_score + age + gender + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + diastolic + (1|fid) ))$coefficient
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
		sumt <- summary(lmer(as.numeric(obs) ~ cognitive_function_score + age + gender + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5 + diastolic + (1|fid) ))$coefficient
		sumt <- unmatrix(sumt, byrow = TRUE)
		sumt <- cbind.data.frame(unit = rownames(obs), t(sumt))
		sumt$unit <- rownames(obs)
		return(sumt)
	}, error = function(e) {
		dummy$unit <- rownames(obs)
		return(dummy)
	})
}, .progress = "none", .parallel = TRUE))

result <- as.data.frame(result)

saveRDS(result, 'rrbs_clean_data1/lmer.RDS')

