setwd("/media/data/project/HUA_methylation_analysis_matie/")
getwd()

library(plyr)
library(doMC)
doMC::registerDoMC(cores = 14)

phenfile <- 'SUA_sample_information.xlsx'
datafile <- 'rrbs_clean_data_matie/predictedMeth_m.RDS'
refactor_obj_name <- 'rrbs_clean_data_matie/refactor_obj.RDS'

library(openxlsx)
phen <- read.xlsx(phenfile)
data_d <- readRDS(datafile)
refactor_obj <- readRDS(refactor_obj_name)

data_d <- as.data.frame(data_d)
rownames(data_d) <- paste0("r", rownames(data_d))

rownames(phen) <- phen$no
phen <- phen[colnames(data_d), ]

SUA <- phen$SUA
age <- phen$age
gender <- as.factor(phen$gender)
GLU <- phen$GLU
CHOL <- phen$CHOL
TG <- phen$TG
HDLC <- phen$HDLC
BMI <- phen$BMI
scoretotal <- phen$scoretotal
fid <- as.factor(phen$family)
t1r <- which(phen$SUA_level == "H")
t2r <- which(phen$SUA_level == "L")
if(!all(phen$family.ID[t1r] == phen$family.ID[t2r])) {
    stop("family ids not maching")
}
PCs <- as.data.frame(refactor_obj$standard_pca)



require(gdata)
require(matie)

uobs <- as.numeric(data_d[1, ])
testdf <- residuals(lm(uobs ~ age + gender + scoretotal + BMI + systolic + GLU + CHOL + TG + HDLC + LDLC + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5))
testdf <- data.frame(meth = testdf[t1r] - testdf[t2r], phen = SUA[t1r] - SUA[t2r])
rawma <- ma(d = testdf)
pma <- ma.test(testdf, rawma)
uobssd <- sd(uobs)
uobssm <- mean(uobs)
uobscv <- uobssd/abs(uobssm)
dummy <- data.frame(unit = NA, pval = pma, mean = uobssm, sd = uobssd, cv = uobscv, rawma)

require(plyr)
require(doMC)
doMC::registerDoMC(cores = 14)
require(gdata)
require(data.table)
result <- rbindlist(alply(data_d, 1, function(obs) {
	tryCatch({
        uobs <- as.numeric(obs)
        testdf <- residuals(lm(uobs ~ age + gender + scoretotal + BMI + GLU + CHOL + TG + HDLC + PCs$PC1 + PCs$PC2 + PCs$PC3 + PCs$PC4 + PCs$PC5))
        testdf <- data.frame(meth = testdf[t1r] - testdf[t2r], phen = SUA[t1r] - SUA[t2r])
        rawma <- ma(d = testdf)
        pma <- ma.test(testdf, rawma)
        uobssd <- sd(uobs)
        uobssm <- mean(uobs)
        uobscv <- uobssd/abs(uobssm)
        sumt <- data.frame(unit = rownames(obs), pval = pma, mean = uobssm, sd = uobssd, cv = uobscv, rawma)
		return(sumt)
	}, error = function(e) {
		dummy$unit <- rownames(obs)
		return(dummy)
	})
}, .progress = "none", .parallel = TRUE))

result <- as.data.frame(result)

saveRDS(result, 'rrbs_clean_data_matie/matie.RDS')

