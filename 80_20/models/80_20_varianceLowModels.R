
# reducing variance
setwd("/Users/beaunorgeot/burchardRotation/80_20")
########### big guy ##################
#meta hits
load("80_20.meta.clean.RData") #metaClean, this is filtered for correct people, 1188
#clinical vars
clinical = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T)
clinical = clinical %>% mutate(personID = paste("X0", ID, sep = "_"), ethnicity = ethnicity.am) %>% select(personID, age, sex, bmicat.fix,ethnicity,bdrGp) %>% filter(personID %in% metaClean$personID)
save(clinical, file = "clinical.RData")
meta_clinical = inner_join(clinical,metaClean)

#local ancestry
localHits = read.table("ancestry_at_hits_collapsed.txt", sep = "\t") # this data from Zach:szpiech@gmail.com
localHits = cbind(personID = rownames(localHits),localHits)
rownames(localHits) = NULL
# convert all of these to factors
localHits = as.data.frame(lapply(localHits, factor))
localHits$personID = as.character(localHits$personID)
localHits = localHits %>% mutate(personID = paste("X0", personID, sep = "_")) 
localHitsTMP = localHits %>% filter(personID %in% clinical$personID)

meta_clinical_local = inner_join(meta_clinical,localHits)
save(meta_clinical_local, file = "meta_clinical_local.RData")
### run this ###########################################################################
library(caret)
library(dplyr)
library(doMC) 
registerDoMC(cores=40)

setwd("/media/BurchardRaid01/LabShare/Home/bnorgeot/Desktop/80_20/models/")
myControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = T, returnResamp = "final")
load("meta_clinical_local")

meta_clinical_local_train = meta_clinical_local %>% select(-personID)

set.seed(1)
meta_clinical_local_500_rf_1 = train(bdrGp ~. , data = aa_pheno_train, method = "rf", metric = "Kappa",tuneLength = 5, ntree = 500, proximity = T, importance = T)
set.seed(2)
meta_clinical_local_500_rf_2 = train(bdrGp ~. , data = aa_pheno_train, method = "rf", metric = "Kappa",tuneLength = 5, ntree = 500, proximity = T, importance = T)
set.seed(3)
meta_clinical_local_500_rf_3 = train(bdrGp ~. , data = aa_pheno_train, method = "rf", metric = "Kappa",tuneLength = 5, ntree = 500, proximity = T, importance = T)


##### small guys ###########


