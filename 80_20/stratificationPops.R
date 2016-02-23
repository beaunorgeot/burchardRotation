
setwd("/Users/beaunorgeot/burchardRotation/80_20")
mex_80_all = allFeaturesResponse %>% filter(ethnicity == "Mexican")
set.seed(42)
inTrain <- createDataPartition(mex_80_all$bdrGp, p = 0.8, list = F)
mex_80_trainSet <- mex_80_all[inTrain,]
mex_80_testSet <- mex_80_all[-inTrain,]

fileConn<-file("mex_80_trainSet_id.txt")
writeLines(as.character(mex_80_trainSet$personID), fileConn)
close(fileConn)
##########################################
aa_80_all = allFeaturesResponse %>% filter(ethnicity == "African American")
set.seed(43)
inTrain <- createDataPartition(aa_80_all$bdrGp, p = 0.8, list = F)
aa_80_trainSet <- aa_80_all[inTrain,]
aa_80_testSet <- aa_80_all[-inTrain,]

fileConn<-file("aa_80_trainSet_id.txt")
writeLines(as.character(aa_80_trainSet$personID), fileConn)
close(fileConn)
####################################
pr_80_all = allFeaturesResponse %>% filter(ethnicity == "Puerto Rican")
set.seed(44)
inTrain <- createDataPartition(pr_80_all$bdrGp, p = 0.8, list = F)
pr_80_trainSet <- pr_80_all[inTrain,]
pr_80_testSet <- pr_80_all[-inTrain,]

fileConn<-file("pr_80_trainSet_id.txt")
writeLines(as.character(pr_80_trainSet$personID), fileConn)
close(fileConn)