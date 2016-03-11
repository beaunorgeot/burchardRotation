#local
setwd("/Users/beaunorgeot/burchardRotation/80_20/pr_pvalue")
#server
#setwd("/media/BurchardRaid01/LabShare/Home/bnorgeot/Desktop/80_20/pr_by_pvalue")
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 54)
#local
genotypes = read.table("../PR_top200_MAFp01.ALL.ind-PR-ALL.replaceNA.rIN", sep = " ",header = T,comment.char = "") 
#server
#genotypes = read.table("/media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/asso/plink_80pct/top_200/plink/PR_top200_MAFp01.ALL.ind-PR-ALL.replaceNA.rIN", sep = " ",header = T,comment.char = "") 
genotypes = as.data.frame(lapply(genotypes, as.character), stringsAsFactors = T) %>% dplyr::mutate(personID = paste("X0",IID,sep="_")) %>% select(-IID)
pvalues = read.table("../Result_pur_1482individuals_80prct_Clean.MAFp01_Top1000.top200.txt",header = T)
load("../clinical_predFEV_train_pr.RData")
load("../clinical_predFEV_test_pr.RData")

genotypes_train = genotypes %>% filter(personID %in% clinical_predFEV_train_pr$personID)  
genotypes_train = inner_join(genotypes_train,clinical_predFEV_train_pr) %>% select(-c(age,Sex,BMI))
genotypes_test = genotypes %>% filter(personID %in% clinical_predFEV_test_pr$personID)
genotypes_test = inner_join(genotypes_test,clinical_predFEV_test_pr) %>% select(-c(age,Sex,BMI))
## 
myControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
##
i = 10
pvalues$SNP = as.character(pvalues$SNP)
snps = pvalues %>% arrange(P)
snps = snps[1:i,]
glmnet_raw_predictors = genotypes_train[,snps$SNP] 
glmnet_response = genotypes_train$bdrGp 
dummies = dummyVars("~.", data=glmnet_raw_predictors)
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
glmnet_raw_predictors$Pre.FEV1.perc.pred = genotypes_train$Pre.FEV1.perc.pred
truth = genotypes_test$bdrGp
test = genotypes_test[,snps$SNP]
test_dummies = dummyVars("~.", data=test)
test = as.data.frame(predict(test_dummies,newdata = test))
test$Pre.FEV1.perc.pred = genotypes_test$Pre.FEV1.perc.pred
enet = paste('pvalue_enet',i,sep=".")
assign(enet,train(glmnet_raw_predictors,glmnet_response, method = "glmnet", family = "binomial",trControl = myControl))



#make predictions
glmnet_clinical_predFEV_predictions_pr = predict(glmnet_clinical_predFEV_pr,testy)
glmnet_clinical_predFEV_confus_pr = confusionMatrix(glmnet_clinical_predFEV_predictions_pr,clinical_predFEV_test_pr$bdrGp) 