setwd("/Users/beaunorgeot/burchardRotation/80_20")
genotypes = read.table("PR_top200_MAFp01.ALL.ind-PR.replaceNA.rIN", sep = " ",header = T,comment.char = "") 
genotypes = as.data.frame(lapply(genotypes, as.character), stringsAsFactors = T) %>% dplyr::mutate(personID = paste("X0",IID,sep="_")) %>% select(-IID)
pvalues = read.table("Result_pur_1482individuals_80prct_Clean.MAFp01_Top1000.top200.txt",header = T)
load("clinical_predFEV_train_pr.RData")
load("clinical_predFEV_test_pr.RData")
load("clinical_predFEV_train_pr.RData")

library(caret)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(cores = 54)

genotypes_train = genotypes %>% filter(personID %in% clinical_predFEV_train_pr$personID)  
genotypes_train = inner_join(genotypes_train,clinical_predFEV_train_pr) %>% select(-c(age,Sex,BMI))

myControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
model_list = list()
set.seed(666)

i = 10
pvalues$SNP = as.character(pvalues$SNP)
snps = pvalues %>% arrange(P)
snps = snps[1:i,]
train_response = genotypes_train$bdrGp
nam = genotypes_train[,snps$SNP]
nam$Pre.FEV1.perc.pred = genotypes_train$Pre.FEV1.perc.pred
forest = paste('pvalue_rf',i,sep=".")
assign(forest,train(nam,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T, trControl = myControl))
save(list=c(forest),file = paste("pvalue_rf",i,"RData", sep = "."))  

########### testing ##################
genotypes_test = genotypes %>% filter(personID %in% clinical_predFEV_test_pr$personID)
genotypes_test = inner_join(genotypes_test,clinical_predFEV_test_pr) %>% select(-c(age,Sex,BMI))
