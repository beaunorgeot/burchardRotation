#local
#setwd("/Users/beaunorgeot/burchardRotation/80_20/pr_pvalue")
#server
setwd("/media/BurchardRaid01/LabShare/Home/bnorgeot/Desktop/80_20/pr_by_pvalue")
library(caret)
library(dplyr)
library(doMC)
registerDoMC(cores = 54)
#local
#genotypes = read.table("../PR_top200_MAFp01.ALL.ind-PR-ALL.replaceNA.rIN", sep = " ",header = T,comment.char = "") 
#server
genotypes = read.table("/media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/asso/plink_80pct/top_200/plink/PR_top200_MAFp01.ALL.ind-PR-ALL.replaceNA.rIN", sep = " ",header = T,comment.char = "") 
genotypes = as.data.frame(lapply(genotypes, as.character), stringsAsFactors = T) %>% dplyr::mutate(personID = paste("X0",IID,sep="_")) %>% select(-IID)
pvalues = read.table("../Result_pur_1482individuals_80prct_Clean.MAFp01_Top1000.top200.txt",header = T)
load("../clinical_predFEV_train_pr.RData")
load("../clinical_predFEV_test_pr.RData")

genotypes_train = genotypes %>% filter(personID %in% clinical_predFEV_train_pr$personID)  
genotypes_train = inner_join(genotypes_train,clinical_predFEV_train_pr) %>% select(-c(age,Sex,BMI))

genotypes_test = genotypes %>% filter(personID %in% clinical_predFEV_test_pr$personID)
genotypes_test = inner_join(genotypes_test,clinical_predFEV_test_pr) %>% select(-c(age,Sex,BMI))

myControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
set.seed(666)

for (i in seq(10,30,10)){
  print(paste("starting forest",i/10,sep = " "))
  pvalues$SNP = as.character(pvalues$SNP)
  snps = pvalues %>% arrange(P)
  snps = snps[1:i,]
  train_response = genotypes_train$bdrGp
  nam = genotypes_train[,snps$SNP]
  nam$Pre.FEV1.perc.pred = genotypes_train$Pre.FEV1.perc.pred
  truth = genotypes_test$bdrGp
  test = genotypes_test[,snps$SNP]
  test$Pre.FEV1.perc.pred = genotypes_test$Pre.FEV1.perc.pred
  forest = paste('pvalue_rf',i,sep=".")
  assign(forest,train(nam,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T, trControl = myControl))
  pred = paste('pred_rf',i,sep=".")
  assign(pred,predict(get(forest),test))
  confus = paste('confus_rf',i, sep = ".")
  assign(confus,confusionMatrix(get(pred),truth))
  save(list=c(forest),file = paste("pvalue_rf",i,"RData", sep = "."))
  save(list=c(pred), file = paste("pred_rf",i,"RData",sep="."))
  save(list=c(confus), file = paste("confus_rf",i,"RData",sep="."))
}


############### inSample plotting ##########
resamps <- resamples(list(V_10 = pvalue_rf.10,V_20 = pvalue_rf.20), V_30 = pvalue_rf.30,V_40 = pvalue_rf.40,V_50 = pvalue_rf.50,
                     V_60 = pvalue_rf.60,V_70 = pvalue_rf.70,V_80 = pvalue_rf.80,V_90 = pvalue_rf.90,V_100 = pvalue_rf.100,
                     V_110 = pvalue_rf.110,V_120 = pvalue_rf.120,V_130 = pvalue_rf.130,V_140 = pvalue_rf.140,V_150 = pvalue_rf.150,
                     V_160 = pvalue_rf.160,V_170 = pvalue_rf.170,V_180 = pvalue_rf.180,V_190 = pvalue_rf.190,V_200 = pvalue_rf.200)
bwplot(resamps, layout = c(2, 1))
trellis.par.set(caretTheme())
dotplot(resamps, metric = "Accuracy")

#can get the mean oob error (mean of all trees)
mean(pvalue_rf.10$finalModel$err.rate[,1])
# get the confusion matrix for the final model. This provides misclaffication rates
pvalue_rf.10$finalModel$confusion
# see the misclassification rates for each class in that model
pvalue_rf.10$finalModel$confusion[,3]

###### out of sample plotting ##################
confuse_acc = paste("conf_acc",i,sep = "_")
confuse_kapp = paste("conf_kap",i,sep = "_")
assign(confuse_acc,confusionMatrix(preds,genotypes_test$bdrGp)$overall[1]) # the first 2 items in $overall are Accuracy(1) and Kappa(2)
assign(confuse_kapp,confusionMatrix(preds,genotypes_test$bdrGp)$overall[2])
i = i/10
confuse_acc[[i]] = confuse_acc
confuse_kapp[[i]] = confuse_kapp