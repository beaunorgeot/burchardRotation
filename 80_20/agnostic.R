library(caret)
library(dplyr)
library(tidyr)
library(doMC)
registerDoMC(cores = 54)

setwd("/media/BurchardRaid01/LabShare/Home/bnorgeot/Desktop/80_20/agnostic")


load("clinical_train.RData")
pr_response = clinical_train %>% filter(ethnicity == "Puerto Rican") %>% separate(col = personID, into = c("crap","IID"), sep = "_") %>% select(c(IID,bdrGp))


myControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3 )

set.seed(666)
model_list = list()
for (i in 1:3124){
  print(paste("starting forest",i,sep = " "))
  nam = read.table(paste('/media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/asso/plink_80pct/var_4000/plink/var.',i,".final.PR.replaceNA.rIN",sep=""), header=T,sep="")
  nam = inner_join(nam, pr_response)
  train_response = nam$bdrGp
  nam = nam %>% select(-c(IID,bdrGp))
  nam = as.data.frame(sapply(nam, as.factor))
  forest = paste('rf_1',i,sep="_")
  assign(forest,train(nam,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T, trControl = myControl))
  model_list[i] = forest
}
# extract a list of the vars by importance as a seperate df for each model
for (i in 1:length(model_list)){
  vars = paste('important_rf_vars',i,sep="_")
  assign(vars,data.frame(varImp(get(model_list[[i]]))$importance))
  print(class(get(vars)))
  tempTab = get(vars)
  tempTab[paste("myVars", i, sep="_")] <- row.names(tempTab)
  colnames(tempTab) = paste(colnames(tempTab),i,sep = "_")
  assign(vars,tempTab)
}

#bind the important vars into 1 df
all_vars = important_rf_vars_1
for (i in 2:length(model_list)){
  vars = paste('important_rf_vars',i,sep="_")
  all_vars = cbind(all_vars, get(vars))
}

# either extract the top n important vars, or remove the top_n() step to simply return a sorted list of vars with their relative importance
# YOU NEED TO CHANGE THE RESPONSE select(contain("your actual variable")),, in the first select from "auto" to whatever you want
# look at all_vars to choose a name
n = 4000
my_avgTop_vars = all_vars %>% select(Vars = myVars_1_1,contains("low")) %>% mutate(avg_imp = rowSums(.[, -1])/length(model_list)) %>% arrange(desc(avg_imp)) %>% select(Vars,avg_imp) %>% top_n(n = n,wt=avg_imp)

############
# Get a final ranking of the top variables, by combining the previous top variables from each run together

# /media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/asso/plink_80pct/var_4000/RData/merged_model_list_1-3124.RData
# name = df
#local location
#top_4000 = read.table("PR_feature4000_MAFp01.ALL.ind-PR.replaceNA.rIN",sep = " ",header = T,comment.char = "") %>% dplyr::mutate(personID = paste("X0",IID,sep="_")) %>% select(-IID)
top_4000 = read.table("/media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/asso/plink_80pct/feature4000/plink/PR_feature4000_MAFp01.ALL.ind-PR.replaceNA.rIN",sep = " ",header = T,comment.char = "") %>% dplyr::mutate(personID = paste("X0",IID,sep="_")) %>% select(-IID)
temp = top_4000$personID
top_4000 = top_4000 %>% select(-personID)
top_4000 = data.frame(lapply(top_4000,as.character), stringsAsFactors = T)
top_4000$personID = temp
load("clinical_predFEV_train_pr.RData")
load("clinical_predFEV_test_pr.RData")
# make a personID for the input nam = as.data.frame(sapply(nam, as.factor))
top_4000_train = top_4000 %>% filter(personID %in% clinical_predFEV_train_pr$personID)
top_4000_train = inner_join(top_4000_train,clinical_predFEV_train_pr) %>% select(-c(age,Sex,BMI,Pre.FEV1.perc.pred))

#top_4000_test = top_4000 %>% filter(personID %in% clinical_predFEV_test_pr$personID)
#top_4000_test  = inner_join(top_4000_test,clinical_predFEV_test_pr) %>% select(-c(age,Sex,BMI,Pre.FEV1.perc.pred))

train_response = top_4000_train$bdrGp
train_4000 = top_4000_train %>% select(-c(bdrGp,personID))

set.seed(42)
myControl = trainControl(method = "repeatedcv", number = 5, repeats = 5)
rf_4000 = train(train_4000,train_response,method = "rf", metric = "Kappa",tuneLength =5,ntree = 1500, importance = T, trControl = myControl)

best_4000 = varImp(rf_4000)$importance
best_4000 = cbind(rownames(best_4000),best_4000)
rownames(best_4000) = NULL
####################
# Build models using the top_n most important snps
set.seed(666)
model_list = list()
for(i in seq(10,200,10)){
  print(paste("starting forest",i/10,sep = " "))
  snps = best_4000 %>% top_n(n = i, wt = ) %>% select(SNP)
  train_response = train_response
  nam = as.data.frame(sapply(nam, as.factor))
  nam$Pre.FEV1.perc.pred = as.numeric(genotypes_train$Pre.FEV1.perc.pred)
  forest = paste('rf',i,sep="_")
  assign(forest,train(nam,train_response, method = "rf", metric = "Kappa", tuneLength = 3, ntree = 1500, importance = T, trControl = myControl))
  model_list[i] = forest
}
save(model_list,file = "model_list.RData")