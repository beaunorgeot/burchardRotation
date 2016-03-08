
library(caret)
library(dplyr)
# reducing variance
setwd("/Users/beaunorgeot/burchardRotation/80_20")
load("80_20.meta.clean.RData") #metaClean this is filtered for correct train people, 1188
people_80 = metaClean %>% select(personID)
load("all_meta_testing.RData") # to get the correct person IDs

################### build a model that should be perfectly predictive #############################
clinical = read.table("merge.wgs.phenotype.1484_2015_1209_all.txt",sep = "\t", header = T)
#which bmi values should be used? 
colnames(clinical)[grep("bmi",colnames(clinical))] # there's missing values in all, add the lite table which doen't have any
lite_c_train = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T) %>% mutate(personID = paste("X0", ID, sep = "_")) %>% filter(personID %in% people_80$personID)
lite_c_test = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T) %>% mutate(personID = paste("X0", ID, sep = "_")) %>% filter(personID %in% all_meta_testing$personID)

clinical_train = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>% filter(personID %in% people_80$personID) %>% mutate(BMI = factor(lite_c_train$bmicat.recode.collapse.fix)) %>% select(personID, Pre.FEV1.Meas,age, Sex, BMI,ethnicity,bdrGp)
clinical_test = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>%  filter(personID %in% all_meta_testing$personID) %>% mutate(BMI = factor(lite_c_test$bmicat.recode.collapse.fix)) %>% select(personID, Pre.FEV1.Meas,age, Sex, BMI,ethnicity,bdrGp)
save(clinical_train,file = "clinical_train.RData") #use clinical_train/test from now on!!! No more metaClean etc!
save(clinical_test,file = "clinical_test.RData")
# just train on these clin vars #######################
glmnet_raw_predictors = clinical_train %>% select(-c(personID,age,bdrGp))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
glmnet_clinical <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_test %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_predictions = predict(glmnet_clinical,testy)
glmnet_clinical_confus = confusionMatrix(glmnet_clinical_predictions,clinical_test$bdrGp) 
# Acc = .62, Kapp = .23
#Sensitivity : 0.6096          
#Specificity : 0.6216

#get the coeffecients of the final model
coef(glmnet_clinical$finalModel, glmnet_clinical$bestTune$.lambda) # this doesn't do what I wanted
varImp(glmnet_clinical) # returns the normal list of ranked vars

###########################################
# make a minimal clinical model
glmnet_raw_predictors = clinical_train %>% select(-c(personID,age,bdrGp,ethnicity))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(68)
glmnet_clinical_minimal <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_test %>% select(-c(bdrGp,age,personID,ethnicity))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_minimal_predictions = predict(glmnet_clinical_minimal,testy)
glmnet_clinical_minimal_confus = confusionMatrix(glmnet_clinical_minimal_predictions,clinical_test$bdrGp)
#Acc = 62, Kapp = .23

###########################################
# make a FEV model
glmnet_raw_predictors = clinical_train %>% select(-c(personID,age,bdrGp,ethnicity,Sex))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(68)
glmnet_fev <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_test %>% select(-c(bdrGp,age,personID,ethnicity,Sex))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_fev_predictions = predict(glmnet_fev,testy)
glmnet_fev_confus = confusionMatrix(glmnet_fev_predictions,clinical_test$bdrGp)
#Acc = 63, Kapp = .27


#####################
# do a forest just for poops
set.seed(999)
rf_clinical = train(dummy_glmnet_train, glmnet_response,method = "rf", metric = "Kappa",tuneLength =3,ntree = 1000, importance = T, type = 1)
# worse than glmnet
### ####################################
# munging for top 35 snps from each population by p-value(meta,aa,mex,pur)
#head -37 top80meta.txt | awk '{print $3}' | tail -n +2 > dSNPs.txt
#head -37 top80aa.txt | awk '{print $2}' | tail -n +3 >> dSNPs.txt
#head -37 top80mex.txt | awk '{print $2}' | tail -n +3 >> dSNPs.txt
#head -37 top80pur.txt | awk '{print $2}' | tail -n +3 >> dSNPs.txt

# 4. extract genotype
#for chr in {1..22}; do plink --bfile /media/BurchardRaid01/LabShare/Data/wgs_shared/plink/joint.1484.2015_1030.chr$chr.fix2.pass.dedup.geno-p1.name --extract dSNPs.txt --recode vcf --out genotypes_o_dSNPs$chr;done
#combine the results
#cat genotypes_o_dSNPs*.vcf | grep "CHROM" | sort | uniq > combine_genotypes_o_dSNPs ; cat genotypes_o_dSNPs*.vcf | grep -ve "^#" >> combine_genotypes_o_dSNPs.txt

combine_genotypes_o_dSNPs = read.table("combine_genotypes_o_dSNPs.txt", sep = "\t",header = T,comment.char = "")
combine_genotypes_o_dSNPs = combine_genotypes_o_dSNPs %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
combine_genotypes_o_dSNPs = t(combine_genotypes_o_dSNPs)
colnames(combine_genotypes_o_dSNPs) = combine_genotypes_o_dSNPs[1,]
combine_genotypes_o_dSNPs = combine_genotypes_o_dSNPs[-1,]
combine_genotypes_o_dSNPs = as.data.frame(combine_genotypes_o_dSNPs)
combine_genotypes_o_dSNPs = cbind(personID = rownames(combine_genotypes_o_dSNPs),combine_genotypes_o_dSNPs)
rownames(combine_genotypes_o_dSNPs) = NULL
#filter for the training cases

combine_genotypes_o_dSNPs_train = combine_genotypes_o_dSNPs %>% filter(personID %in% people_80$personID) #dim = 1188 x 141
combine_genotypes_o_dSNPs_test = combine_genotypes_o_dSNPs %>% filter(personID %in% all_meta_testing$personID) #dim = 294 x 141
geno_clin_train = inner_join(clinical_train,combine_genotypes_o_dSNPs_train)
geno_clin_test = inner_join(clinical_test,combine_genotypes_o_dSNPs_test)

# train them
glmnet_raw_predictors = geno_clin_train %>% select(-c(personID,age,bdrGp))
glmnet_response = geno_clin_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(42)
glmnet_geno_clin_mod <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
save(glmnet_geno_clin_mod, file = "glmnet_geno_clin_mod.RData") #this one takes a while to run on a laptop
#test it out of sample
testy = geno_clin_test %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_geno_clin_mod_predictions = predict(glmnet_geno_clin_mod,testy)
glmnet_geno_clin_mod_confus = confusionMatrix(glmnet_geno_clin_mod_predictions,geno_clin_test$bdrGp)
# Acc = .52, Kapp = .05

########################################
# just clinical and global ancestry
clinical_train_global = clinical_train
clinical_test_global = clinical_test
clinical_train_global$EUR = lite_c_train$EUR
clinical_test_global$EUR = lite_c_test$EUR
clinical_train_global$AFR = lite_c_train$AFR
clinical_test_global$AFR = lite_c_test$AFR
clinical_train_global$NAM = lite_c_train$NAM
clinical_test_global$NAM = lite_c_test$NAM

glmnet_raw_predictors = clinical_train_global %>% select(-c(personID,age,bdrGp,ethnicity))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(70)
glmnet_clinical_global <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_test_global %>% select(-c(bdrGp,age,personID,ethnicity))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_global_predictions = predict(glmnet_clinical_global,testy)
glmnet_clinical_minimal_confus = confusionMatrix(glmnet_clinical_global_predictions,clinical_test$bdrGp)
#Acc = .61, Kapp = .22

########################################
# just clinical and local ancestry
localHits = read.table("ancestry_at_hits_collapsed.txt", sep = "\t") # this data from Zach:szpiech@gmail.com
localHits = cbind(personID = rownames(localHits),localHits)
rownames(localHits) = NULL
# convert all of these to factors
localHits = as.data.frame(lapply(localHits, factor))
localHits$personID = as.character(localHits$personID)
localHits = localHits %>% mutate(personID = paste("X0", personID, sep = "_")) 
local_train = localHits %>% filter(personID %in% clinical_train$personID)
local_test = localHits %>% filter(personID %in% clinical_test$personID)

glmnet_raw_predictors = inner_join(clinical_train,local_train) %>% select(-c(personID,age,bdrGp))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(61)
glmnet_clinical_local <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = inner_join(clinical_test,local_test) %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_local_predictions = predict(glmnet_clinical_local,testy)
glmnet_clinical_confus = confusionMatrix(glmnet_clinical_local_predictions,clinical_test$bdrGp) 
# this one warrants some more examination. Acc = .6, Kapp = .2
# slightly increased sensitivity, at the expense of specifity
#Sensitivity : 0.6438          
#Specificity : 0.5608

#Next######################
#. 0 use % pre.fev b/c it supposedly corrects for stuff
clinical = read.table("merge.wgs.phenotype.1484_2015_1209_all.txt",sep = "\t", header = T)

lite_c_train = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T) %>% mutate(personID = paste("X0", ID, sep = "_")) %>% filter(personID %in% people_80$personID)
lite_c_test = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T) %>% mutate(personID = paste("X0", ID, sep = "_")) %>% filter(personID %in% all_meta_testing$personID)

clinical_predFEV_train = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>% filter(personID %in% people_80$personID) %>% mutate(BMI = factor(lite_c_train$bmicat.recode.collapse.fix)) %>% select(personID, Pre.FEV1.perc.pred,age, Sex, BMI,ethnicity,bdrGp)
clinical_predFEV_test = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>%  filter(personID %in% all_meta_testing$personID) %>% mutate(BMI = factor(lite_c_test$bmicat.recode.collapse.fix)) %>% select(personID, Pre.FEV1.perc.pred,age, Sex, BMI,ethnicity,bdrGp)

glmnet_raw_predictors = clinical_predFEV_train %>% select(-c(personID,age,bdrGp))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
glmnet_clinical_predFEV <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_predFEV_test %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_predFEV_predictions = predict(glmnet_clinical_predFEV,testy)
glmnet_clinical_predFEV_confus = confusionMatrix(glmnet_clinical_predFEV_predictions,clinical_test$bdrGp) 
#Accuracy : 0.6735, Kappa : 0.347
#Sensitivity : 0.6781          
#Specificity : 0.6689

###########################################
# repeat above for just PRs
clinical_predFEV_train_pr = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>% filter(personID %in% people_80$personID) %>% mutate(BMI = factor(lite_c_train$bmicat.recode.collapse.fix)) %>% filter (ethnicity == "Puerto Rican") %>% select(personID, Pre.FEV1.perc.pred,age, Sex, BMI,bdrGp)
clinical_predFEV_test_pr = clinical %>% mutate(personID = paste("X0", NWD_ID, sep = "_"), ethnicity = ethnicity.am) %>%  filter(personID %in% all_meta_testing$personID) %>% mutate(BMI = factor(lite_c_test$bmicat.recode.collapse.fix)) %>% filter (ethnicity == "Puerto Rican") %>% select(personID, Pre.FEV1.perc.pred,age, Sex, BMI,bdrGp)

lite_c_train_pr = lite_c_train %>% filter(ethnicity.am == "Puerto Rican")
lite_c_test_pr = lite_c_test %>% filter(ethnicity.am == "Puerto Rican")
clinical_predFEV_train_pr$EUR = lite_c_train_pr$EUR
clinical_predFEV_test_pr$EUR = lite_c_test_pr$EUR
clinical_predFEV_train_pr$AFR = lite_c_train_pr$AFR
clinical_predFEV_test_pr$AFR = lite_c_test_pr$AFR
clinical_predFEV_train_pr$NAM = lite_c_train_pr$NAM
clinical_predFEV_test_pr$NAM = lite_c_test_pr$NAM

glmnet_raw_predictors = clinical_predFEV_train_pr %>% select(-c(personID,age,bdrGp))
glmnet_response = clinical_predFEV_train_pr$bdrGp 
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
glmnet_clinical_predFEV_pr <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_predFEV_test_pr %>% select(-c(bdrGp,age,personID))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_predFEV_predictions_pr = predict(glmnet_clinical_predFEV_pr,testy)
glmnet_clinical_predFEV_confus_pr = confusionMatrix(glmnet_clinical_predFEV_predictions_pr,clinical_predFEV_test_pr$bdrGp) 
#Accuracy : 0.6907 || w/global same
#Kappa : 0.3819 || w/global .3822
#Sensitivity : 0.7292 || w/global .75
#Specificity : 0.6531 || w/global 0.6327
## using global ancestry provides a small increase in Sensitivity and Kappa, at the expense of a small decrease in Specifiity
## all things even, I would NOT use global

################################################# 
#1. remove pre-fev, leaving BMI and Sex from the minimal clinical model to show effect of pre-fev 
glmnet_raw_predictors = clinical_predFEV_train %>% select(-c(personID,age,bdrGp, Pre.FEV1.perc.pred))
glmnet_response = clinical_train$bdrGp
#cast into dummies
dummies = dummyVars("~.", data=glmnet_raw_predictors)
# turn dummies into new df
dummy_glmnet_train = as.data.frame(predict(dummies,newdata = glmnet_raw_predictors))
#train
myControl1 <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
set.seed(666)
glmnet_clinical_no_FEV <- train(dummy_glmnet_train, glmnet_response, method = "glmnet", family = "binomial", trControl = myControl1)
#test it out of sample
testy = clinical_predFEV_test %>% select(-c(bdrGp,age,personID,Pre.FEV1.perc.pred))
# make dummies dummy
test_dummies = dummyVars("~.", data=testy)
testy = as.data.frame(predict(test_dummies,newdata = testy))
#make predictions
glmnet_clinical_no_FEV_predictions = predict(glmnet_clinical_no_FEV,testy)
glmnet_clinical_no_FEV_confus = confusionMatrix(glmnet_clinical_no_FEV_predictions,clinical_test$bdrGp)
#Accuracy : 0.5204, Kappa : 0.04
# Sensitivity : 0.5068          
# Specificity : 0.5338
##############################################

#2. compare p-value model vs agnostic selection
# Get the top 100 SNPs for PRs by p-value. Or even just the top 10. Build a model to compare to top-N by agnostic