
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
#1. remove bmi and Sex from the minimal clinical model to show effect of pre-fev alone
#2. run rfe() on these genos to see if there's any that can be saved
#3. plot the models, just fev, minimal clinical, local, global, geno side by side to show the effects
#4. Discuss the study design, why ethnicity is not predictive
#5. Make list of all RF and glmnet methods tried from the very beginning with their brief rational