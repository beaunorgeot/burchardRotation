
library(caret)
library(dplyr)
library(doMC) 
registerDoMC(cores=40)

setwd("/media/BurchardRaid01/LabShare/Home/bnorgeot/Desktop/80_20/models/")
myControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, verboseIter = T, returnResamp = "final")
#train(am ~.,data = mtcars,method = "rf", tuneLength =5,proximity = T, importance = T, type = 2)

######### aa ################
load("aa_pheno.RData")
load("aa_pop_pheno.RData")
load("aa_pop_meta_pheno.RData")
set.seed(1)
aa_pheno_rf = train(bdrGp ~. , data = aa_pheno, method = "rf", metric = "Kappa",tuneLength = 5, proximity = T, importance = T)
save(aa_pheno_rf, file = "aa_pheno_rf.RData")
aa_pop_pheno_rf = train(bdrGp ~. , data = aa_pop_pheno, method = "rf", metric = "Kappa",tuneLength = 5, proximity = T, importance = T)
save(aa_pop_pheno_rf, file = "aa_pop_pheno_rf.RData")
aa_pop_meta_pheno_rf = train(bdrGp ~. , data = aa_pop_meta_pheno, method = "rf", metric = "Kappa",tuneLength = 5, proximity = T, importance = T)
save(aa_pop_meta_pheno_rf, file = "aa_pop_meta_pheno_rf.RData")
##############################

######### pr ################
load("pr_pheno.RData")
load("pr_pop_pheno.RData")
load("pr_pop_meta_pheno.RData")
set.seed(2)
pr_pheno_rf = train(bdrGp ~. , data = pr_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(pr_pheno_rf, file = "pr_pheno_rf")
pr_pop_pheno_rf = train(bdrGp ~. , data = pr_pop_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(pr_pop_pheno_rf, file = "pr_pop_pheno_rf")
pr_pop_meta_pheno_rf = train(bdrGp ~. , data = pr_pop_meta_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(pr_pop_meta_pheno_rf, file = "pr_pop_meta_pheno_rf")
##############################

######### mex ################
load("mex_pheno.RData")
load("mex_pop_pheno.RData")
load("mex_pop_meta_pheno.RData")
set.seed(3)
mex_pheno_rf = train(bdrGp ~. , data = mex_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(mex_pheno_rf, file = "mex_pheno_rf")
mex_pop_pheno_rf = train(bdrGp ~. , data = mex_pop_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(mex_pop_pheno_rf, file = "mex_pop_pheno_rf")
mex_pop_meta_pheno_rf = train(bdrGp ~. , data = mex_pop_meta_pheno, method = "rf", metric = "Kappa", tuneLength = 5, proximity = T, importance = T)
save(mex_pop_meta_pheno_rf, file = "mex_pop_meta_pheno_rf")
##############################


###### feature importance ##########
#train the top model w/10 - 20 different random seeds. Everytime a var shows up in the topN increment it's count. Take the topN by count
# do the below extraction into a df, for each random seed

# extract a list of the vars by importance
ImpMeasure<-data.frame(varImp(bob)$importance)
ImpMeasure$Vars<-row.names(ImpMeasure)
#take the top n-variables
n = 20
bob1_vars = ImpMeasure[order(ImpMeasure$Overall,decreasing = T),][1:n,]