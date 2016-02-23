
#FILTER TO TAKE ONLY THE 80% PEOPLE!!!!

#bash command for pruned snps:
#1 Take all SNPS from the 80% that are more significant than 10^-5
# take just the positions: cat top80aa16.txt | awk '{print $2}' |grep -ve SNP > top80aa16.id.txt
# 2 get the genotypes: for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do plink --bfile 
#/media/BurchardRaid01/LabShare/Data/wgs_shared/plink/joint.1484.2015_1030.chr$chr.fix2.pass.dedup.geno-p1.name 
#--extract top80aa16.id.txt --make-bed --out genoAA$chr;done

#3 prune for LD: for chr in {1..22}; do plink --bfile genoMeta$chr --indep 50 5 2  --out pruneMeta$chr; done
# 4. extract genotype for pruned files
#for chr in {1..22}; do plink --bfile /media/BurchardRaid01/LabShare/Data/wgs_shared/plink/joint.1484.2015_1030.chr$chr.fix2.pass.dedup.geno-p1.name 
#--extract pruneAA$chr.prune.in --recode vcf --out pruneAA$chr;done
# combine those: cat pruneAA*.vcf | grep "CHROM" | sort | uniq > pruneAACombined ; cat pruneAA*.vcf | grep -ve "^#" >> pruneAACombined

# MONDay
#PUR and AA are done. Just do Mex and then rsync the pruneMETA/AA/PUR/MEX/Combined files to me
# pruneMetaGenoCombined, pruneAACombined, prunePURCombined, pruneMEXCombined

#TODO: Make a model with just the top(genome wide sig snps). 
#NEXT: Get the genotypes for the 20% holdout set for each population, so can make OOB predictions

setwd("/Users/beaunorgeot/burchardRotation/80_20")
library(dplyr)

# phenotypes
pheno = read.table("../phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T)
pheno = pheno %>% mutate(personID = paste("X0", ID, sep = "_"), ethnicity = ethnicity.am) %>% select(personID, age, sex, bmicat.fix,ethnicity,bdrGp)

AA = read.table("pruneAACombined",sep = "\t", header = T,comment.char = "")
AA = AA %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
AA = t(AA)
colnames(AA) = AA[1,]
AA = AA[-1,]
AA = as.data.frame(AA)
AA = cbind(personID = rownames(AA),AA)
rownames(AA) = NULL
#filter for the training cases
aa80 = read.table("aa_80_trainSet_id.txt", sep = "\t", header = F)
aaClean = AA %>% filter(personID %in% aa80$V1) #dim = 392
save(aaClean, file = "80_20.AA.clean.RData")


PUR = read.table("prunePURCombined",sep = "\t", header = T,comment.char = "")
PUR = PUR %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
PUR = t(PUR)
colnames(PUR) = PUR[1,]
PUR = PUR[-1,]
PUR = as.data.frame(PUR)
PUR = cbind(personID = rownames(PUR),PUR)
rownames(PUR) = NULL
#filter for training data
pur80 = read.table("pr_80_trainSet_id.txt", sep = "\t", header = F)
purClean = MEX %>% filter(personID %in% pur80$V1)
save(purClean, file = "80_20.PUR.clean.RData") #dim = 396

MEX = read.table("pruneMEXCombined",sep = "\t", header = T,comment.char = "")
MEX = MEX %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
MEX = t(MEX)
colnames(MEX) = MEX[1,]
MEX = MEX[-1,]
MEX = as.data.frame(MEX)
MEX = cbind(personID = rownames(MEX),MEX)
rownames(MEX) = NULL
#take only training 80%
mex80 = read.table("mex_80_trainSet_id.txt", sep = "\t", header = F)
mexClean = MEX %>% filter(personID %in% mex80$V1)
save(mexClean, file = "80_20.MEX.clean.RData") # dim = 400


meta = read.table("pruneMetaGenoCombined",sep = "\t", header = T,comment.char = "")
meta = meta %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
meta = t(meta)
colnames(meta) = meta[1,]
meta = meta[-1,]
meta = as.data.frame(meta)
meta = cbind(personID = rownames(meta),meta)
rownames(meta) = NULL
#take only the 80% training set
allTraining =c(as.character(aaClean$personID),as.character(purClean$personID),as.character(mexClean$personID))
metaClean = meta %>% filter(personID %in% as.character(allTraining)) #dim(metaClean) #1188
save(meta, file = "80_20.meta.clean.RData")


########## local ancestry #############
allLocal = read.table("hits.all.subset.txt", sep = "\t",header = T)
localHits = read.table("ancestry_at_hits.txt", sep = "\t", header = T)
#######################

#### data for aa models ########
aa_pheno = pheno %>% filter(ethnicity == "African American", personID %in% aaClean$personID) #dim 392 x 6
#aa_meta = metaClean %>% filter(personID %in% aaClean$personID) #dim 392 x 228; meta not used on it's own
aa_pop_pheno = inner_join(aa_pheno, aaClean) #dim 392 x 16
#aa_meta_pheno = inner_join(aa_pheno,aa_meta) #dim: 392 x 233
aa_pop_meta_pheno = inner_join(aa_meta_pheno,aaClean) #dim: 392 x 243
save(aa_pheno, file = "aa_pheno.RData")
save(aa_pop_pheno, file = "aa_pop_pheno.RData")
save(aa_pop_meta_pheno, file = "aa_pop_meta_pheno.RData")
# CURRENTLY MISSING THE LOCAL ANCESTRY!!!!
################

#### data for mex models ########
mex_pheno = pheno %>% filter(ethnicity == "Mexican", personID %in% mexClean$personID) #dim 400 x 6
#mex_meta = metaClean %>% filter(personID %in% mexClean$personID) #dim 400 x 228; meta not used on it's own
mex_pop_pheno = inner_join(mex_pheno, mexClean) #dim 400 x 10
mex_pop_meta_pheno = inner_join(mex_pop_pheno,mex_meta) #dim: 400 x 237
save(mex_pheno, file = "mex_pheno.RData")
save(mex_pop_pheno, file = "mex_pop_pheno.RData")
save(mex_pop_meta_pheno, file = "mex_pop_meta_pheno.RData")
# CURRENTLY MISSING THE LOCAL ANCESTRY!!!!
################

#### data for pr models ########
pr_pheno = pheno %>% filter(ethnicity == "Puerto Rican", personID %in% purClean$personID) #dim 396 x 6
#pr_meta = metaClean %>% filter(personID %in% purClean$personID) #dim 396 x 228; meta not used on it's own
pr_pop_pheno = inner_join(pr_pheno, purClean) #dim 396 x 10
pr_pop_meta_pheno = inner_join(pr_pop_pheno,pr_meta) #dim: 396 x 237
save(pr_pheno, file = "pr_pheno.RData")
save(pr_pop_pheno, file = "pr_pop_pheno.RData")
save(pr_pop_meta_pheno, file = "pr_pop_meta_pheno.RData")
# CURRENTLY MISSING THE LOCAL ANCESTRY!!!!
################