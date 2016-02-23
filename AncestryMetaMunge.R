# get data from chromosome 5
setwd("/Users/beaunorgeot/burchardRotation")
library(dplyr)
meta = read.table("inputFromSkat/Result_BDR_Joint1484_Meta_3Studies_Hits_PC_Chr5.csv",header = T, sep = ",", skip=1)
View(meta)
str(meta)
metaMeta = meta %>% select(BP,CHR,SNP,A1,P,OR) %>% arrange(P)
#library(tidyr)
#create an easy table to paste results into ensembl
#split the SNP column into all necessary features
metaSplit = separate(metaMeta, into = c("chr","bp","ref","snp") ,col = SNP, sep= ":", remove = T)
ensFormat = data.frame(metaSplit$chr, metaSplit$bp, metaSplit$bp,paste(metaSplit$ref,metaSplit$snp,sep = "/"))
names(ensFormat) = c("chr","bpStart","bpEnd","ref/snp")
ensReady = paste(ensFormat$chr,ensFormat$bpStart,ensFormat$bpEnd,ensFormat$`ref/snp`, sep = " ")
# write to ens Ready format file
fileConn<-file("ensemblReady/ensemblMetaSnps.txt")
writeLines(ensReady, fileConn)
close(fileConn)

# recover RS ID numbers from ensembl Results to use with regulomedb.org
reg = read.table("ensemblResults/topMetaEnsemblTranscripts", sep = "\t", header = T,comment.char = "") #use comment.char="" so that read.table will not ignore lines starting with "#"
names(reg)
regReady = reg %>% select(Existing_variation)
#write regulom ready format to file
fileConn<-file("regulomReady/ensemblMetaSnpsRSIDs.txt")
writeLines(as.character(regReady$Existing_variation), fileConn)
close(fileConn)
#########################

# Get all meta data hits with a p-value better than e-05

sigMeta = read.csv("inputFromSkat/Result_BDR_Joint1484_Meta_3Studies_P_Lower_1e-5.csv/Sheet1-Table 1.csv", header = T)
#apparently there's no data in the other sheets
#sheet2 = read.csv("Sheet2-Table 1.csv", header = T)
#sheet3 = read.csv("Sheet3-Table 1.csv", header = T)
#sigMeta = rbind(sheet1,sheet2,sheet3)
sigMetaSplit = separate(sigMeta, into = c("chr","bp","ref","snp") ,col = SNP, sep= ":", remove = T)
SigensFormat = data.frame(sigMetaSplit$chr, sigMetaSplit$bp, sigMetaSplit$bp,paste(sigMetaSplit$ref,sigMetaSplit$snp,sep = "/"))
names(SigensFormat) = c("chr","bpStart","bpEnd","ref/snp")

#how many chr have significant snps?
metaChrTable = table(SigensFormat$chr) #1  10  11  12  14  15  16  18   2  22   3   4   5   6   8   9


SigensReady = paste(SigensFormat$chr,SigensFormat$bpStart,SigensFormat$bpEnd,SigensFormat$`ref/snp`, sep = " ")
# write ensembl Ready format to file
fileConn<-file("ensemblReady/sigEnsemblMetaSnps.txt")
writeLines(SigensReady, fileConn)
close(fileConn)
#for Dongleis
fileConn<-file("ensemblReady/dongSNPs")
writeLines(as.character(sigMeta$SNP), fileConn)
close(fileConn)

#get the ensemble results
sigReg = read.table("ensemblResults/sigMetaEnsemblTranscripts", sep = "\t", header = T,comment.char = "") #use comment.char="" so that read.table will not ignore lines starting with "#"
#sigReg has 337 observations, but sigMetaEnsemblTranscripts only has 205 (which is the correct number)
names(sigReg)
check = sigReg %>% group_by(Location,Existing_variation) %>% summarise(n = n())
View(check); dim(check) # some locations have multiple entries, sometimes as many as 18

dongleiSNPs = sigReg$Location
#write regulom Ready format (just rsID numbers) to file
fileConn<-file("ensembleReady/dongleiSNPs.txt")
writeLines(as.character(dongleiSNPs), fileConn)
close(fileConn)

sigRegReady = reg %>% select(Existing_variation)

#write regulom Ready format (just rsID numbers) to file
fileConn<-file("regulomReady/sigEnsemblMetaSnpsRSIDs.txt")
writeLines(as.character(sigRegReady$Existing_variation), fileConn)
close(fileConn)

##############################
# read in the significant snps from each population
mex = read.table("inputFromSkat/Result_mex_P_1en5.txt", header = T, sep = "\t")
pr = read.table("inputFromSkat/Result_pur_P_1en5.txt", header = T, sep = "\t")
aa = read.table("inputFromSkat/Result_afr_P_1en5.txt", header = T, sep = "\t")
#split the snps to make ensemble ready
mexSplit = separate(mex, into = c("chr","bp","ref","snp") ,col = SNP, sep= ":", remove = T)
mexEns = data.frame(mexSplit$chr, mexSplit$bp, mexSplit$bp,paste(mexSplit$ref,mexSplit$snp,sep = "/"))
names(mexEns) = c("chr","bpStart","bpEnd","ref/snp")
mexReady = paste(mexEns$chr,mexEns$bpStart,mexEns$bpEnd,mexEns$`ref/snp`, sep = " ")
#take the Location to get donglei genotypes
mexDong = paste(mex$SNP,sep="")
fileConn=file("mexDong.txt")
writeLines(mexDong,fileConn)
close(fileConn)

prSplit = separate(pr, into = c("chr","bp","ref","snp") ,col = SNP, sep= ":", remove = T)
prEns = data.frame(prSplit$chr, prSplit$bp, prSplit$bp,paste(prSplit$ref,prSplit$snp,sep = "/"))
names(prEns) = c("chr","bpStart","bpEnd","ref/snp")
prReady = paste(prEns$chr,prEns$bpStart,prEns$bpEnd,prEns$`ref/snp`, sep = " ")
#take location to get dong
prDong = paste(pr$SNP,sep="")
fileConn=file("prDong.txt")
writeLines(prDong,fileConn)
close(fileConn)

aaSplit = separate(aa, into = c("chr","bp","ref","snp") ,col = SNP, sep= ":", remove = T)
aaEns = data.frame(aaSplit$chr, aaSplit$bp, aaSplit$bp,paste(aaSplit$ref,aaSplit$snp,sep = "/"))
names(aaEns) = c("chr","bpStart","bpEnd","ref/snp")
aaReady = paste(aaEns$chr,aaEns$bpStart,aaEns$bpEnd,aaEns$`ref/snp`, sep = " ")
#take the Location to get donglei genotypes
aaDong = paste(aa$SNP,sep="")
fileConn=file("aaDong.txt")
writeLines(aaDong,fileConn)
close(fileConn)

# write ensembl ready format to file
fileConn<-file("ensemblReady/mexSnps.txt")
writeLines(mexReady, fileConn)
close(fileConn)

fileConn<-file("ensemblReady/prSnps.txt")
writeLines(prReady, fileConn)
close(fileConn)

fileConn<-file("ensemblReady/aaSnps.txt")
writeLines(aaReady, fileConn)
close(fileConn)

#read ensembl results
metaEnsResults = read.table("ensemblResults/sigMetaEnsemblTranscripts", sep = "\t", header = T,comment.char = "")
mexEnsResults = read.table("ensemblResults/mexEnsemblSnps", sep = "\t", header = T,comment.char = "")
aaEnsResults = read.table("ensemblResults/aaEnsemblSnps", sep = "\t", header = T,comment.char = "")
prEnsResults = read.table("ensemblResults/prEnsemblSnps", sep = "\t", header = T,comment.char = "")


#subset desired features
mexEnsResultsF = mexEnsResults %>% select(X.Uploaded_variation, Location,Consequence, Existing_variation,IMPACT,SYMBOL)
metaEnsResultsF = metaEnsResults %>% select(X.Uploaded_variation, Location,Consequence, Existing_variation,IMPACT,SYMBOL)
aaEnsResultsF = aaEnsResults %>% select(X.Uploaded_variation, Location,Consequence, Existing_variation,IMPACT,SYMBOL)
prEnsResultsF = prEnsResults %>% select(X.Uploaded_variation, Location,Consequence, Existing_variation,IMPACT,SYMBOL)


testF = mexEnsResultsF %>% group_by(Location,Existing_variation,X.Uploaded_variation= X.Uploaded_variation,
                                      Consequence=Consequence,IMPACT=IMPACT, SYMBOL=SYMBOL) %>% summarise(n = n())
# take only unique, take top hit
mexEnsResultsF = distinct(mexEnsResultsF,Location); mexEnsResultsF$pValue = mex$P; mexEnsResultsF$oddsR = mex$OR
metaEnsResultsF = distinct(metaEnsResultsF,Location); metaEnsResultsF$pValue = sigMeta$P; metaEnsResultsF$oddsR = sigMeta$OR
aaEnsResultsF = distinct(aaEnsResultsF,Location); aaEnsResultsF$pValue = aa$P; aaEnsResultsF$oddsR = aa$OR
prEnsResultsF = distinct(prEnsResultsF,Location); prEnsResultsF$pValue = pr$P; prEnsResultsF$oddsR = pr$OR

#save these for Marquitta:
save(mexEnsResultsF, file = "mexEnsResultsF.csv")
save(aaEnsResultsF, file = "aaEnsResultsF.csv")
save(prEnsResultsF, file = "prEnsResultsF.csv")
save(metaEnsResultsF, file = "metaEnsResultsF.csv")

#write rsid's to files for use in regulome or snpInfo
sigRegReady = reg %>% select(Existing_variation)

#write regulom Ready format (just rsID numbers) to file
mexFNoRSID = mexEnsResultsF %>% filter(Existing_variation=="-")
mexFRSID = mexEnsResultsF %>% filter(Existing_variation != "-")
fileConn<-file("regulomReady/mexFRSID.txt")
writeLines(as.character(mexFRSID$Existing_variation), fileConn)
close(fileConn)

aaFNoRSID = aaEnsResultsF %>% filter(Existing_variation=="-")
aaFRSID = aaEnsResultsF %>% filter(Existing_variation != "-")
fileConn<-file("regulomReady/aaFRSID.txt")
writeLines(as.character(aaFRSID$Existing_variation), fileConn)
close(fileConn)

prFNoRSID = prEnsResultsF %>% filter(Existing_variation=="-")
prFRSID = prEnsResultsF %>% filter(Existing_variation != "-")
fileConn<-file("regulomReady/prFRSID.txt")
writeLines(as.character(prFRSID$Existing_variation), fileConn)
close(fileConn)

metaFNoRSID = metaEnsResultsF %>% filter(Existing_variation=="-")
metaFRSID = metaEnsResultsF %>% filter(Existing_variation != "-")
fileConn<-file("regulomReady/metaFRSID.txt")
writeLines(as.character(metaFRSID$Existing_variation), fileConn)
close(fileConn)

#why isn't snpInfo running every rsid?
# Get RSID's that were in the input, but weren't returned as output
aaWeirdStuff = read.csv("RSIDs/RSIDresults/aaSNPinfo.csv", header = T)
aaWeirdStuff = aaWeirdStuff %>% select(Existing_variation = rs)
missingAA = anti_join(aaEnsResultsF, aaWeirdStuff)
###TASKS###
# 1.which snps are public/private: do a join between the groups based on rsID number (mexOnly, aaONly, etc)
# RESULT: all ancestral snps are private
aaJpr = rbind(aaFRSID,prFRSID)
mexOnly = setdiff(mexFRSID$Existing_variation,aaJpr$Existing_variation)
length(mexOnly)/length(mexFRSID$Existing_variation) #1
mexJpr = rbind(mexFRSID,prFRSID)
aaOnly = setdiff(aaFRSID$Existing_variation, mexJpr$Existing_variation)
length(aaOnly)/length(aaFRSID$Existing_variation) #1
mexJaa = rbind(mexFRSID,aaFRSID)
prOnly = setdiff(prFRSID$Existing_variation,mexJaa$Existing_variation)
length(prOnly)/length(prFRSID$Existing_variation) #1
# 2. what are the top snps from each population, including meta. arrange by p-value, take top n-number from each group,
# combine those into a new df. 

#read in the plink genotypes for each individual and transpose
dongChr5 = read.table("chr5.vcf", sep="\t", comment.char = "", skip = 6, header = T)
dongChr51 = dongChr5 %>% select(ID, X0_NWD100822, X0_NWD101012)
dong = t(dongChr51)
colnames(dong) = dong[1,]
dong = dong[-1,]
dong = as.data.frame(dong)
dong = cbind(Person = rownames(dong), dong)
rownames(dong) = NULL

######
mergedPeopleMeta = read.table("mergedPeople.vcf", sep = "\t",comment.char = "", header = T)
mergedPeopleMeta = mergedPeopleMeta %>% select(-c(X.CHROM,POS,REF,ALT,QUAL,FILTER,INFO,FORMAT))
mergedPeopleMeta = t(mergedPeopleMeta)
colnames(mergedPeopleMeta) = mergedPeopleMeta[1,]
mergedPeopleMeta = mergedPeopleMeta[-1,]
mergedPeopleMeta = as.data.frame(mergedPeopleMeta)
mergedPeopleMeta = cbind(personID = rownames(mergedPeopleMeta),mergedPeopleMeta)
rownames(mergedPeopleMeta) = NULL
dim(mergedPeopleMeta) #1484 x 206

#there were 2 people excluded for QC issues, must remove them
#NWD425395 ,NWD649684
mergedPeopleMeta = mergedPeopleMeta %>% filter(personID != "X0_NWD425395", personID != "X0_NWD649684") 

# necessary vars:asthma_hospital_12mo, asthma_hospital_12mo_count, asthma_er_12mo, asthma_er_12mo_count, steroids_more2wks,steroids_lastyear
# age, sex, bmi, ethnicity, BDR group
galaP = read.csv("phenotype/GALAII_WGS_Final_Sample_List_2015_12_data.csv", header = T)
sageP = read.csv("phenotype/SAGEII_WGS_Final_Sample_List_2015_12_data.csv", header = T)
#remove the pesky repeated person
#dups = sageP %>% mutate(personID = paste("X0", NWD_ID, sep = "_")) %>% filter(personID %in% mergedPeopleMeta$personID)
sageP = sageP[-243,]
#sagePsm = sageP %>% mutate(personID = paste("X0", NWD_ID, sep = "_")) %>% filter(personID %in% mergedPeopleMeta$personID) %>% select(personID,asthma_hospital_12mo, asthma_hospital_12mo_count, asthma_er_12mo, asthma_er_12mo_count, steroids_more2wks)
#galaPsm = galaP %>% mutate(personID = paste("X0", NWD_ID, sep = "_")) %>% filter(personID %in% mergedPeopleMeta$personID) %>% select(personID,asthma_hospital_12mo, asthma_hospital_12mo_count, asthma_er_12mo, asthma_er_12mo_count, steroids_more2wks)

#get personID and actual BDR response
sageResponse = sageP %>% mutate(personID = paste("X0", NWD_ID, sep = "_")) %>% filter(personID %in% mergedPeopleMeta$personID) %>% select(personID,delta1,delta2) %>% transform(delta2 = ifelse(is.na(delta2),0,delta2))
#sageR = sageResponse %>% mutate(maxDelta = ifelse(delta1 > delta2, delta1,delta2))

galaResponse = galaP %>% mutate(personID = paste("X0", NWD_ID, sep = "_")) %>% filter(personID %in% mergedPeopleMeta$personID) %>% select(personID,delta1,delta2) %>% transform(delta2 = ifelse(is.na(delta2),0,delta2))
#galaR = galaResponse %>% mutate(maxDelta = ifelse(delta1 > delta2, delta1,delta2))
#combine gala and sage vars into 1 table
galaSage = rbind(galaR,sageR) %>% select(personID, maxDelta)

pheno = read.table("phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt",sep = "\t", header = T)
phenoR = pheno %>% mutate(personID = paste("X0", ID, sep = "_")) %>% select(personID, age, sex, bmicat.fix,max.delta)

# join all the pheno type data on personID
#allPheno = inner_join(galaSage,phenoR)
# join the genetic and phenotype data
allFeaturesResponse = inner_join(phenoR, mergedPeopleMeta)

save(allFeaturesResponse, file = "allFeaturesResponse.RData")


#ToDo
#1.check excerbations with Neeta. 2: repeat gala processes on sage and column bind the 2 sets. 3: join the combined gala-sage set with the merged metaPeople on personID
#4. read in phenotype/merge.wgs.phenotype.1484_2015_1116_lite_pc_global_mod-bmicat_ordinal_fix-bmicat-pct_1482.txt and extract bmi etc. then join with complete data from 3, on ID