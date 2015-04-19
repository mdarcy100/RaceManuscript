#Note: 4/2015.  this is the code for the "primary" analysis (finding differentially associated genes)
# 4/3/2013
#read in UNC337, Normal and look at race associated genes
# supervised analysis by race
# will also look at confounding variables for table 1 in adjacent normal, normal

library(Biobase)
library(gdata) #
library(hgug4112a.db)
library(gplots)
library(limma)
library(gage)
library(survival)
library(qvalue)

source("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/pcaFuncs.R")
source("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts.R")

#############################################################################  set working directory
setwd("/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper")
##############################################################################  load UNC337 data that contains race
load("/Users/mdarcy100/Desktop/MTroester/AgePaper/UNC337/all_337_race_eset.RData") ## the UNC 337 samples.R
dim(all.337.race.eset)
data(egSymb)

entrezIds <- mget(featureNames(all.337.race.eset), envir=hgug4112aENTREZID,ifnotfound=NA)
haveEntrezIds <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(all.337.race.eset)) - length(haveEntrezIds)

fData(all.337.race.eset)$entrez<-entrezIds

symbols<-mget(featureNames(all.337.race.eset), envir=hgug4112aSYMBOL,ifnotfound=NA)

fData(all.337.race.eset)$symbol<-symbols
sub.337.eset <- all.337.race.eset[haveEntrezIds,]
dim(sub.337.eset)


table(sub.337.eset$Race)

# AA   C 
# 57 108 


## Regression on RACE in ALL cases
iqr.unc337.eset <- esApply(sub.337.eset,1,IQR)
design.race<- model.matrix(~ factor(sub.337.eset$Race))
colnames(design.race) <- c("Intercept","RACE")

race.limma.fit1 <- lmFit(sub.337.eset,design.race)
race.limma.ebayes <- eBayes(race.limma.fit1)

## Apply IQR cutoff of
race.limma.filtered <- race.limma.ebayes[iqr.unc337.eset>median(iqr.unc337.eset),]
##MD checking dimension
dim(race.limma.filtered)
#[1] 4976    2

#
race.limma.summary <- topTable(race.limma.filtered,coef="RACE",adjust.method="fdr",num=Inf)

colnames(race.limma.summary)
#[1] "entrez"    "symbol"    "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"          
 

## Distribution of p-values
pval.summary(race.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     21     65    266    474    741   1156   1500   1782   2047   4976 

pval.summary(race.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     3      4     12     25     39     98    199    294    449   4976 

q.value.race <- qvalue(p=race.limma.summary$P.Value)
summary(q.value.race,cuts=c(1e-4,0.001,0.01,0.025,0.05,0.10,0.15,0.20,0.25,1))
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
# p-value     21     65   266    474   741 1156  1500 1782  2047 4976
#q-value      3      6    14     27    60  167   324  521   779 4976


# # ########################### 
sig.race <- race.limma.summary[race.limma.summary$adj.P.Val < 0.05,]
class(sig.race)
colnames(sig.race)




########################################################################
# we also want to write the cluster file to know what direction everything is headed
# 4/17/2013
length(intersect(sig.race$entrez,fData(sub.337.eset)$entrez))
fData(sub.337.eset)[which(duplicated(fData(sub.337.eset))),]$entrez
racecluster.337.eset<-sub.337.eset[which(fData(sub.337.eset)$entrez%in%sig.race$entrez),]
#order by race status
racecluster.337.eset<-racecluster.337.eset[,order(racecluster.337.eset$Race)]
colnames(exprs(racecluster.337.eset))<-paste(racecluster.337.eset$Race,colnames(exprs(racecluster.337.eset))," ")
rownames(exprs(racecluster.337.eset))<-fData(racecluster.337.eset)$symbol


write.table(exprs(racecluster.337.eset),"Paper/RaceTumorCluster4172013.txt",sep="\t",col.names=NA)
########################################################################


# WE NEED THE LINE BELOW TO WRITE OUT THE FILE - THERE WERE SOME ISSUES
sig.race<-as.matrix(sig.race)
#sig.race <- lapply(sig.race, function(x) if(is.list(x)) unlist(x))

class(apply(exprs(racecluster.337.eset),1,sd))
write.table(sig.race,"SigGenesRace.txt",sep="\t",col.names=NA)
duplicated(rownames(as.matrix(apply(exprs(racecluster.337.eset),1,sd))))
merge(sig.race,as.matrix(apply(exprs(racecluster.337.eset),1,sd)),by.x='symbol', by.y='row.names')
write.table(cbind(sig.race,as.vector(apply(exprs(racecluster.337.eset),1,sd))),"SigGenesRace.txt",sep="\t",col.names=NA)
write.table(as.matrix(apply(exprs(racecluster.337.eset),1,sd)),"RaceSD.txt",sep="\t",col.names=NA)


# # ########################### # ############################ # ############################ stratify by subtyoe
# ############################  Luminal A
LumA <- sub.337.eset[,which(sub.337.eset$PAM50.Call %in% c("LumA"))]
dim(LumA)
#Features  Samples 
#    9953       68 

## Regression on RACE in Luminal A cases
iqr.luma.eset <- esApply(LumA,1,IQR)
design.race<- model.matrix(~ factor(LumA$Race))
colnames(design.race) <- c("Intercept","RACE")

race.limma.fit1 <- lmFit(LumA,design.race)
race.limma.ebayes <- eBayes(race.limma.fit1)

## Apply IQR cutoff of
race.limma.filtered <- race.limma.ebayes[iqr.luma.eset>median(iqr.luma.eset),]
##MD checking dimension
dim(race.limma.filtered)
#[1] 4976    2

#
race.limma.summary <- topTable(race.limma.filtered,coef="RACE",adjust.method="fdr",num=Inf)

##MD did this to figure out what the data structure looked like:
colnames(race.limma.summary)
#[1] "entrez"    "symbol"    "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"          
 

## Distribution of p-values
pval.summary(race.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     10     39    192    341    523    884   1188   1484   1745   4976 

pval.summary(race.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     1      1      4      6     10     23     62    137    184   4976 

q.value.race <- qvalue(p=race.limma.summary$P.Value)
summary(q.value.race,cuts=c(1e-4,0.001,0.01,0.025,0.05,0.10,0.15,0.20,0.25,1))
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#p-value     10     39   192    341   523  884  1188 1484  1745 4976
#q-value      1      1     4      7    11   47   135  192   270 4976


sig.race <- race.limma.summary[race.limma.summary$adj.P.Val < 0.10,]

# we also want to write the cluster file to know what direction everything is headed
# 4/17/2013
lumAcluster.337.eset<-sub.337.eset[which(fData(sub.337.eset)$symbol%in%sig.race$symbol),]
dim(lumAcluster.337.eset)
#Features  Samples 
#      23      165 

fData(lumAcluster.337.eset)$symbol
#order by race status
lumAcluster.337.eset<-lumAcluster.337.eset[,order(lumAcluster.337.eset$Race)]

###########################4/25/2013 - writing out table for the paper with expression in significant genes
# MD: 6/25/2013 - MD chaged this because was looking at the expression over all tumor types (lumAcluster.337.eset).  - now print out for just LuminalA tumors
#get the overall median expression for each gene
lumA.337.eset<-sub.337.eset[which(fData(sub.337.eset)$symbol%in%sig.race$symbol),which(sub.337.eset$PAM50.Call %in% c("LumA"))]
median_all <- esApply(lumA.337.eset,1,median)
# the median expression in the AA
median_AA<-esApply(lumA.337.eset[,which(pData(lumA.337.eset)$Race%in%("AA"))],1,median)
# the median expression in caucasians
median_C<-esApply(lumA.337.eset[,which(pData(lumA.337.eset)$Race%in%("C"))],1,median)
LumARaceTable<-cbind(fData(lumA.337.eset)$symbol,fData(lumA.337.eset)$entrez,median_all,median_AA,median_C)
colnames(LumARaceTable)<-c("Symbol","Entrez","MedianLog2rg_All","MedianLog2rg_AA","MedianLog2rg_C")
write.table(LumARaceTable,"Tables/LumARaceTable_June252013.txt",sep="\t",col.names=NA)




###########################4/17/2013
colnames(exprs(lumAcluster.337.eset))<-paste(lumAcluster.337.eset$Race,colnames(exprs(lumAcluster.337.eset))," ")
rownames(exprs(lumAcluster.337.eset))<-fData(lumAcluster.337.eset)$symbol
write.table(exprs(lumAcluster.337.eset),"Paper/RaceLumACluster4172013.txt",sep="\t",col.names=NA)
########################################################################

sig.race<-as.matrix(sig.race)
#sig.race <- lapply(sig.race, function(x) if(is.list(x)) unlist(x))
write.table(sig.race,"Paper/SigGenesRaceLumA_Apr92013.txt",sep="\t",col.names=NA)




# ############################ # ############################ # ############################ 
# Basal subtype

Basal<- sub.337.eset[,which(sub.337.eset$PAM50.Call %in% c("Basal"))]
dim(Basal)
#Features  Samples 
#    9953       39 


## Regression on Race in Basal cases
iqr.basal.eset <- esApply(Basal,1,IQR)
design.race<- model.matrix(~ factor(Basal$Race))
colnames(design.race) <- c("Intercept","RACE")

race.limma.fit1 <- lmFit(Basal,design.race)
race.limma.ebayes <- eBayes(race.limma.fit1)

## Apply IQR cutoff of
race.limma.filtered <- race.limma.ebayes[iqr.basal.eset>median(iqr.basal.eset),]
##MD checking dimension
dim(race.limma.filtered)
#[1] 4976    2

#
race.limma.summary <- topTable(race.limma.filtered,coef="RACE",adjust.method="fdr",num=Inf)

##MD did this to figure out what the data structure looked like:
colnames(race.limma.summary)
#[1] "entrez"    "symbol"    "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"          
 

## Distribution of p-values
pval.summary(race.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     2      3     45    100    223    455    711    986   1246   4976 

pval.summary(race.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     0      0      1      2      2      2      2      2      2   4976 

q.value.race <- qvalue(p=race.limma.summary$P.Value)
summary(q.value.race,cuts=c(1e-4,0.001,0.01,0.025,0.05,0.10,0.15,0.20,0.25,1))
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
# p-value      2      3    45    100   223  455   711  986  1246 4976
#q-value      0      0     1      2     2    2     2    2     2 4976

# # ########################### 
sig.race <- race.limma.summary[race.limma.summary$adj.P.Val < 0.10,]
class(sig.race)
colnames(sig.race)

# we also want to write the cluster file to know what direction everything is headed
# 4/17/2013
Basalcluster.337.eset<-sub.337.eset[which(fData(sub.337.eset)$symbol%in%sig.race$symbol),]
dim(Basalcluster.337.eset)
#Features  Samples 
#      2      165 

fData(Basalcluster.337.eset)$symbol
#order by race status
Basalcluster.337.eset<-Basalcluster.337.eset[,order(Basalcluster.337.eset$Race)]


###########################4/25/2013 - writing out table for the paper with expression in significant genes
#get the overall median expression for each gene
# MD: 6/25/2013 - MD chaged this because was looking at the expression over all tumor types (Basalcluster.337.eset).  - now print out for just LuminalA tumors
Basal.337.eset<-sub.337.eset[which(fData(sub.337.eset)$symbol%in%sig.race$symbol),which(sub.337.eset$PAM50.Call %in% c("Basal"))]
median_all <- esApply(Basal.337.eset,1,median)
# the median expression in the AA
median_AA<-esApply(Basal.337.eset[,which(pData(Basal.337.eset)$Race%in%("AA"))],1,median)
# the median expression in caucasians
median_C<-esApply(Basal.337.eset[,which(pData(Basal.337.eset)$Race%in%("C"))],1,median)
BasalRaceTable<-cbind(fData(Basal.337.eset)$symbol,fData(Basal.337.eset)$entrez,median_all,median_AA,median_C)
colnames(BasalRaceTable)<-c("Symbol","Entrez","MedianLog2rg_All","MedianLog2rg_AA","MedianLog2rg_C")
write.table(BasalRaceTable,"Tables/BasalRaceTable_June252013.txt",sep="\t",col.names=NA)
###########################


#########4/17/2013 

colnames(exprs(Basalcluster.337.eset))<-paste(Basalcluster.337.eset$Race,colnames(exprs(Basalcluster.337.eset))," ")
rownames(exprs(Basalcluster.337.eset))<-fData(Basalcluster.337.eset)$symbol
write.table(exprs(Basalcluster.337.eset),"Paper/RaceBasalCluster4172013.txt",sep="\t",col.names=NA)



#######4/9/2013 
sig.race<-as.matrix(sig.race)
#sig.race <- lapply(sig.race, function(x) if(is.list(x)) unlist(x))
write.table(sig.race,"Paper/SigGenesRaceBasal_Apr92013.txt",sep="\t",col.names=NA)
########################################################################




# # ########################### look at confounders (age, bmi, race)
# also see comments in the CRYBB2_Analysis2.R for some of other associations

table(factor(floor(sub.337.eset$Age.x/10)))
# 2  3  4  5  6  7  8 
# 2 15 41 40 22 27 14 

table(sub.337.eset$Race)
# AA   C 
# 57 108

table()
