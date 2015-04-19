#NKI, UNC337 analyses in beginning;
#Survival analysis 463
#METBARIC analyses begin ~ line 590
#METBARIC survival 910
#younger women begin on line 1050
#metbaric old 1256
#supplemental code for the race paper.
# document the association between disparity points and proliferation score and survival 
#i

# as per chuck's comments on 12/20/2014:
#MVA that includes size, node, age, grade, intrinsic subtype, proliferation score, and disparity score
# he also wants it in a table with univariate, multivariate HR and proliferation scores
library(Biobase)
library(gdata) #
library(hgug4112a.db)
library(gage)
library(survival)
library(muhaz)
data(egSymb)

source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts.R")

race_symbols<-c("CRYBB2","MUC1", "ACOX2", "SQLE", "TYMS", "PSPH")
race_entrez<-sym2eg(race_symbols)

#MD - update 11/2014.  calculate progression score 
proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
proliferation_entrez<-sym2eg(proliferationGenes)

#  set working directory
setwd("/Users/mdarcy/Desktop/UNC337_Race/RCode/")

# first do in the training data.

#############################################################
#MD 1/1/2014: find node, size, grade & age variables in both datasets
# evaluate the association with 'proliferation score'
#update: 12/23/2014. 
# adjust for: size, node, age, grade, intrinsic subtype, proliferation score, and disparity score
#############################################################
#first median center each dataset - then combine the datasets
####################################################################################################
# read in UNC337 data
load("/Users/mdarcy/Desktop/MTroester/AgePaper/UNC337/all_337_race_eset.RData") ## the UNC 337 samples.R
dim(all.337.race.eset)
data(egSymb)

entrezIds <- mget(featureNames(all.337.race.eset), envir=hgug4112aENTREZID,ifnotfound=NA)
haveEntrezIds <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(all.337.race.eset)) - length(haveEntrezIds)

fData(all.337.race.eset)$entrez<-entrezIds

symbols<-mget(featureNames(all.337.race.eset), envir=hgug4112aSYMBOL,ifnotfound=NA)

fData(all.337.race.eset)$symbol<-symbols

sub.337.eset <- all.337.race.eset[haveEntrezIds,]
race.337.eset <- sub.337.eset[featureData(sub.337.eset)$symbol %in% unique(c(race_symbols,proliferationGenes)),]
dim(race.337.eset)
#Features  Samples (there are only 5 proliferation genes) 
#       10      165 

rownames(exprs(race.337.eset))<-fData(race.337.eset)$symbol

race.337.eset$race2<-race.337.eset$Race
table(race.337.eset$race2)

#################### remove those samples without missing survival data
notmissing.ind<-which(sapply(race.337.eset$Overall.Survival.Event..0.alive..1.DOD.or.DOC..x,function(x) !is.na(x)))
length(notmissing.ind)
race.337.eset<-race.337.eset[,notmissing.ind]

race.337.eset$S_Event<-race.337.eset$Overall.Survival.Event..0.alive..1.DOD.or.DOC..x
race.337.eset$S_Time<-race.337.eset$Overall.suvival.months.x/12
dim(race.337.eset)


#create size variable
table(race.337.eset$Size..1....2cm..2...2cm.to...5cm..3..5cm..4.any.size.with.direct.extension.to.chest.wall.or.skin..x)
#  1 1.5   2   3   4 
# 42   1  72  24  13 
race.337.eset$Size<-ifelse(race.337.eset$Size..1....2cm..2...2cm.to...5cm..3..5cm..4.any.size.with.direct.extension.to.chest.wall.or.skin..x >= 2,1,
                           ifelse(race.337.eset$Size..1....2cm..2...2cm.to...5cm..3..5cm..4.any.size.with.direct.extension.to.chest.wall.or.skin..x < 2, 0,NA)) #useNA="always")
table(race.337.eset$Size)
#  0   1 
# 43 109 

#create age variable
age_cat_c<-cut(race.337.eset$Age.y,breaks=c(0,30,40,50,60,70,100),right=FALSE) # corresponds to 
f.age.c<-factor(age_cat_c)
race.337.eset$f.age.c<-f.age.c
table(f.age.c)
#  [0,30)  [30,40)  [40,50)  [50,60)  [60,70) [70,100) 
#       2       15       40       36       22       40 

#create node variable
table(race.337.eset$Node.status..1.positive.1.or.more.nodes...0.negative..y)
race.337.eset$Node<-race.337.eset$Node.status..1.positive.1.or.more.nodes...0.negative..y
table(race.337.eset$Node)
# 0  1  2 
#78 74  1 
#combine 1/2 together
race.337.eset$Node<-ifelse(race.337.eset$Node > 1, 1, race.337.eset$Node)

####create more simple ER variable "ER...1.positive..0.negative..x" 
table(race.337.eset$ER...1.positive..0.negative..x)
# 0  1 
#58 94 
race.337.eset$ER<-race.337.eset$ER...1.positive..0.negative..x
table(race.337.eset$ER)

##### create more simple PR variable PGR...1.positive..0.negative.
race.337.eset$PR<-race.337.eset$PGR...1.positive..0.negative.
table(race.337.eset$PR)
# 0  1 
#67 67 

######create grade variable
table(race.337.eset$Grade)
# 1  2  3 
#14 53 80 

table(race.337.eset$PAM50.Call)
# Basal   Her2   LumA   LumB Normal 
#    35     17     66     26     11 
colnames(pData(race.337.eset))

#now create the proliferation score
proliferation.idx<-which(rownames(exprs(race.337.eset))%in%proliferationGenes)
proliferation.idx
prolif_score<-apply(exprs(race.337.eset)[proliferation.idx,],2,function(x){sum(x,na.rm = TRUE)})
pData(race.337.eset)$prolif_score<-prolif_score
hist(pData(race.337.eset)$prolif_score)

# median center the data
medians<- esApply(race.337.eset,1,median)
unc337.exprs.scale <- t(scale(t(exprs(race.337.eset)),center=medians,scale=FALSE))

## Create a new eset
unc337.scale.eset<-race.337.eset
exprs(unc337.scale.eset) <- unc337.exprs.scale
dim(unc337.scale.eset)
##10 155
colnames(pData(race.337.eset))


####################################################################################################
# read in NKI data
load("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/nki295_eset.RData") ## process_NKI295samples.R
dim(nki295.eset)
nki295.eset <- nki295.eset[!is.na(featureData(nki295.eset)$geneid),]
dim(nki295.eset)
#Features  Samples 
#   17600      295 

race.nki295.eset <- nki295.eset[fData(nki295.eset)$symbol %in% unique(c(race_symbols,proliferationGenes)),]
rownames(exprs(race.nki295.eset))
fData(race.nki295.eset)$symbol
#[1] "SQLE"   "MKI67"  "MUC1"   "PTTG1"  "ACOX2"  "PSPH"   "KNTC2"  "UBE2C"  "CCNB1"  "RRM2"   "CEP55"  "RRM2"  
#[13] "TYMS"   "CDCA1"  "BIRC5"  "CRYBB2" "CDC20" 


#there are duplicates
exprs(race.nki295.eset)<-collapseIDs(exprs(race.nki295.eset),fData(race.nki295.eset)$symbol,"mean")
dim(exprs(race.nki295.eset))
#[1]  16 295

colnames(pData(race.nki295.eset))
max.col <- apply(pData(race.nki295.eset)[,37:41], 1, which.max) 
#[16] "diameter(mm)"                           "T1_T2"                                  "Lymph_node_number_postive"             
#[19] "pN_3_classes"                           "Mastectomy"                             "ER"                                    
#[22] "Grade_3_classes"                        "Age(years)"                             "Chemo"       

pData(race.nki295.eset)$subtype<-names(pData(race.nki295.eset)[,37:41])[max.col] #43

#just remove the .Cor from beginning of 
pData(race.nki295.eset)$subtype<-gsub("Cor.","",pData(race.nki295.eset)$subtype)
pData(race.nki295.eset)$subtype<-ifelse(pData(race.nki295.eset)$subtype=="ERBB2","Her2",pData(race.nki295.eset)$subtype) 

race.nki295.eset$S_Time<-race.nki295.eset$"survival(death)"
race.nki295.eset$S_Event<-race.nki295.eset$"event_death" 

race.nki295.eset$race2<-"C"

# categorize age variable
mean(race.nki295.eset$'Age(years)')
age_cat_c<-cut(race.nki295.eset$'Age(years)',breaks=c(0,30,40,50,60,70,100),right=FALSE) # corresponds to 
f.age.c<-factor(age_cat_c)
race.nki295.eset$f.age.c<-f.age.c
table(f.age.c)
# [0,30) [30,40) [40,50) [50,60) 
#      4      59     183      49 

#create size variable
mean(race.nki295.eset$'diameter(mm)')
# [1] 22.53898
race.nki295.eset$Size<-ifelse(race.nki295.eset$'diameter(mm)'>=20,1,0)
table(race.nki295.eset$Size)
#  0   1 
#114 181 

#create node variable
table(race.nki295.eset$Lymph_node_number_postive,useNA="always")
#  0   1   2   3   4   5   6   7   8   9  11  13 <NA>
#151  58  27  21  14   7   4   3   3   4   2   1   0
race.nki295.eset$Node<-ifelse(race.nki295.eset$Lymph_node_number_postive > 0, 1, 0)
table(race.nki295.eset$Node)
#  0   1 
# 151 144 

#create grade variable
table(race.nki295.eset$Grade_3_classes,useNA="always")
# Intermediate  Poorly diff    Well diff         <NA> 
#         101          119           75            0  
race.nki295.eset$Grade<-ifelse(race.nki295.eset$Grade_3_classes=="Intermediate",2,
                               ifelse(race.nki295.eset$Grade_3_classes=="Well diff",1,3))

table(race.nki295.eset$Grade)
#  1   2   3 
# 75 101 119 
dim(race.nki295.eset)
#Features  Samples 
#       16      295 

table(pData(race.nki295.eset)$subtype)
# Basal   Her2   LumA   LumB Normal 
#    46     49     86     81     33

#now create the proliferation score
proliferation.idx<-which(rownames(exprs(race.nki295.eset))%in%proliferationGenes)
proliferation.idx
prolif_score<-apply(exprs(race.nki295.eset)[proliferation.idx,],2,function(x){sum(x,na.rm = TRUE)})
pData(race.nki295.eset)$prolif_score<-prolif_score
hist(pData(race.nki295.eset)$prolif_score)

# median center the data
medians<- esApply(race.nki295.eset,1,median)
nki295.exprs.scale <- t(scale(t(exprs(race.nki295.eset)),center=medians,scale=FALSE))

## Create a new eset
nki295.scale.eset<-race.nki295.eset
exprs(nki295.scale.eset) <- nki295.exprs.scale
dim(nki295.scale.eset)
#Features  Samples 
#      16      295 

####################################################################################################
#create new expression set with combined data, 
# need race variable
# only survival data for now - now also need node, size, age, grade

row.names<-rbind(as.matrix(colnames(exprs(nki295.scale.eset))),as.matrix(colnames(exprs(unc337.scale.eset))))

STime_combined<-rbind(as.matrix(nki295.scale.eset$S_Time),
                      as.matrix(unc337.scale.eset$S_Time))

SEvent_combined<-rbind(as.matrix(nki295.scale.eset$S_Event),
                       as.matrix(unc337.scale.eset$S_Event))

race_combined<-rbind(as.matrix(nki295.scale.eset$race2),
                     as.matrix(unc337.scale.eset$race2))

subtype<-rbind(as.matrix(nki295.scale.eset$subtype),
               as.matrix(unc337.scale.eset$PAM50.Call))

grade<-rbind(as.matrix(nki295.scale.eset$Grade),
             as.matrix(unc337.scale.eset$Grade))

node<-rbind(as.matrix(nki295.scale.eset$Node),
            as.matrix(unc337.scale.eset$Node))

size<-rbind(as.matrix(nki295.scale.eset$Size),
            as.matrix(unc337.scale.eset$Size))

ER<-rbind(as.matrix(nki295.scale.eset$ER),
          as.matrix(unc337.scale.eset$ER))


age_cat<-rbind(as.matrix(nki295.scale.eset$f.age.c),
               as.matrix(unc337.scale.eset$f.age.c))

age<-rbind(as.matrix(nki295.scale.eset$'Age(years)'),
           as.matrix(unc337.scale.eset$Age.y))

prolif_score<-rbind(as.matrix(nki295.scale.eset$prolif_score),
                    as.matrix(unc337.scale.eset$prolif_score))

data_source<-rbind(as.matrix(rep("NKI",295)),as.matrix(rep("UNC337",155)))	


combined<-cbind(data_source,STime_combined,SEvent_combined,race_combined,subtype,grade,node,size,ER,age_cat,age,prolif_score)
rownames(combined)<-row.names
colnames(combined)<-c("Data_Source","S_Time","S_Event","Race","subtype","grade","node","size","ER","age_cat","age","prolif_score")

patient.phenotype <- data.frame(combined)
hist(as.numeric(as.character(patient.phenotype$S_Time)))
combined.phenoData <- new("AnnotatedDataFrame", data = patient.phenotype)

#combined expression data

two.data<-merge(exprs(nki295.scale.eset),exprs(unc337.scale.eset),by='row.names', all.x=TRUE, all.y=TRUE, incomparables=NA)
rownames(two.data)<-two.data$Row.names
head(colnames(two.data))
two.data<-two.data[,-1] # remove the first column


dim(two.data)

#[1] 16 450

race.combinedtraining.eset <- new("ExpressionSet", exprs = data.matrix(two.data),
                                 phenoData = combined.phenoData, annotation = "hgug4112a")

save(race.combinedtraining.eset,file="race_combinedtraining_eset.RData")

########### ########### ########### ########### ########### ########### ###########
#calculate whether each of genes in disparity score
# is relatively up or down regulated
# proliferation score has already been calculated for each of the training dataset
########### ########### ########### ########### ########### ########### ###########
rownames(exprs(race.combinedtraining.eset))
# [1] "ACOX2"  "BIRC5"  "CCNB1"  "CDC20"  "CDCA1"  "CEP55"  "CRYBB2" "KNTC2"  "MKI67"  "MUC1"   "PSPH"   "PTTG1" 
#[13] "RRM2"   "SQLE"   "TYMS"   "UBE2C" 

#[ [1] "ACOX2"  "BIRC5"  "CCNB1"  "CDC20"  "CDCA1"  "CEP55"  "CRYBB2" "KNTC2"  "MKI67"  "MUC1"   "PSPH"   "PTTG1" 
#[13] "RRM2"   "SQLE"   "TYMS"   "UBE2C" 
#create a "ACOX2" "MUC1"  "SQLE"  "TYMS"  variable based on median centered data
# ACOX2 - up is good -1* exprs
# CYRYBB2 - up is bad 1*exprs
# MUC1 - up is good-1*exprs
# SQLE - up is bad 1*exprs
# TYMS - up is bad 1*exprs
#PSPH - up is bad 1*expr

ACOX2<-sapply(exprs(race.combinedtraining.eset)[1,],castDirection)
CRYBB2<-sapply(exprs(race.combinedtraining.eset)[7,],castDirection)
MUC1<-sapply(exprs(race.combinedtraining.eset)[10,],castDirection)
PSPH<-sapply(exprs(race.combinedtraining.eset)[11,],castDirection)
SQLE<-sapply(exprs(race.combinedtraining.eset)[14,],castDirection)
TYMS<-sapply(exprs(race.combinedtraining.eset)[15,],castDirection)

table(ACOX2)
table(CRYBB2)
table(MUC1)
table(SQLE)
table(TYMS)
table(PSPH)
# all are like PSPH
#PSPH
# -1   1 
#224 226 
race.combinedtraining.eset$ACOX2<-ACOX2
race.combinedtraining.eset$CRYBB2<-CRYBB2
race.combinedtraining.eset$MUC1<-MUC1
race.combinedtraining.eset$SQLE<-SQLE
race.combinedtraining.eset$PSPH<-PSPH
race.combinedtraining.eset$TYMS<-TYMS


colnames(pData(race.combinedtraining.eset))
# [1] "Data_Source"  "S_Time"       "S_Event"      "Race"         "subtype"      "grade"        "node"        
# [9] "ER"  "size"         "age_cat"      "prolif_score" "ACOX2"        "CRYBB2"       "MUC1"         "SQLE"    "PSPH"    
#[17]    "TYMS"     


######## now calculate disparity points and look at survival associations.
scale_gene <- matrix(c(-1,1,-1,1,1,1), ncol = 6)
d_score<-scale_gene%*%t(pData(race.combinedtraining.eset)[,which(colnames(pData(race.combinedtraining.eset))%in%c("ACOX2","CRYBB2","MUC1","SQLE","PSPH","TYMS"))])

quantile(d_score,na.rm=TRUE)
#  0%  25%  50%  75% 100% 
#  -6   -2    0    2    6 
#<2, >2 are good cutpoints

###################### cdiscrete score
# first evaluate the association between 'points' and proliferation score.
pData(race.combinedtraining.eset)$d_score<-t(d_score)
table(pData(race.combinedtraining.eset)$d_score)
#-6 -4 -2  0  2  4  6 
#29 63 80 97 97 54 30 
head(pData(race.combinedtraining.eset))

#NKI and UNC337 have different scales for proliferation scores since different number of genes present
#NKI
cor(as.numeric(as.character(pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("NKI"))])$prolif_score)),
    pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("NKI"))])$d_score,use="complete.obs")
# [1,] 0.572546

#UNC337
cor(as.numeric(as.character(pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("UNC337"))])$prolif_score)),
    pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("UNC337"))])$d_score,use="complete.obs")
#[1,] 0.5996354

cor(as.numeric(as.character(pData(race.combinedtraining.eset)$prolif_score)),
    pData(race.combinedtraining.eset)$d_score,use="complete.obs")
#[1,] 0.2883048

points_fit<-lm(as.numeric(as.character(pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("NKI"))])$prolif_score))~
                 pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("NKI"))])$d_score) 

summary(points_fit)
#pData(race.combinedtraining.eset[, which(pData(race.combinedtraining.eset)$Data_Source %in% c("NKI"))])$d_score  0.37143    0.03107   11.95
points_fit<-lm(as.numeric(as.character(pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("UNC337"))])$prolif_score))~
                 pData(race.combinedtraining.eset[,which(pData(race.combinedtraining.eset)$Data_Source%in%c("UNC337"))])$d_score) 

summary(points_fit)
#pData(race.combinedtraining.eset[, which(pData(race.combinedtraining.eset)$Data_Source %in% c("UNC337"))])$d_score   0.78905    0.08513   9.268

points_fit<-lm(as.numeric(as.character(pData(race.combinedtraining.eset)$prolif_score))~
                 pData(race.combinedtraining.eset)$d_score) 

summary(points_fit)
#pData(race.combinedtraining.eset)$d_score  0.51884    0.08141   6.373 4.62e-10 ***

pData(race.combinedtraining.eset)$prolif_score<-as.numeric(as.character(pData(race.combinedtraining.eset)$prolif_score))
dscore_cat_c<-cut(race.combinedtraining.eset$d_score,breaks=c(-7,-2,3,7),right=FALSE) # corresponds to quantiles
f.dscore.c<-factor(dscore_cat_c)

pData(race.combinedtraining.eset)$f.dscore.c<-f.dscore.c

table(pData(race.combinedtraining.eset)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#     92     274      84 

TrainingLumA<-race.combinedtraining.eset[,which(race.combinedtraining.eset$subtype%in%c("LumA"))]
######### what is the association between proliferation score and Race overall and in luminal A tumors
#this is the continuous variable
t.test(pData(TrainingLumA)$d_score~factor(pData(TrainingLumA)$Race))
#t = 7.2599, df = 23.959, p-value = 1.7e-07
#95 percent confidence interval:
#  2.933961 5.265044
#mean in group AA  mean in group C 
#         1.666667        -2.432836 


table(pData(TrainingLumA)$Data_Source)
#AA only in UNC337 data (NKI predominately CAU), proliferation scores not on the same scale
#UNC337
t.test(pData(TrainingLumA[,which(TrainingLumA$Data_Source%in%c("UNC337"))])$prolif_score~factor(pData(TrainingLumA[,which(TrainingLumA$Data_Source%in%c("UNC337"))])$Race))
#t = 4.3175, df = 57.324, p-value = 6.324e-05
#95 percent confidence interval:
#1.563935 4.268828
#mean in group AA  mean in group C 
#       -11.95189        -14.86827 


table(TrainingLumA$Race)
#AA   C 
#18 134 

#how does proliferation score vary overall by race in unc337 dataset
t.test(pData(race.combinedtraining.eset[,which(race.combinedtraining.eset$Data_Source%in%c("UNC337"))])$prolif_score~factor(pData(race.combinedtraining.eset[,which(race.combinedtraining.eset$Data_Source%in%c("UNC337"))])$Race))
#t = 3.6891, df = 139.106, p-value = 0.0003219
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1.095857 3.627085
#sample estimates:
#  mean in group AA  mean in group C 
#      -9.332509       -11.693980 

##########################################################################
#Survival analyses: 
######################################################## 
#Discrete (points) continuous (cox proportional) model
# IN TRAINING DATA: UNC337+ NKI296
#-6 to 6 points
race.combinedtraining.eset$S_Time<-as.numeric(as.character(race.combinedtraining.eset$S_Time))
race.combinedtraining.eset$S_Event<-as.numeric(as.character(race.combinedtraining.eset$S_Event))
all.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ race.combinedtraining.eset$d_score)

summary(all.surv.death)

#race.combinedtraining.eset$d_score 0.17398   1.19003  0.03151 5.521 3.38e-08 ***
#race.combinedtraining.eset$d_score      1.19     0.8403     1.119     1.266

#Likelihood ratio test= 32.01  on 1 df,   p=1.533e-08
#Wald test            = 30.48  on 1 df,   p=3.377e-08
#Score (logrank) test = 31.61  on 1 df,   p=1.88e-08


#HR 6 points/-6 points = exp(12*0.17398) 
##[1] 8.07
L95 = exp(12*0.17398-1.96*.03151)
L95
#[1] 7.583694
U95= exp(12*0.17398+1.96*0.03151)
U95
#[1] 8.580735
#HR 0 points/-6 points = exp(6*0.17398) = [1] 2.840216
L95 = exp(6*0.17398-1.96*0.03151)
L95
#[1] 2.670112
U95= exp(0.17398+1.96*0.03151)
U95
#1.265845

#in categories
all.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ race.combinedtraining.eset$f.dscore.c)
summary(all.surv.death)
#race.combinedtraining.eset$f.dscore.c[-2,3)     2.369     0.4222     1.215     4.617
#race.combinedtraining.eset$f.dscore.c[3,7)      5.273     0.1896     2.607    10.666

#adjust for grade, node, subtype, age, proliferation score, size
#adjust for disparity score and grade
grade.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$d_score+
                               race.combinedtraining.eset$grade)

summary(grade.adj.surv.death)

grade.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                                race.combinedtraining.eset$f.dscore.c+
                                race.combinedtraining.eset$grade)

summary(grade.adj.surv.death)
#race.combinedtraining.eset$f.dscore.c[-2,3)     1.447     0.6913     0.732     2.859
#race.combinedtraining.eset$f.dscore.c[3,7)      2.792     0.3582     1.350     5.772

#node isn't significant in the model
node.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                                race.combinedtraining.eset$d_score+
                                race.combinedtraining.eset$grade+
                                race.combinedtraining.eset$node)

summary(node.adj.surv.death)

# do in categories
node.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$f.dscore.c+
                               race.combinedtraining.eset$grade+
                               race.combinedtraining.eset$node)

summary(node.adj.surv.death)
#race.combinedtraining.eset$f.dscore.c[-2,3)     1.395     0.7170    0.7029     2.768
#race.combinedtraining.eset$f.dscore.c[3,7)      2.746     0.3641    1.3258     5.689

#add in size to the model
size.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$d_score+
                               race.combinedtraining.eset$grade+
                               race.combinedtraining.eset$node+
                               race.combinedtraining.eset$size)

summary(size.adj.surv.death)

#do in categories
size.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$f.dscore.c+
                               race.combinedtraining.eset$grade+
                               race.combinedtraining.eset$node+
                               race.combinedtraining.eset$size)

summary(size.adj.surv.death)
#race.combinedtraining.eset$f.dscore.c[-2,3) 0.32217   1.38011  0.35055 0.919 0.358083    
#race.combinedtraining.eset$f.dscore.c[3,7)  0.89733   2.45306  0.37478 2.394 0.016652 * 

#add in age to the model
race.combinedtraining.eset$age<-as.numeric(as.character(race.combinedtraining.eset$age))
age.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$d_score+
                               race.combinedtraining.eset$grade+
                               race.combinedtraining.eset$node+
                               race.combinedtraining.eset$size+
                               race.combinedtraining.eset$age)

summary(age.adj.surv.death)
#do in categories
age.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                              race.combinedtraining.eset$f.dscore.c+
                              race.combinedtraining.eset$grade+
                              race.combinedtraining.eset$node+
                              race.combinedtraining.eset$size+
                              race.combinedtraining.eset$age)

summary(age.adj.surv.death)


#add in proliferation score to the model

ps.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                              race.combinedtraining.eset$d_score+
                              race.combinedtraining.eset$grade+
                              race.combinedtraining.eset$node+
                              race.combinedtraining.eset$size+
                              race.combinedtraining.eset$age+
                              race.combinedtraining.eset$prolif_score)

summary(ps.adj.surv.death)

#do in categories
ps.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                             race.combinedtraining.eset$f.dscore.c+
                             race.combinedtraining.eset$grade+
                             race.combinedtraining.eset$node+
                             race.combinedtraining.eset$size+
                             race.combinedtraining.eset$age+
                             race.combinedtraining.eset$prolif_score)

summary(ps.adj.surv.death)
#race.combinedtraining.eset$f.dscore.c[-2,3)    1.3399     0.7463    0.6626     2.710
#race.combinedtraining.eset$f.dscore.c[3,7)     2.3646     0.4229    1.1024     5.072

# add in subtype
full.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                              race.combinedtraining.eset$d_score+
                              race.combinedtraining.eset$grade+
                              race.combinedtraining.eset$node+
                              race.combinedtraining.eset$size+
                              race.combinedtraining.eset$age+
                              race.combinedtraining.eset$subtype+
                              race.combinedtraining.eset$prolif_score)

summary(full.adj.surv.death)
# in categories
full.adj.surv.death <- coxph(Surv(race.combinedtraining.eset$S_Time,race.combinedtraining.eset$S_Event) ~ 
                               race.combinedtraining.eset$f.dscore.c+
                               race.combinedtraining.eset$grade+
                               race.combinedtraining.eset$node+
                               race.combinedtraining.eset$size+
                               race.combinedtraining.eset$age+
                               race.combinedtraining.eset$subtype+
                               race.combinedtraining.eset$prolif_score)

summary(full.adj.surv.death)


colnames(pData(race.combinedtraining.eset))



#now add clinical features: size, grade, node status
#################### remove those samples clinical data
#MD HERE NOW
table(race.combinedtraining.eset$subtype)

TrainingLumAB<-race.combinedtraining.eset[,which(race.combinedtraining.eset$subtype%in%c("LumA","LumB"))]
table(TrainingLumAB$subtype)
#unadjusted survival association in training dataset with disparity score (in points)
lum.surv.death <- coxph(Surv(TrainingLumAB$S_Time,
                             TrainingLumAB$S_Event) ~ TrainingLumAB$d_score)
summary(lum.surv.death)
#  n= 259, number of events= 41 

#coef exp(coef) se(coef)     z Pr(>|z|)   
#TrainingLumAB$d_score 0.12949   1.13824  0.04852 2.669  0.00762 **

#                      exp(coef) exp(-coef) lower .95 upper .95
#TrainingLumAB$d_score     1.138     0.8785     1.035     1.252

#calculate HR and 95% CI comparing +6 points to -6 points

#adjusted for clinical factors
adj.lum.surv.death <- coxph(Surv(TrainingLumAB$S_Time,TrainingLumAB$S_Event) ~ 
                              TrainingLumAB$d_score+
                              TrainingLumAB$grade+
                            TrainingLumAB$size+
                            TrainingLumAB$node)
summary(adj.lum.surv.death)

#coef exp(coef) se(coef)     z Pr(>|z|)   
#TrainingLumAB$d_score 0.03487   1.03548  0.05315 0.656  0.51175   
#TrainingLumAB$size1   0.68074   1.97533  0.42815 1.590  0.11185   
#TrainingLumAB$node1   0.73065   2.07643  0.34248 2.133  0.03289 * 
#TrainingLumAB$grade2  1.05105   2.86065  0.64441 1.631  0.10288   
#TrainingLumAB$grade3  1.69051   5.42224  0.64444 2.623  0.00871 **
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#exp(coef) exp(-coef) lower .95 upper .95
#TrainingLumAB$d_score     1.035     0.9657    0.9331     1.149
#TrainingLumAB$size1       1.975     0.5062    0.8535     4.572
#TrainingLumAB$node1       2.076     0.4816    1.0612     4.063
#TrainingLumAB$grade2      2.861     0.3496    0.8090    10.115
#TrainingLumAB$grade3      5.422     0.1844    1.5333    19.175

#unadjusted survival association in training dataset with disparity score (in category)
f.lum.surv.death <- coxph(Surv(TrainingLumAB$S_Time,TrainingLumAB$S_Event) ~ TrainingLumAB$f.dscore.c)
summary(f.lum.surv.death)
#                                coef exp(coef) se(coef)     z Pr(>|z|)
#TrainingLumAB$f.dscore.c[-2,3) 0.3783    1.4598   0.4063 0.931    0.352
#TrainingLumAB$f.dscore.c[3,7)  0.8088    2.2453   0.5005 1.616    0.106



adj.all.surv.death <- coxph(Surv(TrainingLumAB$S_Time,TrainingLumAB$S_Event) 
                        ~ TrainingLumAB$d_score + TrainingLumAB$size +
                          TrainingLumAB$node + TrainingLumAB$grade)

summary(adj.all.surv.death)




#######the association between race and proliferation score, disparity points
#this is the discrete variable
t.test(pData(race.combinedtraining.eset)$d_score~factor(pData(race.combinedtraining.eset)$Race))
#t = 7.969, df = 85.194, p-value = 6.451e-12
#95 percent confidence interval:
# 2.046963 3.407911 
#       2.4150943       -0.3123426 


############################################################################
#MD 11/2014: find node, size, grade in test dataset.  look at association between
# disparity points and survival after adjusting for these clinical features
# NOW DO ASSOCIATIONS IN TEST DATASET: METBARIC (CALDAS)
##############################################################################
####################################
# read in dataset
load(file="/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/caldas_race_eset.RData")

dim(caldas.race.eset) # this contains the proliferation genes
#Features  Samples 
#       14     1584 

##################################################################
# quick check of the tumor characteristics available
##################################################################
colnames(pData(caldas.race.eset))

table(caldas.race.eset$PAM50)
# Basal   Her2   LumA   LumB Normal 
#   339    223    401    397    224 
table(caldas.race.eset$grade)
#  1   2   3 
# 126 618 840 
table(caldas.race.eset$Node)
#  0   1   2 
# 708 585 291 

table(caldas.race.eset$size)
# these are strange values - need to look up metric
hist(as.numeric(caldas.race.eset$size))
mean(as.numeric(as.character(caldas.race.eset$size)))
#[1] 26.62723

caldas.race.eset$size_group<-ifelse(as.numeric(as.character(caldas.race.eset$size))>=20,1,0)
table(caldas.race.eset$size_group)

mean(as.numeric(as.character(caldas.race.eset$age_at_diagnosis)))
caldas.race.eset$age<-as.numeric(as.character(caldas.race.eset$age_at_diagnosis))
#[1] 60.90307

# create a few more variables for age and size
age_cat_c<-cut(caldas.race.eset$age,breaks=c(0,30,40,50,60,70,80,100),right=FALSE) #  
f.age.c<-factor(age_cat_c)
caldas.scale.eset$f.age.c<-f.age.c
table(f.age.c)
#f.age.c
#  [0,30)  [30,40)  [40,50)  [50,60)  [60,70)  [70,80) [80,100) 
#      11       89      239      370      452      334       89 


##################################################################
# now set direction and create 'points'
##################################################################
# now create proliferation score - sum up rows 5,7-14
proliferation.idx<-which(rownames(exprs(caldas.race.eset))%in%proliferationGenes)
proliferation.idx
prolif_score<-apply(exprs(caldas.race.eset)[proliferation.idx,],2,function(x){sum(x,na.rm = TRUE)})
pData(caldas.race.eset)$prolif_score<-prolif_score

colnames(pData(caldas.race.eset))
rownames(exprs(caldas.race.eset))


medians.caldas <- esApply(caldas.race.eset,1,median)
caldas.exprs.scale <- t(scale(t(exprs(caldas.race.eset)),center=medians.caldas,scale=FALSE))

## Create a new eset
caldas.scale.eset <-caldas.race.eset
exprs(caldas.scale.eset) <- caldas.exprs.scale
dim(caldas.scale.eset)
featureNames(caldas.race.eset)
#[1] "CRYBB2" "MUC1"   "PSPH"   "SQLE"   "TYMS"   "ACOX2"  "BIRC5"  "CCNB1" 
# [9] "CDC20"  "MKI67"  "RRM2"   "PTTG1"  "UBE2C"  "CEP55"
# 1,-1,1,1,1,-1
#create a "ACOX2" "MUC1"  "SQLE"  "TYMS"  variable based on median centered data
# ACOX2 - up is good -1* exprs
# CYRYBB2 - up is bad 1*exprs
# MUC1 - up is good-1*exprs
# SQLE - up is bad 1*exprs
# TYMS - up is bad 1*exprs
#PSPH - up is bad 1*expr



rownames(exprs(caldas.race.eset))
#[1] "CRYBB2" "MUC1"   "PSPH"   "SQLE"   "TYMS"   "ACOX2"  "BIRC5"  "CCNB1" 
# [9] "CDC20"  "MKI67"  "RRM2"   "PTTG1"  "UBE2C"  "CEP55" 
CRYBB2<-sapply(exprs(caldas.scale.eset)[1,],castDirection)
MUC1<-sapply(exprs(caldas.scale.eset)[2,],castDirection)
PSPH<-sapply(exprs(caldas.scale.eset)[3,],castDirection)
SQLE<-sapply(exprs(caldas.scale.eset)[4,],castDirection)
TYMS<-sapply(exprs(caldas.scale.eset)[5,],castDirection)
ACOX2<-sapply(exprs(caldas.scale.eset)[6,],castDirection)

table(CRYBB2)
table(MUC1)
table(PSPH)
table(SQLE)
table(TYMS)
table(ACOX2)

caldas.scale.eset$ACOX2<-ACOX2
caldas.scale.eset$CRYBB2<-CRYBB2
caldas.scale.eset$MUC1<-MUC1
caldas.scale.eset$SQLE<-SQLE
caldas.scale.eset$PSPH<-PSPH
caldas.scale.eset$TYMS<-TYMS

##################scale the results - generate the 'points'
library(psych)
table(caldas.scale.eset$ER_IHC_status)
# neg  pos 
#364 1220 
colnames(pData(caldas.scale.eset))
#[43] "f.age.c"                    "prolif_score"              
#[45] "ACOX2"                      "CRYBB2"                    
#[47] "MUC1"                       "SQLE"                      
#[49] "PSPH"                       "TYMS" 
scale_gene <- matrix(c(-1,1,-1,1,1,1), ncol = 6)
d_score<-scale_gene%*%t(pData(caldas.scale.eset)[,45:50])

quantile(d_score)

pData(caldas.scale.eset)$d_score<-t(d_score)
class(pData(caldas.scale.eset)$d_score)

##################### cdiscrete score
# first evaluate the association between 'points' and proliferation score.
table(pData(caldas.scale.eset)$d_score)
cor(pData(caldas.scale.eset)$prolif_score,pData(caldas.scale.eset)$d_score,use="complete.obs")
#          [,1]
#[1,] 0.5932023
points_fit<-lm(pData(caldas.scale.eset)$prolif_score~pData(caldas.scale.eset)$d_score) 

summary(points_fit)
#pData(caldas.scale.eset)$d_score  2.44783    0.08352  29.308   <2e-16 ***

########does the mean score vary by size_group
t.test(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$size_group)
#t = -2.4352, df = 881.145, p-value = 0.01508
# -0.78953674 -0.08482783
#mean in group 0 mean in group 1 
#     -0.2993763       0.1378060 

describeBy(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$size_group)
#INDICES: 0
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 481 -0.3 3.33      0   -0.35 2.97  -6   6    12 0.15     -0.8 0.15
#------------------------------------------------------------------ 
#INDICES: 1
#   var    n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 1103 0.14 3.19      0    0.12 2.97  -6   6    12 0.03    -0.86 0.1

#####################################################
#does the mean score vary by grade
aov(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$grade)
fit_full<- aov(d_score~grade,data=pData(caldas.scale.eset))


summary(aov(d_score~grade,data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -2.3810     0.2658  -8.958  < 2e-16 ***
#factor(grade)2   1.3162     0.2916   4.513 6.86e-06 ***
#factor(grade)3   3.5310     0.2850  12.387  < 2e-16 ***
#F-statistic: 141.9 on 2 and 1581 DF,  p-value: < 2.2e-16

#does the mean score vary by node status
aov(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$Node)
fit_full<- aov(d_score~factor(Node),data=pData(caldas.scale.eset))


summary(aov(d_score~factor(Node),data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -0.2373     0.1215  -1.953   0.0510 .
#Node1         0.4048     0.1806   2.241   0.0252 *
#Node2         0.5053     0.2251   2.245   0.0249 *

#F-statistic: 3.691 on 2 and 1581 DF,  p-value: 0.02516

dscore_cat_c<-cut(caldas.scale.eset$d_score,breaks=c(-7,-2,3,7),right=FALSE) # corresponds to quantiles
f.dscore.c<-factor(dscore_cat_c)

pData(caldas.scale.eset)$f.dscore.c<-f.dscore.c
table(pData(caldas.scale.eset)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#    309     941     334 


lumA<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("LumA"))]

lumAB<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("LumA","LumB"))]

table(pData(lumA)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#    177     210      14 

table(pData(lumAB)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#    227     485      86 
##########################################################################
#Survival analyses
######################################################## 
#Discrete (points) continuous (cox proportional) model
#Chris Fan told me in an email that these two variables are the survival variables
#RFSE10, #RFS10yr

# points in continuous form
colnames(pData(caldas.scale.eset))

caldas.scale.eset$RFS10yr<-as.numeric(as.character(caldas.scale.eset$RFS10yr))
caldas.scale.eset$RFSE10<-as.numeric(as.character(caldas.scale.eset$RFSE10))
all.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score)
summary(all.surv.death)
#pData(caldas.scale.eset)$d_score 0.0589      1.06   0.0127 4.63 3.7e-06
#pData(caldas.scale.eset)$d_score     1.061     0.9428     1.035     1.087

#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.0589) = [1] 2.027493
##[1] 2.078595
L95 = exp(12*0.0589-1.96*0.0127)
#[1] 1.977647
U95= exp(12*0.0589+1.96*0.0127)
#[1] 2.078595
#HR 0 points/-6 points = exp(6*0.0589) = [1] 1.423901

L95 = exp(6*0.0589-1.96*0.0127)
#[1] 1.388894
U95= exp(6*0.0589+1.96*0.0127)
#[1] 1.459789

#now in categories
all.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c)
summary(all.surv.death)
#caldas.scale.eset$f.dscore.c[-2,3) 0.5054    1.6577   0.1220 4.143 3.43e-05 ***
#caldas.scale.eset$f.dscore.c[3,7)  0.6286    1.8749   0.1392 4.516 6.29e-06 ***
#adjusted analyses over all tumors (adjusting for size, grade, Noe)
#adjust for grade, node, subtype, age, proliferation score, size
#adjust for disparity score and grade
table(caldas.race.eset$PAM50)
table(caldas.race.eset$grade)
table(caldas.race.eset$Node)
table(caldas.race.eset$size_group)
mean(caldas.race.eset$prolif_score)
hist(caldas.race.eset$prolif_score)


#now adjust for grade
grade.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                                caldas.scale.eset$grade)

summary(grade.adj.surv.death)
#in categories
grade.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                                caldas.scale.eset$grade)

summary(grade.adj.surv.death)

#node isn't significant in the model
node.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node)

summary(node.adj.surv.death)
#in categories
node.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node)

summary(node.adj.surv.death)


#add in size to the model
size.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node+
                               caldas.scale.eset$size_group)

summary(size.adj.surv.death)
#in categories
size.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node+
                               caldas.scale.eset$size_group)

summary(size.adj.surv.death)
#add in age to the model

age.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                              caldas.scale.eset$grade+
                              caldas.scale.eset$Node+
                              caldas.scale.eset$size_group+
                              caldas.scale.eset$age)

summary(age.adj.surv.death)
#in categories
age.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                              caldas.scale.eset$grade+
                              caldas.scale.eset$Node+
                              caldas.scale.eset$size_group+
                              caldas.scale.eset$age)

summary(age.adj.surv.death)

#add in proliferation score to the model

ps.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                             caldas.scale.eset$grade+
                             caldas.scale.eset$Node+
                             caldas.scale.eset$size_group+
                             caldas.scale.eset$age+
                             caldas.scale.eset$prolif_score)

summary(ps.adj.surv.death)
#in categories
ps.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                             caldas.scale.eset$grade+
                             caldas.scale.eset$Node+
                             caldas.scale.eset$size_group+
                             caldas.scale.eset$age+
                             caldas.scale.eset$prolif_score)

summary(ps.adj.surv.death)
#with subtype
full.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$d_score+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node+
                               caldas.scale.eset$size_group+
                               caldas.scale.eset$age+
                               caldas.scale.eset$prolif_score+
                               caldas.scale.eset$PAM50)

summary(full.adj.surv.death)
#in categories
full.adj.surv.death <- coxph(Surv(caldas.scale.eset$RFS10yr,caldas.scale.eset$RFSE10) ~ caldas.scale.eset$f.dscore.c+
                               caldas.scale.eset$grade+
                               caldas.scale.eset$Node+
                               caldas.scale.eset$size_group+
                               caldas.scale.eset$age+
                               caldas.scale.eset$prolif_score+
                               caldas.scale.eset$PAM50)

summary(full.adj.surv.death)
###################################################################################
#
# NOW LOOK AT younger women since these ladies are so old (relative to training data)
#
###################################################################################
metbaric_young<-caldas.scale.eset[,which(caldas.scale.eset$age<60)]
dim(metbaric_young)
table(pData(metbaric_young)$f.dscore.c)
colnames(pData(metbaric_young))
metbaric_old<-caldas.scale.eset[,which(caldas.scale.eset$age>=60)]
dim(metbaric_old)
all.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score)
summary(all.surv.death)
#metbaric_young$d_score 0.10710   1.11305  0.02176 4.921 8.59e-07 ***
#metbaric_young$d_score     1.113     0.8984     1.067     1.162

#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.10710) = [1] 2.027493
# 3.615391
L95 = exp(12*0.10710-1.96*0.02176)
L95
#[1] 1.977647
U95= exp(12*0.10710+1.96*0.02176)
U95
#[1] 3.772921
#HR 0 points/-6 points = exp(6*0.10710) = [1] 1.901418
exp(6*0.10710)
L95 = exp(6*0.10710-1.96*0.02176)
#[1] 1.822028
U95= exp(6*0.10710+1.96*0.02176)
#[1] 1.984267

#discrete score
d.all.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c)
summary(d.all.surv.death)
#metbaric_young$f.dscore.c[-2,3) 1.1764    3.2427   0.3023 3.891 9.97e-05 ***
#metbaric_young$f.dscore.c[3,7)  1.4267    4.1649   0.3139 4.545 5.50e-06 ***

#metbaric_young$f.dscore.c[-2,3)     3.243     0.3084     1.793     5.864
#metbaric_young$f.dscore.c[3,7)      4.165     0.2401     2.251     7.706
#now adjust for grade
grade.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                                metbaric_young$grade)

summary(grade.adj.surv.death)
#metbaric_young$d_score 0.05243   1.05383  0.02320 2.260  0.02385 * 
#metbaric_young$d_score     1.054     0.9489    1.0070     1.103

#discrete + grade
d.grade.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                                metbaric_young$grade)

summary(d.grade.adj.surv.death)
#etbaric_young$f.dscore.c[-2,3)     2.651     0.3772    1.4613     4.811
#metbaric_young$f.dscore.c[3,7)      2.527     0.3957    1.3446     4.750

# grade 2 NS  
node.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                               metbaric_young$grade+
                               metbaric_young$Node)

summary(node.adj.surv.death)
#metbaric_young$d_score 0.04571   1.04677  0.02391 1.912  0.05588 . 
#metbaric_young$d_score     1.047     0.9553    0.9989     1.097

#discrete + grade + Node
d.node.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                                  metbaric_young$grade+
                                  metbaric_young$Node)

summary(d.node.adj.surv.death)
#metbaric_young$f.dscore.c[-2,3)     2.237     0.4471    1.2268     4.077
#metbaric_young$f.dscore.c[3,7)      2.126     0.4704    1.1235     4.023

#add in size to the model
size.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                               metbaric_young$grade+
                               metbaric_young$Node+
                               metbaric_young$size_group)

summary(size.adj.surv.death)
#metbaric_young$d_score    0.04483   1.04585  0.02398 1.870  0.06151 .  
#metbaric_young$d_score        1.046     0.9562    0.9978     1.096
# size and grade2 NS
#discrete + grade + Node + size
d.size.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                                 metbaric_young$grade+
                                 metbaric_young$Node+
                                 metbaric_young$size_group)

summary(d.size.adj.surv.death)
#metbaric_young$f.dscore.c[-2,3)     2.216     0.4513    1.2145     4.043
#metbaric_young$f.dscore.c[3,7)      2.103     0.4754    1.1105     3.984

#add in age to the model

age.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                              metbaric_young$grade+
                              metbaric_young$Node+
                              metbaric_young$size_group+
                              metbaric_young$age)

summary(age.adj.surv.death)
#metbaric_young$d_score       1.0441     0.9577    0.9958     1.095
#metbaric_young$d_score     0.043172  1.044117  0.024154  1.787  0.07388 . 

#discrete + grade + Node + size + age
d.age.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                                 metbaric_young$grade+
                                 metbaric_young$Node+
                                 metbaric_young$size_group+
                                metbaric_young$age)

summary(d.age.adj.surv.death)
#metbaric_young$f.dscore.c[-2,3)    2.2004     0.4545    1.2055     4.016
#metbaric_young$f.dscore.c[3,7)     2.0614     0.4851    1.0857     3.914


#add in proliferation score to the model

ps.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                             metbaric_young$grade+
                             metbaric_young$Node+
                             metbaric_young$size_group+
                             metbaric_young$age+
                             metbaric_young$prolif_score)

summary(ps.adj.surv.death)
#metbaric_young$d_score         0.9910     1.0091    0.9376     1.047 
#metbaric_young$prolif_score    1.0267     0.9740    1.0116     1.042
#metbaric_young$age             0.9976     1.0024    0.9804     1.015
#metbaric_young$d_score      -0.009039  0.991001  0.028246 -0.320 0.748950  
#metbaric_young$prolif_score  0.026350  1.026700  0.007567  3.482 0.000497 ***


#discrete + grade + Node + size + age + proliferation score
d.ps.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                                metbaric_young$grade+
                                metbaric_young$Node+
                                metbaric_young$size_group+
                                metbaric_young$age+
                                metbaric_young$prolif_score)

summary(d.ps.adj.surv.death)
#metbaric_young$f.dscore.c[-2,3)    1.7903     0.5586    0.9710     3.301
#metbaric_young$f.dscore.c[3,7)     1.2758     0.7838    0.6427     2.533
#metbaric_young$grade2              1.6613     0.6020    0.5894     4.682
#metbaric_young$grade3              2.5883     0.3864    0.9107     7.356
#metbaric_young$Node1               1.3214     0.7568    0.9281     1.881
#metbaric_young$Node2               3.1540     0.3171    2.1772     4.569
#metbaric_young$size_group          1.0997     0.9093    0.7895     1.532
#metbaric_young$age                 0.9968     1.0032    0.9797     1.014
#metbaric_young$prolif_score        1.0281     0.9726    1.0135     1.043

full.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$d_score+
                               metbaric_young$grade+
                               metbaric_young$Node+
                               metbaric_young$size_group+
                               metbaric_young$age+
                               metbaric_young$prolif_score+
                               metbaric_young$PAM50)

summary(full.adj.surv.death)
#                                 coef exp(coef)  se(coef)      z Pr(>|z|)    
#metbaric_young$d_score      -0.018592  0.981580  0.029911 -0.622  0.53422    
#metbaric_young$grade2        0.368678  1.445821  0.533610  0.691  0.48962    
#metbaric_young$grade3        0.714899  2.043980  0.539063  1.326  0.18478    
#metbaric_young$Node1         0.228662  1.256917  0.181673  1.259  0.20816    
#metbaric_young$Node2         1.158669  3.185689  0.189575  6.112 9.84e-10 ***
#metbaric_young$size_group    0.107053  1.112994  0.168954  0.634  0.52633    
#metbaric_young$age          -0.003282  0.996723  0.009037 -0.363  0.71647    
#metbaric_young$prolif_score  0.012317  1.012393  0.009308  1.323  0.18574    
#metbaric_young$PAM50Her2     0.169415  1.184612  0.198004  0.856  0.39221    
#metbaric_young$PAM50LumA    -1.265954  0.281970  0.407498 -3.107  0.00189 ** 
#metbaric_young$PAM50LumB    -0.147012  0.863284  0.234010 -0.628  0.52985    
#metbaric_young$PAM50Normal  -0.435561  0.646902  0.338425 -1.287  0.19809    

#discrete + grade + Node + size + age + proliferation score +subtype
d.subtype.adj.surv.death <- coxph(Surv(metbaric_young$RFS10yr,metbaric_young$RFSE10) ~ metbaric_young$f.dscore.c+
                               metbaric_young$grade+
                               metbaric_young$Node+
                               metbaric_young$size_group+
                               metbaric_young$age+
                               metbaric_young$prolif_score+
                               metbaric_young$PAM50)

summary(d.subtype.adj.surv.death)
#                                 coef exp(coef)  se(coef)      z Pr(>|z|)    
#metbaric_young$f.dscore.c[-2,3)  0.453948  1.574517  0.315235  1.440  0.14986    
#metbaric_young$f.dscore.c[3,7)   0.144875  1.155895  0.357088  0.406  0.68495
#metbaric_young$grade2            0.320988  1.378489  0.533508  0.602  0.54740  
#metbaric_young$grade3            0.675692  1.965392  0.538277  1.255  0.20937 
#metbaric_young$Node1             0.229129  1.257504  0.181502  1.262  0.20680   
#metbaric_young$Node2             1.140847  3.129418  0.189427  6.023 1.72e-09 ***
#metbaric_young$size_group        0.092435  1.096842  0.169243  0.546  0.58495   
#metbaric_young$age              -0.003057  0.996948  0.009044 -0.338  0.73535   
#metbaric_young$prolif_score      0.013829  1.013925  0.009266  1.492  0.13560  
#metbaric_young$PAM50Her2         0.114545  1.121363  0.198108  0.578  0.56313  
#metbaric_young$PAM50LumA        -1.201637  0.300701  0.399131 -3.011  0.00261 ** 
#metbaric_young$PAM50LumB        -0.178646  0.836402  0.229713 -0.778  0.43675     
#metbaric_young$PAM50Normal      -0.380455  0.683550  0.331863 -1.146  0.25162 

###############################################################################
#METBARIC OLD
###############################################################################

all.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score)
summary(all.surv.death)
#metbaric_old$d_score 0.04583   1.04690  0.01612 2.843  0.00447 **
#metbaric_old$d_score     1.047     0.9552     1.014      1.08

#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.04583) = [1] 2.027493
# 3.615391
L95 = exp(12*0.04583-1.96*0.01612)
L95
#[1] 1.679279
U95= exp(12*0.04583+1.96*0.01612)
U95
#[1] 1.788818


#discrete score
d.all.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c)
summary(d.all.surv.death)
#metbaric_old$f.dscore.c[-2,3) 0.3812    1.4640   0.1359 2.805  0.00503 **
#metbaric_old$f.dscore.c[3,7)  0.4680    1.5968   0.1680 2.786  0.00533 **

#metbaric_old$f.dscore.c[-2,3)     1.464     0.6831     1.122     1.911
#metbaric_old$f.dscore.c[3,7)      1.597     0.6263     1.149     2.219
#now adjust for grade
grade.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                                metbaric_old$grade)

summary(grade.adj.surv.death)
#metbaric_young$d_score 0.05243   1.05383  0.02320 2.260  0.02385 * 
#metbaric_young$d_score     1.054     0.9489    1.0070     1.103

#discrete + grade
d.grade.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                                  metbaric_old$grade)

summary(d.grade.adj.surv.death)
#metbaric_old$f.dscore.c[-2,3)     1.371     0.7294    1.0427     1.803
#metbaric_old$f.dscore.c[3,7)      1.442     0.6936    1.0231     2.031

# grade 2 NS  
node.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                               metbaric_old$grade+
                               metbaric_old$Node)

summary(node.adj.surv.death)
#metbaric_old$d_score 0.04552   1.04657  0.01709 2.664  0.00773 ** 
#metbaric_old$d_score     1.047     0.9555    1.0121     1.082

#discrete + grade + Node
d.node.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                                 metbaric_old$grade+
                                 metbaric_old$Node)

summary(d.node.adj.surv.death)
#metbaric_old$f.dscore.c[-2,3)     1.445     0.6920    1.0974     1.903
#metbaric_old$f.dscore.c[3,7)      1.612     0.6202    1.1437     2.273

#add in size to the model
size.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                               metbaric_old$grade+
                               metbaric_old$Node+
                               metbaric_old$size_group)

summary(size.adj.surv.death)
#metbaric_old$d_score    0.04372   1.04469  0.01720 2.542 0.011007 * 
#metbaric_old$d_score        1.045     0.9572    1.0101     1.080
# size and grade2 NS
#discrete + grade + Node + size
d.size.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                                 metbaric_old$grade+
                                 metbaric_old$Node+
                                 metbaric_old$size_group)

summary(d.size.adj.surv.death)
#metbaric_old$f.dscore.c[-2,3)     1.422     0.7032    1.0804     1.872
#metbaric_old$f.dscore.c[3,7)      1.598     0.6259    1.1338     2.251

#add in age to the model

age.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                              metbaric_old$grade+
                              metbaric_old$Node+
                              metbaric_old$size_group+
                              metbaric_old$age)

summary(age.adj.surv.death)
#metbaric_old$d_score        1.051     0.9516    1.0159     1.087
#metbaric_old$d_score    0.049560  1.050808 0.017224 2.877 0.004011 **  

#discrete + grade + Node + size + age
d.age.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                                metbaric_old$grade+
                                metbaric_old$Node+
                                metbaric_old$size_group+
                                metbaric_old$age)

summary(d.age.adj.surv.death)
#metbaric_old$f.dscore.c[-2,3)     1.431     0.6988    1.0866     1.885
#metbaric_old$f.dscore.c[3,7)      1.613     0.6199    1.1436     2.276


#add in proliferation score to the model

ps.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                             metbaric_old$grade+
                             metbaric_old$Node+
                             metbaric_old$size_group+
                             metbaric_old$age+
                             metbaric_old$prolif_score)

summary(ps.adj.surv.death)
#metbaric_old$d_score          1.014     0.9865    0.9746     1.054
#metbaric_old$grade2           1.119     0.8937    0.6969     1.797
#metbaric_old$grade3           1.223     0.8174    0.7520     1.990
#metbaric_old$Node1            1.562     0.6403    1.2250     1.991
#metbaric_old$Node2            2.749     0.3638    2.1051     3.590
#metbaric_old$size_group       1.379     0.7251    1.0575     1.799
#metbaric_old$age              1.056     0.9467    1.0406     1.072
#metbaric_old$prolif_score     1.020     0.9803    1.0084     1.032

#discrete + grade + Node + size + age + proliferation score
d.ps.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                               metbaric_old$grade+
                               metbaric_old$Node+
                               metbaric_old$size_group+
                               metbaric_old$age+
                               metbaric_old$prolif_score)

summary(d.ps.adj.surv.death)
#metbaric_old$f.dscore.c[-2,3)     1.239     0.8070    0.9297     1.652
#metbaric_old$f.dscore.c[3,7)      1.172     0.8534    0.7973     1.722
#metbaric_old$grade2               1.087     0.9202    0.6757     1.748
#metbaric_old$grade3               1.177     0.8496    0.7217     1.920
#metbaric_old$Node1                1.553     0.6441    1.2177     1.980
#metbaric_old$Node2                2.734     0.3658    2.0925     3.571
#metbaric_old$size_group           1.376     0.7266    1.0550     1.795
#metbaric_old$age                  1.056     0.9469    1.0404     1.072
#metbaric_old$prolif_score         1.020     0.9802    1.0090     1.032


full.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$d_score+
                               metbaric_old$grade+
                               metbaric_old$Node+
                               metbaric_old$size_group+
                               metbaric_old$age+
                               metbaric_old$prolif_score+
                               metbaric_old$PAM50)

summary(full.adj.surv.death)
#                                 coef exp(coef)  se(coef)      z Pr(>|z|)    
#metbaric_old$d_score      -0.001047  0.998954  0.021008 -0.050  0.96025    
#metbaric_old$grade2        0.129421  1.138169  0.242468  0.534  0.59350    
#metbaric_old$grade3        0.144353  1.155291  0.251480  0.574  0.56596    
#metbaric_old$Node1         0.498994  1.647064  0.125314  3.982 6.83e-05 ***
#  metbaric_old$Node2         1.086189  2.962961  0.138109  7.865 3.66e-15 ***
#  metbaric_old$size_group    0.346804  1.414539  0.136189  2.546  0.01088 *  
#  metbaric_old$age           0.057367  1.059044  0.007681  7.468 8.13e-14 ***
#  metbaric_old$prolif_score  0.019281  1.019468  0.007371  2.616  0.00890 ** 
#  metbaric_old$PAM50Her2    -0.266155  0.766320  0.189889 -1.402  0.16103    
#metbaric_old$PAM50LumA    -0.526213  0.590838  0.224472 -2.344  0.01907 *  
#  metbaric_old$PAM50LumB    -0.608491  0.544172  0.166592 -3.653  0.00026 ***
#  metbaric_old$PAM50Normal  -0.215247  0.806343  0.269591 -0.798  0.42463   

#discrete + grade + Node + size + age + proliferation score +subtype
d.subtype.adj.surv.death <- coxph(Surv(metbaric_old$RFS10yr,metbaric_old$RFSE10) ~ metbaric_old$f.dscore.c+
                                    metbaric_old$grade+
                                    metbaric_old$Node+
                                    metbaric_old$size_group+
                                    metbaric_old$age+
                                    metbaric_old$prolif_score+
                                    metbaric_old$PAM50)

summary(d.subtype.adj.surv.death)
#                                 coef exp(coef)  se(coef)      z Pr(>|z|)    
#metbaric_old$f.dscore.c[-2,3)  0.189647  1.208823  0.148967  1.273 0.202990    
#metbaric_old$f.dscore.c[3,7)  -0.004365  0.995645  0.206068 -0.021 0.983100    
#metbaric_old$grade2            0.096169  1.100945  0.243276  0.395 0.692616    
#metbaric_old$grade3            0.096082  1.100850  0.252771  0.380 0.703859    
#metbaric_old$Node1             0.494046  1.638935  0.125271  3.944 8.02e-05 ***
#metbaric_old$Node2             1.080220  2.945326  0.138323  7.809 5.77e-15 ***
#metbaric_old$size_group        0.349394  1.418208  0.136323  2.563 0.010378 *  
#metbaric_old$age               0.057614  1.059306  0.007707  7.476 7.68e-14 ***
#metbaric_old$prolif_score      0.019966  1.020167  0.007276  2.744 0.006066 ** 
#metbaric_old$PAM50Her2        -0.301175  0.739949  0.190419 -1.582 0.113732    
#metbaric_old$PAM50LumA        -0.534553  0.585931  0.221438 -2.414 0.015778 *  
#metbaric_old$PAM50LumB        -0.653899  0.520014  0.168636 -3.878 0.000105 ***
#metbaric_old$PAM50Normal      -0.255306  0.774680  0.268480 -0.951 0.341640    





##########################points in categorical f.dscore.c

######################################################################
# Luminal AB associations
######################################################################
# unadjusted, points from -6 -> +6

#just luminal A
lumA.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ 
                            pData(lumA)$d_score)
summary(lumA.surv.death)
#                        coef exp(coef) se(coef)     z Pr(>|z|)  
#pData(lumA)$d_score 0.04708   1.04820  0.03658 1.287    0.198

#                     exp(coef) exp(-coef) lower .95 upper .95
#pData(lumA)$d_score     1.048      0.954    0.9757     1.126
#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.04708)
##[1] 1.759377
L95 = exp(12*0.04708-1.96*0.03658)
L95
#[1] 1.637651
U95= exp(12*0.04708+1.96*0.03658)
U95
#[1] 1.890151

#luminal A adjusted
adj.lumA.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ 
                           pData(lumA)$d_score+ pData(lumA)$grade + 
                           pData(lumA)$size_group + pData(lumA)$Node)
summary(adj.lumA.surv.death)
#                          coef exp(coef) se(coef)     z Pr(>|z|)    
#pData(lumA)$d_score    0.02537   1.02569  0.03803 0.667 0.504676    
#pData(lumA)$grade2     0.08901   1.09310  0.29221 0.305 0.760656    
#pData(lumA)$grade3     0.17642   1.19293  0.32160 0.549 0.583315    
#pData(lumA)$size_group 0.80891   2.24546  0.24459 3.307 0.000942 ***
#  pData(lumA)$Node1      0.30138   1.35172  0.22508 1.339 0.180570    
#pData(lumA)$Node2      0.80355   2.23345  0.24714 3.251 0.001148 ** 

#                       exp(coef) exp(-coef) lower .95 upper .95
#pData(lumA)$d_score        1.026     0.9750    0.9520     1.105

#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.02537)
##[1] 1.355866
L95 = exp(12*0.02537-1.96*0.03803)
L95
#[1] 1.258476
U95= exp(12*0.02537+1.96*0.03803)
U95
#[1] 1.460792

#luminal A abd B

lumAB.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ 
                                pData(lumAB)$d_score)
summary(lumAB.surv.death)
#                        coef exp(coef) se(coef)     z Pr(>|z|)  
#pData(lumAB)$d_score 0.03908   1.03985  0.02001 1.953   0.0508 .

#                     exp(coef) exp(-coef) lower .95 upper .95
#pData(lumAB)$d_score      1.04     0.9617    0.9999     1.081

adj.lumAB.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ 
                            pData(lumAB)$d_score + pData(lumAB)$grade + 
                            pData(lumAB)$size_group + pData(lumAB)$Node)

summary(adj.lumAB.surv.death)
#coef exp(coef) se(coef)     z Pr(>|z|)   
#pData(lumAB)$d_score    0.02442   1.02472  0.02106 1.160 0.246096    
#pData(lumAB)$grade2     0.10519   1.11093  0.23774 0.442 0.658141    
#pData(lumAB)$grade3     0.24901   1.28275  0.24389 1.021 0.307261    
#pData(lumAB)$size_group 0.60271   1.82707  0.15936 3.782 0.000155 ***
#  pData(lumAB)$Node1      0.57446   1.77617  0.14521 3.956 7.62e-05 ***
#  pData(lumAB)$Node2      1.01185   2.75069  0.16398 6.171 6.80e-10 ***

#                                    exp(coef) exp(-coef) lower .95 upper .95
#pData(lumAB)$d_score        1.025     0.9759    0.9833     1.068
#pData(lumAB)$grade2         1.111     0.9002    0.6971     1.770
#pData(lumAB)$grade3         1.283     0.7796    0.7953     2.069
#pData(lumAB)$size_group     1.827     0.5473    1.3369     2.497
#pData(lumAB)$Node1          1.776     0.5630    1.3362     2.361
#pData(lumAB)$Node2          2.751     0.3635    1.9946     3.793
#HR 6 points/-6 points = exp(12*0.0589) = [1] 2.027493
##[1] 2.078595
L95 = exp(12*0.0589-1.96*0.0127)
#[1] 1.977647
U95= exp(12*0.0589+1.96*0.0127)
##########################points in categorical f.dscore.c
lumAB.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ 
                            pData(lumAB)$f.dscore.c)
summary(lumAB.surv.death)
#                               coef exp(coef) se(coef)     z Pr(>|z|)   
#pData(lumAB)$f.dscore.c[-2,3) 0.4657    1.5932   0.1485 3.135  0.00172 **
#  pData(lumAB)$f.dscore.c[3,7)  0.4061    1.5009   0.2257 1.799  0.07197 . 

#                        coef exp(coef) se(coef)     z Pr(>|z|)  
#pData(lumAB)$d_score 0.03908   1.03985  0.02001 1.953   0.0508 .

#                     exp(coef) exp(-coef) lower .95 upper .95
#pData(lumAB)$f.dscore.c[-2,3)     1.593     0.6277    1.1907     2.132
#pData(lumAB)$f.dscore.c[3,7)      1.501     0.6663    0.9644     2.336

#adjusted luminal A/B associations
adj.lumAB.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ 
                            pData(lumAB)$f.dscore.c+ pData(lumAB)$grade + 
                            pData(lumAB)$size_group + pData(lumAB)$Node)

#                                 coef exp(coef) se(coef)     z Pr(>|z|)    
#pData(lumAB)$f.dscore.c[-2,3) 0.37965   1.46177  0.15222 2.494 0.012631 *  
#  pData(lumAB)$f.dscore.c[3,7)  0.32763   1.38767  0.22970 1.426 0.153769    
#pData(lumAB)$grade2           0.07159   1.07422  0.23779 0.301 0.763358    
#pData(lumAB)$grade3           0.19414   1.21427  0.24327 0.798 0.424836    
#pData(lumAB)$size_group       0.59920   1.82066  0.15943 3.758 0.000171 ***
#  pData(lumAB)$Node1            0.58824   1.80081  0.14535 4.047 5.18e-05 ***
#  pData(lumAB)$Node2            0.99777   2.71222  0.16428 6.074 1.25e-09 ***

#                              exp(coef) exp(-coef) lower .95 upper .95
#pData(lumAB)$f.dscore.c[-2,3)     1.462     0.6841    1.0847     1.970
#pData(lumAB)$f.dscore.c[3,7)      1.388     0.7206    0.8846     2.177
#pData(lumAB)$grade2               1.074     0.9309    0.6740     1.712
#pData(lumAB)$grade3               1.214     0.8235    0.7538     1.956
#pData(lumAB)$size_group           1.821     0.5493    1.3320     2.489
#pData(lumAB)$Node1                1.801     0.5553    1.3544     2.394
#pData(lumAB)$Node2                2.712     0.3687    1.9656     3.742


