# this code will create an aggregate score of race/survival associated genes:
# CYRYBB2 - up is bad 1*exprs HR = 1.4
# ACOX2 - up is good -1* exprs HR = 0.65
# SQLE - up is bad 1*exprs HR = 1.98
# MUC1 - up is good-1*exprs HR = 0.65
# TYMS - up is bad 1*exprs HR = 2.67
# PSPH - up is bad 1*exprs HR = 1.65

# it also evaluates point differences by subtype: see line 270


# evaluate both discrete score and score.
# make sure that the associations hold for all the genes (CRYBB2, ACOX1 etc)

library(Biobase)
library(gdata) #
library(hgug4112a.db)
library(gplots)
library(limma)
library(gage)
library(survival)
library(muhaz)

data(egSymb)



race_symbols<-c("CRYBB2","MUC1", "ACOX2", "SQLE", "TYMS", "PSPH")
#MD - update 1/2014.  calculate progression score 
proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
race_entrez<-sym2eg(race_symbols)

source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/pcaFuncs.R")
source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts.R")

#  set working directory
setwd("/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper")

#############################################################
#MD 1/1/2014: find node, size, grade & age variables in both datasets
# evaluate the association with 'progression points'
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
#Features  Samples 
#       10      165 
rownames(exprs(race.337.eset))<-fData(race.337.eset)$symbol

race.337.eset$race2<-race.337.eset$Race

#################### remove those samples without missing survival data
notmissing.ind<-which(sapply(race.337.eset$Overall.Survival.Event..0.alive..1.DOD.or.DOC..x,function(x) !is.na(x)))
length(notmissing.ind)
race.337.eset<-race.337.eset[,notmissing.ind]

race.337.eset$S_Event<-race.337.eset$Overall.Survival.Event..0.alive..1.DOD.or.DOC..x
race.337.eset$S_Time<-race.337.eset$Overall.suvival.months.x/12

#LumA <- sub.337.eset[,which(sub.337.eset$PAM50.Call %in% c("LumA"))]
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

#create grade variable
table(race.337.eset$Grade)
# 1  2  3 
#14 53 80 

table(race.337.eset$PAM50.Call)
# Basal   Her2   LumA   LumB Normal 
#    35     17     66     26     11 

colnames(pData(race.337.eset))
# median center the data
medians<- esApply(race.337.eset,1,median)
unc337.exprs.scale <- t(scale(t(exprs(race.337.eset)),center=medians,scale=FALSE))


## Create a new eset
unc337.scale.eset<-race.337.eset
exprs(unc337.scale.eset) <- unc337.exprs.scale
dim(unc337.scale.eset)
##10 155

####################################################################################################
# read in NKI data
load("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/nki295_eset.RData") ## process_NKI295samples.R
dim(nki295.eset)
nki295.eset <- nki295.eset[!is.na(featureData(nki295.eset)$geneid),]
dim(nki295.eset)
#Features  Samples 
#   17600      295 
#length(intersect(fData(nki295.eset)$symbol,age$symbols))
race.nki295.eset <- nki295.eset[fData(nki295.eset)$symbol %in% unique(c(race_symbols,proliferationGenes)),]
rownames(exprs(race.nki295.eset))
 #[1] "SQLE"   "MKI67"  "MUC1"   "PTTG1"  "ACOX2"  "PSPH"   "KNTC2"  "UBE2C"  "CCNB1"  "RRM2"   "CEP55"  "RRM2"  
#[13] "TYMS"   "CDCA1"  "BIRC5"  "CRYBB2" "CDC20" 

#there are duplicates
exprs(race.nki295.eset)<-collapseIDs(exprs(race.nki295.eset),fData(race.nki295.eset)$symbol,"mean")
dim(exprs(race.nki295.eset))


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
#table(race.nki295.eset$Size)
  0   1 
114 181 

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
#dim(race.nki295.eset)
#Features  Samples 
#       6      295 

pData(race.nki295.eset)$subtype

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
#create new expression set with combined data, then median center
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
				
age_cat<-rbind(as.matrix(nki295.scale.eset$f.age.c),
				as.matrix(unc337.scale.eset$f.age.c))
				
data_source<-rbind(as.matrix(rep("NKI",295)),as.matrix(rep("UNC337",155)))	

				

combined<-cbind(data_source,STime_combined,SEvent_combined,race_combined,subtype,grade,node,size,age_cat)
rownames(combined)<-row.names
colnames(combined)<-c("Data_Source","S_Time","S_Event","Race","subtype","grade","node","size","age_cat")

patient.phenotype <- data.frame(combined)
hist(as.numeric(as.character(patient.phenotype$S_Time)))
combined.phenoData <- new("AnnotatedDataFrame", data = patient.phenotype)

#combined expression data

two.data<-merge(exprs(nki295.scale.eset),exprs(unc337.scale.eset),by='row.names', all.x=TRUE, all.y=TRUE, incomparables=NA)
rownames(two.data)<-two.data$Row.names
head(colnames(two.data))
two.data<-two.data[,-1] # remove the first column


dim(two.data)

#[1] 6 450


tumor6.race.combined.eset <- new("ExpressionSet", exprs = data.matrix(two.data),
                   phenoData = combined.phenoData, annotation = "hgug4112a")


save(tumor6.race.combined.eset,file="tumor6_race_combined_eset.RData")
####################################################################################################
####################################################################################################
rownames(exprs(tumor6.race.combined.eset))
# [1] "ACOX2"  "BIRC5"  "CCNB1"  "CDC20"  "CDCA1"  "CEP55"  "CRYBB2" "KNTC2"  "MKI67"  "MUC1"   "PSPH"   "PTTG1" 
#[13] "RRM2"   "SQLE"   "TYMS"   "UBE2C" 
proliferation.idx<-which(rownames(exprs(tumor6.race.combined.eset))%in%proliferationGenes)
proliferation.idx
prolif_score<-apply(exprs(tumor6.race.combined.eset)[proliferation.idx,],2,sum)
pData(tumor6.race.combined.eset)$prolif_score<-prolif_score

rownames(exprs(tumor6.race.combined.eset))
#[1][1] "ACOX2"  "CRYBB2" "MUC1"   "PSPH"   "SQLE"   "TYMS" 
#create a "ACOX2" "MUC1"  "SQLE"  "TYMS"  variable based on median centered data
# ACOX2 - up is good -1* exprs
# CYRYBB2 - up is bad 1*exprs
# MUC1 - up is good-1*exprs
# SQLE - up is bad 1*exprs
# TYMS - up is bad 1*exprs
#PSPH - up is bad 1*expr

ACOX2<-sapply(exprs(tumor6.race.combined.eset)[1,],castDirection)
CRYBB2<-sapply(exprs(tumor6.race.combined.eset)[2,],castDirection)
MUC1<-sapply(exprs(tumor6.race.combined.eset)[3,],castDirection)
PSPH<-sapply(exprs(tumor6.race.combined.eset)[4,],castDirection)
SQLE<-sapply(exprs(tumor6.race.combined.eset)[5,],castDirection)
TYMS<-sapply(exprs(tumor6.race.combined.eset)[6,],castDirection)

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

tumor6.race.combined.eset$ACOX2<-ACOX2
tumor6.race.combined.eset$CRYBB2<-CRYBB2
tumor6.race.combined.eset$MUC1<-MUC1
tumor6.race.combined.eset$SQLE<-SQLE
tumor6.race.combined.eset$PSPH<-PSPH
tumor6.race.combined.eset$TYMS<-TYMS

scale_gene <- matrix(c(-1,1,-1,1,1,1), ncol = 6)

colnames(pData(tumor6.race.combined.eset))

# [1] "Data_Source" "S_Time"      "S_Event"     "Race"        "subtype"     "grade"       "node"        "size"       
#  [8] "size"         "age_cat"      "prolif_score" "ACOX2"        "CRYBB2"       "MUC1"         "SQLE"        
#[15] "PSPH"         "TYMS"  
d_score<-scale_gene%*%t(pData(tumor6.race.combined.eset)[,11:16])

quantile(d_score)
#  0%  25%  50%  75% 100% 
#  -6   -2    0    2    6 
#<-2, >3 are good cutpoints

pData(tumor6.race.combined.eset)$d_score<-t(d_score)
dscore_cat_c<-cut(tumor6.race.combined.eset$d_score,breaks=c(-7,-2,3,7),right=FALSE) # corresponds to quantiles
f.dscore.c<-factor(dscore_cat_c)

pData(tumor6.race.combined.eset)$f.dscore.c<-f.dscore.c

table(pData(tumor6.race.combined.eset)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#     92     274      84 


score<-scale_gene%*%exprs(tumor6.race.combined.eset)
quantile(score)
#         0%         25%         50%         75%        100% 
#-8.3890000 -0.8682958 -0.0236250  1.1136555  9.2280000  
pData(tumor6.race.combined.eset)$score<-t(score)

score_cat_c<-cut(tumor6.race.combined.eset$score,breaks=c(-10,-0.09,1.11,10),right=FALSE) # corresponds to quantiles
f.score.c<-factor(score_cat_c)
pData(tumor6.race.combined.eset)$f.score.c<-f.score.c

#this is the continuous variable
t.test(pData(tumor6.race.combined.eset)$score~factor(pData(tumor6.race.combined.eset)$Race))
#t = 8.7407, df = 60.206, p-value = 2.628e-12
#95 percent confidence interval:
# 2.813511 4.483248 
#        3.327660        -0.320719 

#this is the discrete variable
t.test(pData(tumor6.race.combined.eset)$d_score~factor(pData(tumor6.race.combined.eset)$Race))
#t = 7.969, df = 85.194, p-value = 6.451e-12
#95 percent confidence interval:
# 2.046963 3.407911 
#       2.4150943       -0.3123426 

table(pData(tumor6.race.combined.eset)$Race)
#"t=7.969, DF=85.194, p-value=6.451e-12"
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_RaceScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(Race),data=pData(tumor6.race.combined.eset), col=c("grey","white"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


######################## look at average score by subtype and then by race in Luminal A tumors and then in Luminal A/B tumors, and overall
lumA<-tumor6.race.combined.eset[,which(subtype%in%c("LumA"))]
#Features  Samples 
#       6      152 
lumAB<-tumor6.race.combined.eset[,which(subtype%in%c("LumA","LumB"))]
#Features  Samples 
#       6      259


#this is the continuous variable
t.test(pData(lumA)$score~factor(pData(lumA)$Race))
#t = 7.3119, df = 21.92, p-value = 2.598e-07
#95 percent confidence interval:
#  3.058206 5.480548
#        2.383056        -1.886322 

table(lumA$Race)
# AA   C 
# 18 134 

#this is the discrete variable
mean(pData(lumA)$d_score)
#[1] -1.947368
#make a boxplot for AA vs Caucasians
t.test(pData(lumA)$d_score~factor(pData(lumA)$Race))
#t = 7.2599, df = 23.959, p-value = 1.7e-07
#95 percent confidence interval:
# 2.933961 5.265044 
#         1.666667        -2.432836 

#"t=7.26, DF=23.96, p-value=1.7e-07"
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_LumARaceScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(Race),data=pData(lumA), col=c("grey","white"))
abline(h = 0, col = "black", lty="dotted")
dev.off()


#this is the discrete variable for Luminal AB/B
mean(pData(lumAB)$d_score)
#[1] -1.947368
#make a boxplot for AA vs Caucasians
t.test(pData(lumAB)$d_score~factor(pData(lumAB)$Race))
#t = 5.9428, df = 45.123, p-value = 3.762e-07
#95 percent confidence interval:
#  1.727189 3.497919 
#         1.642857        -0.969697 
#"t=5.94, DF=45.12, p-value=3.8e-07"
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_LumABRaceScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(Race),data=pData(lumAB), col=c("grey","white"),ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


# ######################## 
#first evaluate the association between 'points' and proliferation score.
table(pData(tumor6.race.combined.eset)$d_score)
cor(pData(tumor6.race.combined.eset)$prolif_score,pData(tumor6.race.combined.eset)$d_score,use="complete.obs")
#          [,1]
#[1,] 0.5593362
points_fit<-lm(pData(tumor6.race.combined.eset)$prolif_score~pData(tumor6.race.combined.eset)$d_score) 

summary(points_fit)
#pData(tumor6.race.combined.eset)$d_score  0.44603    0.03862  11.550   <2e-16 ***
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_ProliferationPointsCorrelation.tiff',width=350,height=350,res=100,antialias="none")
plot(pData(tumor6.race.combined.eset)$d_score,pData(tumor6.race.combined.eset)$prolif_score)
abline(points_fit)
dev.off()


######################## look at average score by subtype
aov(pData(tumor6.race.combined.eset)$d_score~pData(tumor6.race.combined.eset)$subtype)
fit_full<- aov(d_score~factor(subtype),data=pData(tumor6.race.combined.eset))


summary(aov(d_score~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#does the mean score vary by subtype


table(tumor6.race.combined.eset$subtype)
# Basal   Her2   LumA   LumB Normal 
#    81     66    152    107     44 
summary(aov(d_score~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             2.7654     0.3030   9.126  < 2e-16 ***
#factor(subtype)Her2    -2.5836     0.4522  -5.713 2.03e-08 ***
#factor(subtype)LumA    -4.7128     0.3752 -12.562  < 2e-16 ***
#factor(subtype)LumB    -1.6626     0.4017  -4.139 4.16e-05 ***
#factor(subtype)Normal  -3.9927     0.5107  -7.818 3.93e-14 ***
#F-statistic: 46.87 on 4 and 445 DF,  p-value: < 2.2e-16 
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_SubtypeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(subtype),data=pData(tumor6.race.combined.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score,tumor6.race.combined.eset$subtype)
#INDICES: Basal
#   var  n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 81 2.77 2.18      2    2.74 2.97  -2   6     8 0.07    -1.13 0.24

#INDICES: Her2
#   var  n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 66 0.18 2.84      0    0.19 2.97  -6   6    12  0.1    -0.96 0.35

#INDICES: LumA
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 152 -1.95 2.88     -2   -2.08 2.97  -6   6    12 0.42    -0.49 0.23

#INDICES: LumB
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 107  1.1 2.76      2    1.17 2.97  -6   6    12 -0.25    -0.23 0.27

#INDICES: Normal
#   var  n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 44 -1.23 2.8     -2   -1.39 2.97  -6   6    12 0.45    -0.51 0.42


# now do proliferation score
aov(pData(tumor6.race.combined.eset)$prolif_score~tumor6.race.combined.eset$subtype)
fit_full<- aov(prolif_score~factor(subtype),data=pData(tumor6.race.combined.eset))
summary(aov(prolif_score~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             2.7129     0.2369  11.453  < 2e-16 ***
#factor(subtype)Her2    -2.7113     0.3298  -8.220 6.85e-15 ***
#factor(subtype)LumA    -3.5777     0.2935 -12.191  < 2e-16 ***
#factor(subtype)LumB    -1.9148     0.2966  -6.456 4.52e-10 ***
#factor(subtype)Normal  -4.3133     0.3665 -11.769  < 2e-16 ***

#F-statistic: 50.93 on 4 and 290 DF,  p-value: < 2.2e-16
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SubtypeProliferationScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(prolif_score~factor(subtype),data=pData(tumor6.race.combined.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="proliferation score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$prolif_score,pData(tumor6.race.combined.eset)$subtype)
#group: Basal
#  var   n  mean   sd median trimmed  mad    min   max range  skew kurtosis   se
#1   1 46 2.71 1.46    3.2    2.82 1.34 -1.23 5.33  6.56 -0.7    -0.18 0.22
#------------------------------------------------------------ 
#group: Her2
#  var   n mean   sd median trimmed  mad    min  max range skew kurtosis   se
#1   1 49    0 1.69    0.2    0.14 1.71 -4.82 2.73  7.55 -0.84     0.38 0.24
#------------------------------------------------------------ 
#group: LumA
#  var   n   mean   sd median trimmed  mad    min  max range  skew kurtosis   se
#1   1 86 -0.86 1.71  -0.85   -0.88 1.38 -4.66 4.71  9.37 0.19     0.64 0.18
#------------------------------------------------------------ 
#group: LumB
#  var   n mean   sd median trimmed mad    min   max range skew kurtosis  se
#1   1 81  0.8 1.63   0.91    0.83 1.98 -2.5 4.21  6.71 -0.1    -0.87 0.18
#------------------------------------------------------------ 
#group: Normal
#  var   n   mean   sd median trimmed  mad    min   max range skew kurtosis   se
#1   1 33 -1.6 1.31  -1.83   -1.65 1.09 -4.12 1.74  5.87 0.39        0 0.23


######################## look at average score by grade

aov(pData(tumor6.race.combined.eset)$d_score~pData(tumor6.race.combined.eset)$grade)
fit_full<- aov(d_score~factor(grade),data=pData(tumor6.race.combined.eset))

summary(aov(d_score~factor(grade),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#does the mean score vary by subtype
#(Intercept)     -1.9775     0.3141  -6.296 7.41e-10 ***
#factor(grade)2   1.3671     0.3945   3.465 0.000582 ***
#factor(grade)3   3.3745     0.3778   8.931  < 2e-16 ***
#F-statistic:  45.2 on 2 and 439 DF,  p-value: < 2.2e-16

table(tumor6.race.combined.eset$grade)
#  1   2   3  
#  89 154 199 

describeBy(pData(tumor6.race.combined.eset)$d_score,tumor6.race.combined.eset$grade)
#INDICES: 1
#   var   n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 89 -1.98 2.8     -2   -2.14 2.97  -6   6    12 0.57    -0.28 0.3
#INDICES: 2
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 154 -0.61 3.18      0   -0.69 2.97  -6   6    12 0.25    -0.57 0.26
#-------------------------------------------------------------------------------------------- 
#INDICES: 3
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 199  1.4 2.86      2    1.53 2.97  -6   6    12 -0.44     -0.2 0.2
#####################################################
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_GradeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(grade),data=pData(tumor6.race.combined.eset), col=c("grey"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

######################## look at average score by size
t.test(d_score~size,data=pData(tumor6.race.combined.eset))
#t = -4.2726, df = 318.728, p-value = 2.554e-05
#95 percent confidence interval:
# -1.9686474 -0.7272454
#      -0.8789809       0.4689655 

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_SizeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~size,data=pData(tumor6.race.combined.eset), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score,tumor6.race.combined.eset$size)
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 157 -0.88 3.19      0   -0.96 2.97  -6   6    12 0.16    -0.68 0.25

#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 290 0.47 3.17      0    0.48 2.97  -6   6    12 -0.1    -0.81 0.19


######################## look at average score by age

aov(pData(tumor6.race.combined.eset)$d_score~pData(tumor6.race.combined.eset)$age_cat)
fit_full<- aov(d_score~factor(age_cat),data=pData(tumor6.race.combined.eset))

summary(aov(d_score~factor(age_cat),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#does the mean score vary by subtype
#                        Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.3333     1.3214   0.252    0.801
#factor(age_cat)[30,40)    0.3153     1.3739   0.230    0.819
#factor(age_cat)[40,50)   -0.3782     1.3390  -0.282    0.778
#factor(age_cat)[50,60)   -0.4039     1.3672  -0.295    0.768
#factor(age_cat)[60,70)   -0.4242     1.4907  -0.285    0.776
#factor(age_cat)[70,100)  -1.0333     1.4170  -0.729    0.466
#F-statistic: 1.001 on 5 and 444 DF,  p-value: 0.4168
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_AgeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~age_cat,data=pData(tumor6.race.combined.eset), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score,tumor6.race.combined.eset$age_cat)
#INDICES: [0,30)
#   var n mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 6 0.33 3.2      0    0.33 4.45  -4   4     8 0.02    -1.82 1.31
#---------------------------------------------------------------------------------------- 
#INDICES: [30,40)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 74 0.65 2.91      0     0.7 2.97  -6   6    12 -0.12    -0.56 0.34
#---------------------------------------------------------------------------------------- 
#INDICES: [40,50)
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 223 -0.04 3.26      0   -0.12 2.97  -6   6    12 0.09    -0.81 0.22
#---------------------------------------------------------------------------------------- 
#INDICES: [50,60)
#   var  n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 85 -0.07 3.21      0    0.03 2.97  -6   6    12 -0.22    -0.65 0.35
#---------------------------------------------------------------------------------------- 
#INDICES: [60,70)
#   var  n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 22 -0.09 3.24      0       0 4.45  -6   4    10 -0.19    -1.39 0.69
#---------------------------------------------------------------------------------------- 
#INDICES: [70,100)
#   var  n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 40 -0.7 3.72     -1   -0.81 4.45  -6   6    12 0.22     -1.2 0.59
#

######################## look at average score by node
#remove those with value = 2; not sure what that is
node.idx<-which(as.numeric(as.character(tumor6.race.combined.eset$node))>1)
t.test(d_score~node,data=pData(tumor6.race.combined.eset[,-node.idx]))
#t = -0.3218, df = 440.951, p-value = 0.7477
#95 percent confidence interval:
#  -0.6983883  0.5018417
#    -0.05240175      0.04587156 

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/NKIUNC337_NodeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~node,data=pData(tumor6.race.combined.eset[,-node.idx]), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()


describeBy(pData(tumor6.race.combined.eset)$d_score,tumor6.race.combined.eset$node)
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 229 -0.05 3.46      0    0.01 2.97  -6   6    12 -0.11    -0.98 0.23

#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 218 0.05 2.99      0   -0.02 2.97  -6   6    12 0.13    -0.62 0.2

#INDICES: 2
#   var n mean sd median trimmed mad min max range skew kurtosis se
#V1   1 1    6 NA      6       6   0   6   6     0   NA       NA NA

####################### survival by score in quantiles in ALL


all.surv.death <- coxph(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ f.dscore.c)

summary(all.surv.death)
#f.dscore.c[-2,3)     2.369     0.4222     1.215     4.617
#f.dscore.c[3,7)      5.273     0.1896     2.607    10.666

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/All_TrainingDataSurvival.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
 
dev.off()


######################## survival by score in quantiles in Luminal A/Luminal B


all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ lumAB$f.dscore.c)

summary(all.surv.death)
#lumAB$f.dscore.c[-2,3) 0.3783    1.4598   0.4063 0.931    0.352
#lumAB$f.dscore.c[3,7)  0.8088    2.2453   0.5005 1.616    0.106

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/LumAB_TrainingDataSurvival.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ lumAB$f.dscore.c),
     lty = 1:3,col=c("darkblue","black","darkred"), mark.time = TRUE)

#legend = c(paste("Bottom 25% (ref), N = ",length(which(factor(tumor6.race.combined.eset$f.dscore.c)=='[-7,-2)')),sep=""), 
#paste("Middle 50%, N = ",length(which(factor(tumor6.race.combined.eset$f.dscore.c)=='[-2,3)')),"; HR=", round(summary(all.surv.death)$coef[3],2)," P=", round(summary(all.surv.death)$coef[9],2),sep=""),
#paste("Top 25%, N = ",length(which(factor(tumor6.race.combined.eset$f.dscore.c)=='[3,7)')),"; HR=", round(summary(all.surv.death)$coef[4],2)," P=", round(summary(all.surv.death)$coef[10],2),sep="")),
# fill = c("darkblue","black","darkred"),lty = c(1,2,3), bty = "n",cex=0.7)    
dev.off()





#now stratify by race
cau<-tumor6.race.combined.eset[,which(tumor6.race.combined.eset$Race%in%c("C"))]
aa<-tumor6.race.combined.eset[,which(tumor6.race.combined.eset$Race%in%c("AA"))]

#calculate correlation of gene expression overall
library(psych)
corr.test(t(exprs((tumor6.race.combined.eset))))
#Correlation matrix 
#       ACOX2 CRYBB2  MUC1  PSPH  SQLE  TYMS
#ACOX2   1.00  -0.09  0.29 -0.06 -0.21 -0.25
#CRYBB2 -0.09   1.00 -0.19  0.27  0.14  0.12
#MUC1    0.29  -0.19  1.00 -0.06 -0.03 -0.41
#PSPH   -0.06   0.27 -0.06  1.00  0.23  0.20
#SQLE   -0.21   0.14 -0.03  0.23  1.00  0.36
#TYMS   -0.25   0.12 -0.41  0.20  0.36  1.00

#       ACOX2 CRYBB2 MUC1 PSPH SQLE TYMS
#ACOX2   0.00   0.18 0.00 0.64 0.00 0.00
#CRYBB2  0.05   0.00 0.00 0.00 0.02 0.04
#MUC1    0.00   0.00 0.00 0.64 0.64 0.00
#PSPH    0.23   0.00 0.21 0.00 0.00 0.00
#SQLE    0.00   0.00 0.59 0.00 0.00 0.00
#TYMS    0.00   0.01 0.00 0.00 0.00 0.00

#calculate correlation of gene expression in CAU
corr.test(t(exprs((cau))))
#Correlation matrix 
#       ACOX2 CRYBB2  MUC1  PSPH  SQLE  TYMS
#ACOX2   1.00  -0.05  0.31 -0.03 -0.24 -0.23
#CRYBB2 -0.05   1.00  0.01  0.21  0.10  0.06
#MUC1    0.31   0.01  1.00 -0.03 -0.05 -0.49
#PSPH   -0.03   0.21 -0.03  1.00  0.17  0.20
#SQLE   -0.24   0.10 -0.05  0.17  1.00  0.36
#TYMS   -0.23   0.06 -0.49  0.20  0.36  1.00

#Probability values (Entries above the diagonal are adjusted for multiple tests.) 
#       ACOX2 CRYBB2 MUC1 PSPH SQLE TYMS
#ACOX2   0.00   1.00 0.00    1 0.00    0
#CRYBB2  0.37   0.00 1.00    0 0.36    1
#MUC1    0.00   0.86 0.00    1 1.00    0
#PSPH    0.53   0.00 0.62    0 0.01    0
#SQLE    0.00   0.05 0.34    0 0.00    0
#TYMS    0.00   0.27 0.00    0 0.00    0



#calculate correlation of gene expression in AA
corr.test(t(exprs((aa))))
#Correlation matrix 
#       ACOX2 CRYBB2  MUC1  PSPH  SQLE  TYMS
#ACOX2   1.00  -0.04  0.17  0.02 -0.07 -0.24
#CRYBB2 -0.04   1.00 -0.12 -0.11 -0.11 -0.10
#MUC1    0.17  -0.12  1.00  0.25  0.26 -0.04
#PSPH    0.02  -0.11  0.25  1.00  0.16 -0.10
#SQLE   -0.07  -0.11  0.26  0.16  1.00  0.21
#TYMS   -0.24  -0.10 -0.04 -0.10  0.21  1.00

#Probability values (Entries above the diagonal are adjusted for multiple tests.) 
#       ACOX2 CRYBB2 MUC1 PSPH SQLE TYMS
#ACOX2   0.00   1.00 1.00 1.00 1.00    1
#CRYBB2  0.80   0.00 1.00 1.00 1.00    1
#MUC1    0.23   0.38 0.00 1.00 0.83    1
#PSPH    0.86   0.45 0.07 0.00 1.00    1
#SQLE    0.60   0.44 0.06 0.24 0.00    1
#TYMS    0.08   0.49 0.80 0.48 0.13    0


######################## survival by score in quantiles in Caucasians
cau.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ cau$f.dscore.c)

summary(cau.surv.death)
#cau$f.dscore.c[-2,3)     2.167     0.4614     1.103     4.256
#cau$f.dscore.c[3,7)      4.772     0.2096     2.307     9.869

#cau$f.dscore.c[-2,3) 0.7734    2.1671   0.3444 2.246   0.0247 *  
#cau$f.dscore.c[3,7)  1.5627    4.7716   0.3708 4.215  2.5e-05 ***

dim(cau)
#Features  Samples 
#       6      397 
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/CAU_All_TrainingDataSurvival.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ cau$f.dscore.c),
      lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
      dev.off()

######################## survival by score in quantiles in AA
aa.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ aa$f.dscore.c)
table(aa$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#      0      32      21 
summary(aa.surv.death)
#cau$f.dscore.c[-2,3)     2.167     0.4614     1.103     4.256
#cau$f.dscore.c[3,7)      4.772     0.2096     2.307     9.869


tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CAU_Survival.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ cau$f.dscore.c),
     lty = 1:3,col=c("darkblue","black","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="Caucasian")
     
     
legend(0, 0.35, 
legend = c(paste("Bottom 25% (ref), N = ",length(which(factor(cau$f.dscore.c)=='[-7,-2)')),sep=""), 
paste("Middle 50%, N = ",length(which(factor(cau$f.dscore.c)=='[-2,3)')),"; HR=", round(summary(cau.surv.death)$coef[3],2)," P=", round(summary(cau.surv.death)$coef[9],2),sep=""),
paste("Top 25%, N = ",length(which(factor(cau$f.dscore.c)=='[3,7)')),"; HR=", round(summary(cau.surv.death)$coef[4],2)," P=", round(summary(cau.surv.death)$coef[10],2),sep="")),

 fill = c("darkblue","black","darkred"),lty = c(1,2,3), bty = "n",cex=0.7)     

dev.off()



######################## do a little bit of testing.

hist(pData(tumor6.race.combined.eset[,which(tumor6.race.combined.eset$Race%in%c("C"))])$score)
mean(pData(tumor6.race.combined.eset[,which(tumor6.race.combined.eset$Race%in%c("C"))])$score)



cyrbb2.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$CRYBB2)) #+all_caldas$AGE)
                    
summary(cyrbb2.surv.death)
#                                               exp(coef) exp(-coef) lower .95 upper .95
#factor(cau$CRYBB2)1     1.366     0.7323    0.9035     2.064

########### seems to be no effect in African Americans - just overall poor survival
cyrbb2.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +factor(aa$CRYBB2)) #+all_caldas$AGE)
                    
summary(cyrbb2.surv.death)
#factor(aa$CRYBB2)1    0.8331        1.2    0.2698     2.572



###############seems like being high is protective in caucasian, but not in AA.
acox2.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$ACOX2)) #+all_caldas$AGE)
  
summary(acox2.surv.death)
#factor(cau$ACOX2)1    0.6553      1.526    0.4325     0.993

#NOW look in AA - it is the exact opposite phenomena
acox2.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +factor(aa$ACOX2)) #+all_caldas$AGE)
                    
#make sure we are doing everything correctly. (look at #ACOX2)
summary(acox2.surv.death)
#factor(aa$ACOX2)1     0.637       1.57    0.2236     1.815



tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/AA_ACOX2_Survival.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$ACOX2)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by ACOX2 in AA")

dev.off()


##############MUC1 in Caucasians/AA

muc1.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$MUC1))
                    
#make sure we are doing everything correctly. (look at #MUC1)
summary(muc1.surv.death)
#factor(cau$MUC1)1    0.6871      1.455    0.4398     1.073



muc1.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +factor(aa$MUC1)) #+all_caldas$AGE)
                    
#same result in AA
summary(muc1.surv.death)
#factor(aa$MUC1)1    0.6499      1.539    0.2109     2.002

#############SQLE in Caucasians/AA

sqle.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$SQLE)) 
                    
#make sure we are doing everything correctly. (look at SQLE)
summary(sqle.surv.death)
#factor(cau$SQLE)1     1.707     0.5859     1.082     2.691


#the effect is 2* as strong in AA
sqle.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +factor(aa$SQLE)) 
                    
#make sure we are doing everything correctly. (look at SQLE)
summary(sqle.surv.death)
#                 exp(coef) exp(-coef) lower .95 upper .95
#factor(aa$SQLE)1     3.163     0.3161     1.193     8.385



tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/AA_SQLE_Survival.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$SQLE)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by SQLE in AA")

dev.off()

#############TYMS in Caucasians/AA

tyms.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$TYMS)) #
                    
#make sure we are doing everything correctly. (look at #TYMS)
summary(tyms.surv.death)
#factor(cau$TYMS)1     1.862     0.5371     1.132     3.063

#can't do this among AA since everyone falls in one group

#there appears to be some EMM - genes act differently by race. (ACOX2 in particular), enhancement in SQLE (same direction), no effect in crybb2

# now look at by sum of genes
cau.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +factor(cau$d_score))
                    
#this was an accident, but interesting
summary(cau.surv.death)
#                      exp(coef) exp(-coef) lower .95 upper .95
#factor(cau$d_score)-3    0.8481     1.1791    0.1895     3.795
#factor(cau$d_score)-1    3.1603     0.3164    0.9181    10.878
#factor(cau$d_score)1     1.7421     0.5740    0.5112     5.937
#factor(cau$d_score)3     3.0432     0.3286    0.9268     9.993
#factor(cau$d_score)5     4.4464     0.2249    1.3141    15.045


cau.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +cau$f.dscore.c)

summary(cau.surv.death)
# exp(coef) exp(-coef) lower .95 upper .95
#cau$f.dscore.c[1,6)     1.524     0.6563    0.9487     2.447

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/Caucasian_TwoGroups_Discrete_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ factor(cau$f.dscore.c)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by discrete score, two groups, Caucasian")

dev.off()




cau.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +cau$f.score.c)

summary(cau.surv.death)
# 
#                        exp(coef) exp(-coef) lower .95 upper .95
#cau$f.score.c[-0.04,19)     1.794     0.5574     1.109     2.902

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/Caucasian_TwoGroups_Continuous_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ factor(cau$f.score.c)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by continuous score, two groups, Caucasian")

dev.off()



##################################################################################################
# now look at African Americans - generate sum score with genes that are important
#ACOX2 (up=bad), MUC1 (up = good), SQLE (up = bad)
###################################################################################################


scale_gene <- matrix(c(1,-1,1), ncol = 3)

colnames(pData(aa))
d_score<-scale_gene%*%t(pData(aa)[,c(6,8,9)])
cor(aa$ACOX2,aa$MUC1)
cor(aa$ACOX2,aa$SQLE)


quantile(d_score)
#  0%  25%  50%  75% 100% 
#   -3   -1   -1    1    1 

pData(aa)$d_score<-t(d_score)
dscore_cat_c<-cut(aa$d_score,breaks=c(-4,-1,3),right=TRUE) # corresponds to 50%
f.dscore.c<-factor(dscore_cat_c)

pData(aa)$f.dscore.c<-f.dscore.c

table(pData(aa)$f.dscore.c)
#(-4,-1]  (-1,3] 
#     37      16 

score<-scale_gene%*%exprs(aa[c(1,3,4)])
quantile(score)
#         0%         25%         50%         75%        100% 
#				-6.892176 -2.607176 -1.400176 -0.167176  2.461824 

cor(aa$score,aa$d_score)
pData(aa)$score<-t(score)

score_cat_c<-cut(aa$score,breaks=c(-7,-1.3,3),right=FALSE) # corresponds to quantiles
f.score.c<-factor(score_cat_c)
pData(aa)$f.score.c<-f.score.c
table(pData(aa)$f.score.c)

#now look at survival in aa with their important genes

aa.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +aa$f.dscore.c)

summary(aa.surv.death)
# exp(coef) exp(-coef) lower .95 upper .95
#aa$f.dscore.c(-1,3]     3.962     0.2524      1.52     10.33

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/AA_TwoGroups_Discrete_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$f.dscore.c)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by discrete score, two groups, AA")

dev.off()




aa.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +aa$f.score.c)

summary(aa.surv.death)
# 
#                        exp(coef) exp(-coef) lower .95 upper .95
#cau$f.score.c[-0.04,19)     1.794     0.5574     1.109     2.902

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/AA_TwoGroups_Continuous_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$f.score.c)),
     lty = 1:2,col=c("darkblue","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by continuous score, two groups, AA")

dev.off()




####################################################################################################
#survival analysis, stratified by quantile score, among all

####################################################################################################
#survival analysis, stratified by quantile score, among Luminal A
#tumor10.race.scale.combined.eset<-tumor.race.scale.combined.eset[,as.numeric(as.character(tumor.race.scale.combined.eset$S_Time))<=10]

aa<-tumor.race.scale.combined.eset[,which(tumor.race.scale.combined.eset$Race%in%c("AA"))]

aa.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +aa$f.score.c) #+all_caldas$AGE)
summary(aa.surv.death)

#tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/ACOX2_Survival_Race.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile (AA)")


all.surv.death <- coxph(Surv(as.numeric(as.character(tumor.race.scale.combined.eset$S_Time)),as.numeric(as.character(tumor.race.scale.combined.eset$S_Event))) ~
                    +tumor.race.scale.combined.eset$f.score.c) #+all_caldas$AGE)
summary(all.surv.death)

#tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/ACOX2_Survival_Race.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(tumor.race.scale.combined.eset$S_Time)),as.numeric(as.character(tumor.race.scale.combined.eset$S_Event))) ~ factor(tumor.race.scale.combined.eset$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile")
     
     
     
     
     
     
#legend(1, 0.4, 
#legend = c(paste("ACOX2 -/Caucasian, N = ",length(which(factor(tumor.race.scale.combined.eset$race_acox2_status)==1)),sep=""), 
#paste("ACOX2 -/AA, N = ",length(which(factor(tumor.race.scale.combined.eset$race_acox2_status)==2)),sep=""),
#paste("ACOX2 +/Caucasian, N = ",length(which(factor(tumor.race.scale.combined.eset$race_acox2_status)==3)),sep=""),
#paste("ACOX2 +/AA, N = ",length(which(factor(tumor.race.scale.combined.eset$race_acox2_status)==4)),sep="")),
# fill = c("darkblue","blue","red","darkred"),lty = c(1,2), bty = "n",cex=0.7)
       

		text(14,0.35,paste("HR: 2:1 = ", round(summary(all.surv.death)$coef[4],4)," P: ", round(summary(all.surv.death)$coef[13],4),sep=' '),cex=0.8,col="blue")
		text(14,0.31,paste("HR: 3:1 = ", round(summary(all.surv.death)$coef[5],4)," P: ", round(summary(all.surv.death)$coef[14],4),sep=' '),cex=0.8,col="red")
		text(14,0.27,paste("HR: 4:1 = ", round(summary(all.surv.death)$coef[6],4)," P: ", round(summary(all.surv.death)$coef[15],4),sep=' '),cex=0.8,col="darkred")
      dev.off()



####################################################################################################
#survival analysis, stratified by quantile score, among African Americans
#tumor10.race.scale.combined.eset<-tumor.race.scale.combined.eset[,as.numeric(as.character(tumor.race.scale.combined.eset$S_Time))<=10]

aa<-tumor.race.scale.combined.eset[,which(tumor.race.scale.combined.eset$Race%in%c("AA"))]

aa.surv.death <- coxph(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~
                    +aa$f.score.c) #+all_caldas$AGE)
summary(aa.surv.death)

#tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/ACOX2_Survival_Race.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa$S_Time)),as.numeric(as.character(aa$S_Event))) ~ factor(aa$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile (AA)")



####################################################################################################
#survival analysis, stratified by quantile score, among Caucasiana


cau<-tumor.race.scale.combined.eset[,which(tumor.race.scale.combined.eset$Race%in%c("C"))]

cau.surv.death <- coxph(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~
                    +cau$f.score.c) #+all_caldas$AGE)
summary(cau.surv.death)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/Score_Survival_Cau.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(cau$S_Time)),as.numeric(as.character(cau$S_Event))) ~ factor(cau$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile (Caucasian)")

dev.off()

















#################### create median centered dataset
#median center Luminal A separatelty
lumA<-tumor5.race.combined.eset[,which(tumor5.race.combined.eset$subtype%in%c("LumA"))]
medians.lumA <- esApply(lumA,1,median)
lumA.exprs.scale <- t(scale(t(exprs(lumA)),center=medians.lumA,scale=FALSE))

## Create a new eset
lumA.scale.eset<-lumA
exprs(lumA.scale.eset) <- lumA.exprs.scale 
dim(lumA.scale.eset)


ACOX2<-sapply(exprs(lumA.scale.eset)[1,],castDirection)
MUC1<-sapply(exprs(lumA.scale.eset)[2,],castDirection)
SQLE<-sapply(exprs(lumA.scale.eset)[3,],castDirection)
TYMS<-sapply(exprs(lumA.scale.eset)[4,],castDirection)

table(ACOX2)
table(MUC1)
table(SQLE)
table(TYMS)

lumA.scale.eset$ACOX2<-ACOX2
lumA.scale.eset$MUC1<-MUC1
lumA.scale.eset$SQLE<-SQLE
lumA.scale.eset$TYMS<-TYMS





rownames(exprs(lumA.scale.eset))
#[1] "ACOX2"  "CRYBB2" "MUC1"   "SQLE"   "TYMS" 
#create a "ACOX2" "MUC1"  "SQLE"  "TYMS"  variable based on median centered data
# ACOX2 - up is good -1* exprs
# CYRYBB2 - up is bad 1*exprs
# MUC1 - up is good-1*exprs
# SQLE - up is bad 1*exprs
# TYMS - up is bad 1*exprs

scale_gene <- matrix(c(-1,-1,1,1), ncol = 4)
colnames(pData(lumA.scale.eset))


#discrete score (should vary between -4->4)
d_score<-scale_gene%*%t(pData(lumA.scale.eset)[,6:9])
quantile(d_score)
#  0%  25%  50%  75% 100% 
#  -4   -4    0    4    4 

pData(lumA.scale.eset)$d_score<-t(d_score)

t.test(pData(lumA.scale.eset)$score~factor(pData(lumA.scale.eset)$Race))


table(pData(lumA.scale.eset)$d_score,pData(lumA.scale.eset)$Race)
#    
#     AA  C
#  -4  6 40
#  -2  9 11
#  0   3 11
#  2   0 32
#  4   0 40

#extract out those individuals with -4, 4
#extremes:
lumA.x.scale.eset<-lumA.scale.eset[,which(abs(d_score)>3)]

all.surv.death <- coxph(Surv(as.numeric(as.character(lumA.x.scale.eset$S_Time)),as.numeric(as.character(lumA.x.scale.eset$S_Event))) ~
                    +factor(lumA.x.scale.eset$d_score)) #+all_caldas$AGE)
                    
#make sure we are doing everything correctly. (look at #ACOX2)
summary(all.surv.death)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/LumA5_Survival_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(lumA.scale.eset$S_Time)),as.numeric(as.character(lumA.scale.eset$S_Event))) ~ factor(lumA.scale.eset$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile, Luminals")

dev.off()




# now look at continuous score.

score<-scale_gene%*%exprs(lumA.scale.eset)
quantile(score)
#         0%         25%         50%         75%        100% 
#-12.9423410  -7.6900910   0.0128155   0.9764888   2.9512310 
pData(lumA.scale.eset)$score<-t(score)

score_cat_c<-cut(lumA.scale.eset$score,breaks=c(-13,0.013,3),right=FALSE) # corresponds to quantiles
f.score.c<-factor(score_cat_c)

pData(lumA.scale.eset)$f.score.c<-f.score.c


dim(lumA.scale.eset)
t.test(pData(lumA.scale.eset)$score~factor(pData(lumA.scale.eset)$Race))
#t = -2.807, df = 44.562, p-value = 0.00739
#mean in group AA  mean in group C 
#        -4.672730        -2.735255 



all.surv.death <- coxph(Surv(as.numeric(as.character(lumA.scale.eset$S_Time)),as.numeric(as.character(lumA.scale.eset$S_Event))) ~
                    +lumA.scale.eset$f.score.c) #+all_caldas$AGE)
summary(all.surv.death)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/LumA5_Survival_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(lumA.scale.eset$S_Time)),as.numeric(as.character(lumA.scale.eset$S_Event))) ~ factor(lumA.scale.eset$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile, Luminals")

dev.off()


lumA.AA.scale.eset<-lumA.scale.eset[,which(lumA.scale.eset$Race%in%c("AA"))]
dim(lumA.AA.scale.eset)

aa.surv.death <- coxph(Surv(as.numeric(as.character(lumA.AA.scale.eset$S_Time)),as.numeric(as.character(lumA.AA.scale.eset$S_Event))) ~
                    +lumA.AA.scale.eset$f.score.c) #+all_caldas$AGE)
summary(aa.surv.death)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/LumA5_Survival_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(aa.surv.death$S_Time)),as.numeric(as.character(aa.surv.death$S_Event))) ~ factor(aa.surv.death$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile, Luminals, AA")

dev.off()

# Caucasian Luminals

lumA.C.scale.eset<-lumA.scale.eset[,which(lumA.scale.eset$Race%in%c("C"))]
dim(lumA.C.scale.eset)

c.surv.death <- coxph(Surv(as.numeric(as.character(lumA.C.scale.eset$S_Time)),as.numeric(as.character(lumA.C.scale.eset$S_Event))) ~
                    +lumA.C.scale.eset$f.score.c) #+all_caldas$AGE)
summary(c.surv.death)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Images/LumA5_Survival_All.tiff',antialias="none")

plot(survfit(Surv(as.numeric(as.character(c.surv.death$S_Time)),as.numeric(as.character(c.surv.death$S_Event))) ~ factor(c.surv.death$f.score.c)),
     lty = 1:4,col=c("darkblue","blue","red","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by score quantile, Luminals,Caucasian")

dev.off()





