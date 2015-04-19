# this is extension of RaceAggregateScore.R
# this is the extended analysis in the UNC337+NKI295; 
# evaluate association of points with tumor characteristics.  analog to Caldas1584Analysis_Extended.R
# analysis begins line 281
# subtype analyses begin 361
# age analysis on line 484
#size associations line 630
# grade associations line 698
# Node association line 840
# survival analyses line 963
library(Biobase)
library(gdata) #
library(hgug4112a.db)
library(gplots)
library(limma)
library(gage)
library(survival)
library(muhaz)
library(psych)
data(egSymb)



race_symbols<-c("CRYBB2","MUC1", "ACOX2", "SQLE", "TYMS", "PSPH")
race_entrez<-sym2eg(race_symbols)

#MD - update 1/2014.  calculate progression score 
proliferationGenes<-c("CCNB1","UBE2C","BIRC5","KNTC2","CDC20","PTTG1","RRM2","MKI67","TYMS","CEP55","CDCA1")
proliferation_entrez<-sym2eg(proliferationGenes)

source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/pcaFuncs.R")
source("/Users/mdarcy/Desktop/MTroester/AgePaper/AgeManuscript_062011/age_paper_scripts.R")

#  set working directory
setwd("/Users/mdarcy/Desktop/UNC337_Race/")

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
#Features  Samples (there are only 5 progression genes) 
#       10      165 
rownames(exprs(race.337.eset))<-fData(race.337.eset)$symbol

race.337.eset$race2<-race.337.eset$Race
table(race.337.eset$Race)
head(race.337.eset$GEO.array.names)
#write out the GEO codes and race
write.table(pData(race.337.eset)[,c(1,30)],file="/Users/mdarcy/Desktop/UNC337_Race/UNC337.Race.txt",sep="\t",col.names=TRUE)

colnames(pData(race.337.eset))

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

#create grade variable
table(race.337.eset$Grade)
# 1  2  3 
#14 53 80 

table(race.337.eset$PAM50.Call)
# Basal   Her2   LumA   LumB Normal 
#    35     17     66     26     11 

# we will have to look at [32] "RFS_event"                                                                                                                                                     
# [33] "RFS_months"   instead of overall survival  - this may explain the difference
#[34] "Dead_of_Disease"                                                                                                                                               
#[35] "OS_event"  
#[36] "DOD_and_OS_months" 
table(race.337.eset$Dead_of_Disease)
#  0   1 
#126  29 

table(race.337.eset$S_Event)
#  0   1 
#126  29 

table(race.337.eset$RFS_event)
# 0   1 
#109  46 

table(race.337.eset$Dead_of_Disease,race.337.eset$S_Event)
#      0   1
#  0 122   4
#  1   4  25

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

race.nki295.eset <- nki295.eset[fData(nki295.eset)$symbol %in% unique(c(race_symbols,proliferationGenes)),]
rownames(exprs(race.nki295.eset))
 #[1] "SQLE"   "MKI67"  "MUC1"   "PTTG1"  "ACOX2"  "PSPH"   "KNTC2"  "UBE2C"  "CCNB1"  "RRM2"   "CEP55"  "RRM2"  
#[13] "TYMS"   "CDCA1"  "BIRC5"  "CRYBB2" "CDC20" 

#there are duplicates
exprs(race.nki295.eset)<-collapseIDs(exprs(race.nki295.eset),fData(race.nki295.eset)$symbol,"mean")
dim(exprs(race.nki295.eset))
#[1]  16 295


# need to discuss with MT which survival to use
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
#       6      295 

table(pData(race.nki295.eset)$subtype)
# Basal   Her2   LumA   LumB Normal 
#    46     49     86     81     33

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

#[1] 16 450

tumor6.race.combined.eset <- new("ExpressionSet", exprs = data.matrix(two.data),
                   phenoData = combined.phenoData, annotation = "hgug4112a")

save(tumor6.race.combined.eset,file="tumor6_race_combined_eset.RData")

########### ########### ########### ########### ########### ########### ###########
# now create 1/-1 variables for 
########### ########### ########### ########### ########### ########### ###########
rownames(exprs(tumor6.race.combined.eset))
# [1] "ACOX2"  "BIRC5"  "CCNB1"  "CDC20"  "CDCA1"  "CEP55"  "CRYBB2" "KNTC2"  "MKI67"  "MUC1"   "PSPH"   "PTTG1" 
#[13] "RRM2"   "SQLE"   "TYMS"   "UBE2C" 
proliferation.idx<-which(rownames(exprs(tumor6.race.combined.eset))%in%proliferationGenes)
proliferation.idx
prolif_score<-apply(exprs(tumor6.race.combined.eset)[proliferation.idx,],2,sum)
pData(tumor6.race.combined.eset)$prolif_score<-prolif_score


rownames(exprs(tumor6.race.combined.eset))
#[ [1] "ACOX2"  "BIRC5"  "CCNB1"  "CDC20"  "CDCA1"  "CEP55"  "CRYBB2" "KNTC2"  "MKI67"  "MUC1"   "PSPH"   "PTTG1" 
#[13] "RRM2"   "SQLE"   "TYMS"   "UBE2C" 
#create a "ACOX2" "MUC1"  "SQLE"  "TYMS"  variable based on median centered data
# ACOX2 - up is good -1* exprs
# CYRYBB2 - up is bad 1*exprs
# MUC1 - up is good-1*exprs
# SQLE - up is bad 1*exprs
# TYMS - up is bad 1*exprs
#PSPH - up is bad 1*expr

ACOX2<-sapply(exprs(tumor6.race.combined.eset)[1,],castDirection)
CRYBB2<-sapply(exprs(tumor6.race.combined.eset)[7,],castDirection)
MUC1<-sapply(exprs(tumor6.race.combined.eset)[10,],castDirection)
PSPH<-sapply(exprs(tumor6.race.combined.eset)[11,],castDirection)
SQLE<-sapply(exprs(tumor6.race.combined.eset)[14,],castDirection)
TYMS<-sapply(exprs(tumor6.race.combined.eset)[15,],castDirection)

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


colnames(pData(tumor6.race.combined.eset))
# [1] "Data_Source"  "S_Time"       "S_Event"      "Race"         "subtype"      "grade"        "node"        
# [8] "size"         "age_cat"      "prolif_score" "ACOX2"        "CRYBB2"       "MUC1"         "SQLE"        
#[15] "PSPH"         "TYMS"     

######################################
#group 1: d_score for crybb2 + psph 
#(these are the genes differentially expressed in normal tissue too)
#######################################
scale_gene_grp1 <- matrix(c(1,1), ncol = 2)
d_score_grp1<-scale_gene_grp1%*%t(pData(tumor6.race.combined.eset)[,c(12,15)])
quantile(d_score_grp1)
#  0%  25%  50%  75% 100% 
#  -2   -2    0    2    2 
#<=-2, >=2 are good cutpoints

pData(tumor6.race.combined.eset)$d_score_grp1<-t(d_score_grp1)

######################################
#group 2: d_score for acox2 + muc1 + tyms + psph
#######################################
## [8] "size"         "age_cat"      "prolif_score" "ACOX2"        "CRYBB2"       "MUC1"         "SQLE"        
#[15] "PSPH"         "TYMS"  
scale_gene_grp2 <- matrix(c(-1,-1,1,1), ncol = 4)
d_score_grp2<-scale_gene_grp2%*%t(pData(tumor6.race.combined.eset)[,c(11,13,14,16)])
quantile(d_score_grp2)
#  0%  25%  50%  75% 100% 
#  -4   -2    0    2    4 
#<-2, >2 are good cutpoints

pData(tumor6.race.combined.eset)$d_score_grp2<-t(d_score_grp2)

######################################
#group 3: d_score for proliferation genes: tyms + psph
#######################################
scale_gene_prol <- matrix(c(1,1), ncol = 2)
d_score_prol<-scale_gene_prol%*%t(pData(tumor6.race.combined.eset)[,c(14,16)])

quantile(d_score_prol)
#  0%  25%  50%  75% 100% 
#  -2   -2    0    2    2 
#<=-2, >=2 are good cutpoints

pData(tumor6.race.combined.eset)$d_score_prol<-t(d_score_prol)



##################################################################
# what is the association between subpoints and proliferation score
# first look at correlation with proliferation score
##################################################################
# group 1
cor(pData(tumor6.race.combined.eset)$prolif_score,pData(tumor6.race.combined.eset)$d_score_grp1,use="complete.obs")
#          [,1]
#[1,][1,] 0.1718086
points_fit<-lm(pData(tumor6.race.combined.eset)$prolif_score~pData(tumor6.race.combined.eset)$d_score_grp1) 

summary(points_fit)
#pData(tumor6.race.combined.eset)$d_score_grp1  0.24020    0.08046   2.985  0.00307 **

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_ProliferationPoints_Grp1Correlation.tiff',width=350,height=350,res=100,antialias="none")
plot(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$prolif_score)
abline(points_fit)
dev.off()


aov(pData(tumor6.race.combined.eset)$prolif_score~factor(pData(tumor6.race.combined.eset)$d_score_grp1))
fit_full<- aov(prolif_score~factor(d_score_grp1),data=pData(tumor6.race.combined.eset))

summary(aov(prolif_score~factor(d_score_grp1),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            -0.2664     0.2287  -1.165  0.24494  
#factor(d_score_grp1)0   0.4708     0.2905   1.621  0.10612  
#factor(d_score_grp1)2   0.9607     0.3224   2.980  0.00312 **
#F-statistic: 4.442 on 2 and 292 DF,  p-value: 0.01258

################ group 2
cor(pData(tumor6.race.combined.eset)$prolif_score,pData(tumor6.race.combined.eset)$d_score_grp2,use="complete.obs")
#          [,1]
#[1,] 0.6125468
points_fit<-lm(pData(tumor6.race.combined.eset)$prolif_score~pData(tumor6.race.combined.eset)$d_score_grp2) 
summary(points_fit)
#pData(tumor6.race.combined.eset)$d_score_grp2  0.49391    0.03723  13.265   <2e-16 ***

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_ProliferationPoints_Grp2Correlation.tiff',width=350,height=350,res=100,antialias="none")
plot(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$prolif_score)
abline(points_fit)
dev.off()


aov(pData(tumor6.race.combined.eset)$prolif_score~factor(pData(tumor6.race.combined.eset)$d_score_grp2))
fit_full<- aov(prolif_score~factor(d_score_grp2),data=pData(tumor6.race.combined.eset))
summary(aov(prolif_score~factor(d_score_grp2),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             -1.4525     0.2484  -5.849 1.34e-08 ***
#factor(d_score_grp2)-2   0.4916     0.3178   1.547    0.123    
#factor(d_score_grp2)0    1.5097     0.3178   4.750 3.20e-06 ***
#factor(d_score_grp2)2    2.7174     0.3178   8.550 7.17e-16 ***
#factor(d_score_grp2)4    3.7552     0.3512  10.692  < 2e-16 ***
#F-statistic: 44.89 on 4 and 290 DF,  p-value: < 2.2e-16

################ group 3
cor(pData(tumor6.race.combined.eset)$prolif_score,pData(tumor6.race.combined.eset)$d_score_prol,use="complete.obs")
#          [,1]
#[1,] 0.6915943
points_fit<-lm(pData(tumor6.race.combined.eset)$prolif_score~pData(tumor6.race.combined.eset)$d_score_prol) 
summary(points_fit)
#pData(tumor6.race.combined.eset)$d_score_prol  0.87949    0.05366  16.390   <2e-16 ***

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_ProliferationPoints_ProlifGenesCorrelation.tiff',width=350,height=350,res=100,antialias="none")
plot(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$prolif_score)
abline(points_fit)
dev.off()

aov(pData(tumor6.race.combined.eset)$prolif_score~factor(pData(tumor6.race.combined.eset)$d_score_prol))
fit_full<- aov(prolif_score~factor(d_score_prol),data=pData(tumor6.race.combined.eset))
summary(aov(prolif_score~factor(d_score_prol),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            -1.7139     0.1507  -11.37   <2e-16 ***
#factor(d_score_prol)0   2.2398     0.2131   10.51   <2e-16 *** 
#factor(d_score_prol)2   3.5196     0.2126   16.56   <2e-16 ***

#F-statistic: 140.4 on 2 and 292 DF,  p-value: < 2.2e-16

##################################################################
#does the mean score vary by subtype for each of the 3 groups
##################################################################

# 1) d_score_grp1
aov(pData(tumor6.race.combined.eset)$d_score_grp1~pData(tumor6.race.combined.eset)$subtype)
fit_full<- aov(d_score_grp1~factor(subtype),data=pData(tumor6.race.combined.eset))


summary(aov(d_score_grp1~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)            0.27160    0.16552   1.641  0.10152   
#factor(subtype)Her2   -0.05948    0.24702  -0.241  0.80982  
#factor(subtype)LumA   -0.56108    0.20493  -2.738  0.00643 **
#factor(subtype)LumB   -0.06600    0.21940  -0.301  0.76370   
#factor(subtype)Normal -0.49888    0.27898  -1.788  0.07443 . 

#F-statistic: 3.204 on 4 and 445 DF,  p-value: 0.01304
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SubtypeScore_Grp1.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp1~factor(subtype),data=pData(tumor6.race.combined.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$subtype)
#INDICES: Basal
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 81 0.27 1.47      0    0.34 2.97  -2   2     4 -0.21    -1.17 0.16

#INDICES: Her2
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 66 0.21 1.57      0    0.26 2.97  -2   2     4 -0.18    -1.39 0.19

#INDICES: LumA
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 152 -0.29 1.48      0   -0.36 2.97  -2   2     4 0.23    -1.17 0.12

#INDICES: LumB
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 107 0.21 1.43      0    0.25 2.97  -2   2     4 -0.15    -1.05 0.14

#INDICES: Normal
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 44 -0.23 1.57      0   -0.28 2.97  -2   2     4 0.19    -1.39 0.24

table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$subtype)
#     Basal Her2 LumA LumB Normal
#  -2    17   17   54   22     16
#  0     36   25   66   52     17
#  2     28   24   32   33     11

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$subtype))
#X-squared = 15.0452, df = 8, p-value = 0.05827


#d_score_grp2
aov(pData(tumor6.race.combined.eset)$d_score_grp2~pData(tumor6.race.combined.eset)$subtype)
fit_full<- aov(d_score_grp2~factor(subtype),data=pData(tumor6.race.combined.eset))
summary(aov(d_score_grp2~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)

#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             2.4938     0.2334  10.683  < 2e-16 ***
#ffactor(subtype)Her2    -2.5241     0.3484  -7.245 1.92e-12 ***
#factor(subtype)LumA    -4.1517     0.2890 -14.365  < 2e-16 ***
#factor(subtype)LumB    -1.5966     0.3094  -5.160 3.73e-07 ***
#factor(subtype)Normal  -3.4938     0.3935  -8.880  < 2e-16 ***

#F-statistic: 59.57 on 4 and 445 DF,  p-value: < 2.2e-16
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SubtypeScore_Grp2.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp2~factor(subtype),data=pData(tumor6.race.combined.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$subtype)
#INDICES: Basal
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 81 2.49 1.53      2    2.68 2.97  -2   4     6 -0.77     0.08 0.17

#INDICES: Her2
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 66 -0.03 2.2      0       0 2.97  -4   4     8 -0.11    -0.72 0.27

#INDICES: LumA
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 152 -1.66 2.19     -2   -1.87 2.97  -4   4     8 0.65    -0.57 0.18

#INDICES: LumB
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 107  0.9 2.24      2    0.99 2.97  -4   4     8 -0.29    -0.72 0.22

#INDICES: Normal
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 44   -1 2.18     -1   -1.11 1.48  -4   4     8 0.32    -0.56 0.33

table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$subtype)
#     Basal Her2 LumA LumB Normal
#  -4     0    7   50    5      9
#  -2     2   14   53   17     13
#  0     10   23   25   31     15
#  2     35   17   21   33      5
#  4     34    5    3   21      2

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$subtype))
#X-squared = 183.7477, df = 16, p-value < 2.2e-16



########d_score_prol
aov(pData(tumor6.race.combined.eset)$d_score_prol~pData(tumor6.race.combined.eset)$subtype)
fit_full<- aov(d_score_prol~factor(subtype),data=pData(tumor6.race.combined.eset))
summary(aov(d_score_prol~factor(subtype),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)

#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             1.0123     0.1541   6.567 1.43e-10 ***
#factor(subtype)Her2    -0.9517     0.2301  -4.137 4.21e-05 ***
#factor(subtype)LumA    -1.8545     0.1909  -9.717  < 2e-16 ***
#factor(subtype)LumB    -0.1899     0.2043  -0.929    0.353   
#factor(subtype)Normal  -1.9669     0.2598  -7.570 2.17e-13 ***

#F-statistic: 39.42 on 4 and 445 DF,  p-value: < 2.2e-16
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SubtypeScore_Prol.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_prol~factor(subtype),data=pData(tumor6.race.combined.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$subtype)
#INDICES: Basal
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 81 1.01 1.05      2    1.05   0  -2   2     4 -0.28    -1.38 0.12

#INDICES: Her2
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 66 0.06 1.57      0    0.07 2.97  -2   2     4 -0.05     -1.4 0.19

#INDICES: LumA
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 152 -0.84 1.47     -2   -1.05   0  -2   2     4 0.83     -0.7 0.12

#INDICES: LumB
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 107 0.82 1.4      2    1.01   0  -2   2     4 -0.75    -0.68 0.14

#INDICES: Normal
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 44 -0.95 1.33     -2   -1.17   0  -2   2     4 0.85    -0.47 0.2


table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$subtype)
#     Basal Her2 LumA LumB Normal
#  -2     1   19   86   13     25
#  0     38   26   44   37     15
#  2     42   21   22   57      4

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$subtype))
#X-squared = 125.8889, df = 8, p-value < 2.2e-16


#########################################################
# does the mean score vary by age
#########################################################
#d_score_grp1 (CRYBB2, PSPH)
aov(pData(tumor6.race.combined.eset)$d_score_grp1~pData(tumor6.race.combined.eset)$age_cat)
fit_full<- aov(d_score_grp1~age_cat,data=pData(tumor6.race.combined.eset))
summary(aov(d_score_grp1~age_cat,data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      -0.6667     0.6149  -1.084    0.279
#age_cat[30,40)    0.8829     0.6394   1.381    0.168 
#age_cat[40,50)    0.6667     0.6231   1.070    0.285 
#age_cat[50,60)    0.6196     0.6363   0.974    0.331
#age_cat[60,70)    0.8485     0.6937   1.223    0.222
#age_cat[70,100)   0.4667     0.6594   0.708    0.480

#F-statistic: 0.7586 on 5 and 444 DF,  p-value: 0.5801
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_AgeScoreGrp1.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp1~age_cat,data=pData(tumor6.race.combined.eset),names=c("<30","30-39","40-49","50-59","60-69","70+"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$age_cat) 
#INDICES: [0,30)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 6 -0.67 1.63     -1   -0.67 1.48  -2   2     4 0.48    -1.58 0.67
#-------------------------------------------------------------------------------------------- 
#INDICES: [30,40)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 74 0.22 1.57      0    0.27 2.97  -2   2     4 -0.19    -1.38 0.18
#-------------------------------------------------------------------------------------------- 
#INDICES: [40,50)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 223    0 1.46      0       0 2.97  -2   2     4    0    -1.13 0.1
#-------------------------------------------------------------------------------------------- 
#INDICES: [50,60)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 85 -0.05 1.54      0   -0.06 2.97  -2   2     4 0.04    -1.34 0.17
#-------------------------------------------------------------------------------------------- 
#INDICES: [60,70)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 22 0.18 1.37      0    0.22   0  -2   2     4 -0.1    -0.97 0.29
#-------------------------------------------------------------------------------------------- 
#INDICES: [70,100)
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 40 -0.2 1.62      0   -0.25 2.97  -2   2     4 0.17    -1.49 0.26


table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$age_cat)
#     [0,30) [30,40) [40,50) [50,60) [60,70) [70,100)
#  -2      3      19      59      26       4       15
#  0       2      28     105      35      12       14
#  2       1      27      59      24       6       11

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$age_cat))
#X-squared = 8.6586, df = 10, p-value = 0.5648



#d_score_grp2 (ACOX2,MUC1,TYMS,SQLE)
aov(pData(tumor6.race.combined.eset)$d_score_grp2~pData(tumor6.race.combined.eset)$age_cat)
fit_full<- aov(d_score_grp2~age_cat,data=pData(tumor6.race.combined.eset))
summary(aov(d_score_grp2~age_cat,data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       1.0000     1.0584   0.945    0.345
#age_cat[30,40)   -0.5676     1.1004  -0.516    0.606
#age_cat[40,50)   -1.0448     1.0725  -0.974    0.330
#age_cat[50,60)   -1.0235     1.0951  -0.935    0.350 
#age_cat[60,70)   -1.2727     1.1940  -1.066    0.287
#age_cat[70,100)  -1.5000     1.1350  -1.322    0.187 

#F-statistic: 0.9514 on 5 and 444 DF,  p-value: 0.4475
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_AgeScoreGrp2.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp2~age_cat,data=pData(tumor6.race.combined.eset),names=c("<30","30-39","40-49","50-59","60-69","70+"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$age_cat) 
#INDICES: [0,30)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 6    1 2.76      1       1 4.45  -2   4     6    0    -2.06 1.13
#-------------------------------------------------------------------------------------------- 
#INDICES: [30,40)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 74 0.43 2.48      0    0.53 2.97  -4   4     8 -0.2    -0.95 0.29
#-------------------------------------------------------------------------------------------- 
#INDICES: [40,50)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 223 -0.04 2.56      0   -0.06 2.97  -4   4     8 -0.02    -1.14 0.17
#-------------------------------------------------------------------------------------------- 
#INDICES: [50,60)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 85 -0.02 2.63      0   -0.03 2.97  -4   4     8 -0.07    -1.15 0.29
#-------------------------------------------------------------------------------------------- 
#INDICES: [60,70)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 22 -0.27 2.78      0   -0.33 2.97  -4   4     8 0.03    -1.38 0.59
#-------------------------------------------------------------------------------------------- 
#INDICES: [70,100)
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 40 -0.5 2.75      0   -0.62 2.97  -4   4     8 0.27     -1.2 0.43
#-------------------------------------------------------------------------------------------- 

table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$age_cat)
#     [0,30) [30,40) [40,50) [50,60) [60,70) [70,100)
#  -4      0       8      34      15       5        9
#  -2      2      13      54      16       4       10
#  0       1      21      47      21       5        9
#  2       1      19      59      21       5        6
#  4       2      13      29      12       3        6

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$age_cat))
#X-squared = 12.2075, df = 20, p-value = 0.9087



#########d_score_prol (TYMS,SQLE)
aov(pData(tumor6.race.combined.eset)$d_score_prol~pData(tumor6.race.combined.eset)$age_cat)
fit_full<- aov(d_score_prol~age_cat,data=pData(tumor6.race.combined.eset))
summary(aov(d_score_prol~age_cat,data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      0.33333    0.65533   0.509    0.611
#age_cat[30,40)  -0.03604    0.68138  -0.053    0.958
#age_cat[40,50)  -0.39611    0.66409  -0.596    0.551
#age_cat[50,60)  -0.19216    0.67806  -0.283    0.777
#age_cat[60,70)  -0.78788    0.73931  -1.066    0.287
#age_cat[70,100) -0.53333    0.70276  -0.759    0.448

#F-statistic: 1.233 on 5 and 444 DF,  p-value: 0.2924
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_AgeScoreProl.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_prol~age_cat,data=pData(tumor6.race.combined.eset),names=c("<30","30-39","40-49","50-59","60-69","70+"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$age_cat) 
#INDICES: [0,30)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 6 0.33 1.97      1    0.33 1.48  -2   2     4 -0.25    -2.08 0.8
#-------------------------------------------------------------------------------------------- 
#INDICES: [30,40)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 74  0.3 1.58      0    0.37 2.97  -2   2     4 -0.26    -1.37 0.18
#-------------------------------------------------------------------------------------------- 
#INDICES: [40,50)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 223 -0.06 1.63      0   -0.08 2.97  -2   2     4 0.06    -1.49 0.11
#-------------------------------------------------------------------------------------------- 
#INDICES: [50,60)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 85 0.14 1.6      0    0.17 2.97  -2   2     4 -0.12    -1.44 0.17
#-------------------------------------------------------------------------------------------- 
#INDICES: [60,70)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 22 -0.45 1.74     -1   -0.56 1.48  -2   2     4 0.42     -1.6 0.37
#-------------------------------------------------------------------------------------------- 
#INDICES: [70,100)
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 40 -0.2 1.42      0   -0.25 1.48  -2   2     4 0.13    -1.06 0.22
#-------------------------------------------------------------------------------------------- 

table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$age_cat)
#     [0,30) [30,40) [40,50) [50,60) [60,70) [70,100)
#  -2      2      18      77      24      11       12
#  0       1      27      76      31       5       20
#  2       3      29      70      30       6        8


chisq.test(table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$age_cat))
#X-squared = 12.5949, df = 10, p-value = 0.2472


################################################
########does the mean score vary by size_group
################################################
#d_score_grp1
t.test(pData(tumor6.race.combined.eset)$d_score_grp1~pData(tumor6.race.combined.eset)$size)
#t = -2.6521, df = 314.94, p-value = 0.008403
# -0.6858794 -0.1016453
#mean in group 0 mean in group 1 
#     -0.2420382       0.1517241 

table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$size)
#       0   1
# -2  55  69
#  0  36 68
# 2   66 130

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$size))
#X-squared = 7.378, df = 2, p-value = 0.025


describeBy(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$size)
#INDICES: 0
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 157 -0.24 1.51      0    -0.3 2.97  -2   2     4  0.2    -1.23 0.12
#------------------------------------------------------------------ 
#INDICES: 1
#   var    n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 290 0.15 1.48      0    0.19 2.97  -2   2     4 -0.12    -1.18 0.09

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SizeScoreGrp1.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp1~size,data=pData(tumor6.race.combined.eset),names=c("< 2cm",">= 2cm"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

#d_score_grp2
t.test(pData(tumor6.race.combined.eset)$d_score_grp2~pData(tumor6.race.combined.eset)$size)
#t = -3.7797, df = 321.68, p-value = 0.0001871
# -1.4508397 -0.4575284
#mean in group 0 mean in group 1 
#     -0.6369427       0.3172414 

table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$size)
#       0   1
#  -4 34 37
#  -2 42 57
#  	0  36 68
# 	2  30 79
#  	4  15 49

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$size))
#X-squared = 17.6111, df = 2, p-value = 0.0001499

describeBy(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$size)
#INDICES: 0
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 157 -0.64 2.54      0   -0.79 2.97  -4   4     8 0.25    -1.04 0.2
#------------------------------------------------------------------ 
#INDICES: 1
#   var    n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 290 0.32 2.56      0     0.4 2.97  -4   4     8 -0.18    -1.05 0.15

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SizeScoreGrp2.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp2~size,data=pData(tumor6.race.combined.eset),names=c("< 2cm",">= 2cm"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

#d_score_prol
t.test(pData(tumor6.race.combined.eset)$d_score_prol~pData(tumor6.race.combined.eset)$size)
#t = -3.9743, df = 311.429, p-value = 8.778e-05
# -0.9394117 -0.3172542
#mean in group 0 mean in group 1 
#     -0.4076433       0.2206897 

table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$size)
#       0   1
#  -2  70  74
#  0   49 110
#  2   38 106

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$size))
#X-squared = 17.6111, df = 2, p-value = 0.0001499

describeBy(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$size)
#INDICES: 0
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 157 -0.41 1.61      0    -0.5 2.97  -2   2     4 0.38    -1.37 0.13
#------------------------------------------------------------------ 
#INDICES: 1
#   var    n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 290 0.22 1.56      0    0.28 2.97  -2   2     4 -0.19    -1.35 0.09

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_SizeScoreProl.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_prol~size,data=pData(tumor6.race.combined.eset),names=c("< 2cm",">= 2cm"), col=c("grey"),ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


#####################################################
#does the mean score vary by grade
#####################################################
# grp1 (CRYBB2, PSPH)
aov(pData(tumor6.race.combined.eset)$d_score_grp1~factor(pData(tumor6.race.combined.eset)$grade))
fit_full<- aov(d_score_grp1~factor(grade),data=pData(tumor6.race.combined.eset))
summary(aov(d_score_grp1~factor(grade),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -0.3146     0.1579  -1.992  0.04698 * 
#factor(grade)2   0.1977     0.1984   0.997  0.31946  
#factor(grade)3   0.5860     0.1900   3.084  0.00217 **
#F-statistic: 5.726 on 2 and 439 DF,  p-value: 0.003507

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp1,factor(pData(tumor6.race.combined.eset)$grade)))
#X-squared = 11.2419, df = 4, p-value = 0.02398
table(pData(tumor6.race.combined.eset)$d_score_grp1,factor(pData(tumor6.race.combined.eset)$grade))
#      1  2  3
#  -2 32 48 43
#  0  39 67 86
#  2  18 39 70
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_GradeScoreGrp1.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp1~factor(grade),data=pData(tumor6.race.combined.eset), col=c("grey"),names=c("1","2","3"), ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp1,pData(tumor6.race.combined.eset)$grade) 
#INDICES: 1
#   var   n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 89 -0.31 1.47      0   -0.38 2.97  -2   2     4 0.25    -1.15 0.16
#INDICES: 2
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 154 -0.12 1.5      0   -0.15 2.97  -2   2     4 0.09    -1.24 0.12
#-------------------------------------------------------------------------------------------- 
#INDICES: 3
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 199 0.27 1.49      0    0.34 2.97  -2   2     4 -0.22    -1.18 0.11


# grp2 (ACOX2, MUC1, TYMS, SQLE)
aov(pData(tumor6.race.combined.eset)$d_score_grp2~factor(pData(tumor6.race.combined.eset)$grade))
fit_full<- aov(d_score_grp2~factor(grade),data=pData(tumor6.race.combined.eset))
summary(aov(d_score_grp2~factor(grade),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -1.6629     0.2503  -6.643 9.11e-11 ***
#factor(grade)2   1.1694     0.3144   3.719 0.000226 ***
#factor(grade)3   2.7885     0.3011   9.260  < 2e-16 ***
#F-statistic: 48.03 on 2 and 439 DF,  p-value: < 2.2e-16

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_grp2,factor(pData(tumor6.race.combined.eset)$grade)))
#X-squared = 89.4804, df = 8, p-value = 5.928e-16
table(pData(tumor6.race.combined.eset)$d_score_grp2,factor(pData(tumor6.race.combined.eset)$grade))
#     1  2  3
#  -4 31 29 10
#  -2 30 44 24
#  0  14 36 51
#  2  10 26 72
#  4   4 19 42
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_GradeScoreGrp2.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp2~factor(grade),data=pData(tumor6.race.combined.eset), col=c("grey"),names=c("1","2","3"), ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_grp2,pData(tumor6.race.combined.eset)$grade) 
#INDICES: 1
#   var   n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 89 -1.66 2.32     -2   -1.92 2.97  -4   4     8  0.8    -0.31 0.25
#INDICES: 2
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 154 -0.49 2.57      0   -0.61 2.97  -4   4     8 0.28    -1.01 0.21
#-------------------------------------------------------------------------------------------- 
#INDICES: 3
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 199 1.13 2.21      2    1.28 2.97  -4   4     8 -0.52    -0.42 0.16


## proliferation group ( TYMS, SQLE)
aov(pData(tumor6.race.combined.eset)$d_score_prol~factor(pData(tumor6.race.combined.eset)$grade))
fit_full<- aov(d_score_prol~factor(grade),data=pData(tumor6.race.combined.eset))
summary(aov(d_score_prol~factor(grade),data=pData(tumor6.race.combined.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -0.9888     0.1573  -6.286 7.85e-10 ***
#factor(grade)2   0.7810     0.1976   3.953 9.00e-05 ***
#factor(grade)3   1.6219     0.1892   8.572  < 2e-16 ***
#F-statistic: 39.37 on 2 and 439 DF,  p-value: < 2.2e-16

chisq.test(table(pData(tumor6.race.combined.eset)$d_score_prol,factor(pData(tumor6.race.combined.eset)$grade)))
#X-squared = 75.8612, df = 4, p-value = 1.31e-15
table(pData(tumor6.race.combined.eset)$d_score_prol,factor(pData(tumor6.race.combined.eset)$grade))

#      1  2  3
#  -2 56 58 27
#  0  21 54 82
#  2  12 42 90
  tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_GradeScoreProl.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_prol~factor(grade),data=pData(tumor6.race.combined.eset), col=c("grey"),names=c("1","2","3"), ylab="disparity score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(tumor6.race.combined.eset)$d_score_prol,pData(tumor6.race.combined.eset)$grade) 
#INDICES: 1
#   var   n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 89 -0.99 1.45     -2   -1.21   0  -2   2     4 1.04    -0.37 0.15
#INDICES: 2
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 154 -0.21 1.6      0   -0.26 2.97  -2   2     4 0.19    -1.43 0.13
#-------------------------------------------------------------------------------------------- 
#INDICES: 3
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 199 0.63 1.4      0    0.78 2.97  -2   2     4 -0.52    -0.88 0.1

#################################################
# Now look at node status
#################################################
node.idx<-which(as.numeric(as.character(tumor6.race.combined.eset$node))>1)
#group1
t.test(d_score_grp1~node,data=pData(tumor6.race.combined.eset[,-node.idx]))
#t = -1.204, df = 444.393, p-value = 0.2292
#95 percent confidence interval:
# -0.4495727  0.1079999
#     -0.0698690       0.1009174 


fisher.test(table(pData(tumor6.race.combined.eset)$d_score_grp1,factor(pData(tumor6.race.combined.eset)$node)))
#p-value = 0.4063
table(pData(tumor6.race.combined.eset)$d_score_grp1,factor(pData(tumor6.race.combined.eset)$node))

#      0  1  2
#  -2 69 55  0
#  0  99 97  0
#  2  61 66  1

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_NodeScoreGrp1.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp1~node,data=pData(tumor6.race.combined.eset), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()


describeBy(pData(tumor6.race.combined.eset)$d_score_grp1,tumor6.race.combined.eset$node)
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 229 -0.07 1.51      0   -0.09 2.97  -2   2     4 0.06    -1.25 0.1

#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 218  0.1 1.49      0    0.12 2.97  -2   2     4 -0.08     -1.2 0.1

#INDICES: 2
#   var n mean sd median trimmed mad min max range skew kurtosis se
#V1   1 1    2 NA      2       2   0   2   2     0   NA       NA NA



# group2
t.test(d_score_grp2~node,data=pData(tumor6.race.combined.eset[,-node.idx]))
#t = 0.2962, df = 444.329, p-value = 0.7672
#95 percent confidence interval:
# -0.4085977  0.5536239
#     0.01746725     -0.05504587 


fisher.test(table(pData(tumor6.race.combined.eset)$d_score_grp2,factor(pData(tumor6.race.combined.eset)$node)))
#pp-value = 0.1999
table(pData(tumor6.race.combined.eset)$d_score_grp2,factor(pData(tumor6.race.combined.eset)$node))

#      0  1  2
#  -4 40 31  0
#  -2 50 49  0
#  0  44 60  0
#  2  58 51  0
#  4  37 27  1

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_NodeScoreGrp2.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_grp2~node,data=pData(tumor6.race.combined.eset), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()


describeBy(pData(tumor6.race.combined.eset)$d_score_grp2,tumor6.race.combined.eset$node)
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 229 0.02 2.7      0    0.02 2.97  -4   4     8 -0.05    -1.24 0.18

#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 218 -0.06 2.47      0   -0.07 2.97  -4   4     8 -0.01    -0.99 0.17

#INDICES: 2
#   var n mean sd median trimmed mad min max range skew kurtosis se
#V1   1 1    4 NA      4       4   0   4   4     0   NA       NA NA



#group3 - proliferation
t.test(d_score_prol~node,data=pData(tumor6.race.combined.eset[,-node.idx]))
#t = -1.3539, df = 442.823, p-value = 0.1765
#95 percent confidence interval:
# -0.50543692  0.09311369
#    -0.09606987      0.11009174 


fisher.test(table(pData(tumor6.race.combined.eset)$d_score_prol,factor(pData(tumor6.race.combined.eset)$node)))
#p-value = 0.3266
table(pData(tumor6.race.combined.eset)$d_score_prol,factor(pData(tumor6.race.combined.eset)$node))

#      0  1  2
#  -2 78 66  0
#  0  84 74  0
#  2  67 78  1

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_NodeScoreProl.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score_prol~node,data=pData(tumor6.race.combined.eset), col=c("grey"))
abline(h = 0, col = "black", lty="dotted")
dev.off()


describeBy(pData(tumor6.race.combined.eset)$d_score_prol,tumor6.race.combined.eset$node)
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 229 -0.07 1.51      0   -0.09 2.97  -2   2     4 0.06    -1.25 0.1

#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis  se
#V1   1 218  0.1 1.49      0    0.12 2.97  -2   2     4 -0.08     -1.2 0.1

#INDICES: 2
#   var n mean sd median trimmed mad min max range skew kurtosis se
#V1   1 1    2 NA      2       2   0   2   2     0   NA       NA NA



####################### survival by score in quantiles in ALL

table(pData(tumor6.race.combined.eset)$d_score_prol)
# -2   0   2 
#144 160 146  
table(pData(tumor6.race.combined.eset)$d_score_grp1)
# -2   0   2 
#126 196 128 

table(pData(tumor6.race.combined.eset)$d_score_grp2)
# -4  -2   0   2   4 
# 71  99 104 111  65 


dscore_cat_c<-cut(tumor6.race.combined.eset$d_score_grp2,breaks=c(-5,-1,2,5),right=FALSE) # corresponds to quantiles
f.dscore.c<-factor(dscore_cat_c)
table(f.dscore.c)
#[-5,-1)  [-1,2)   [2,5) 
#    170     104     176 

pData(tumor6.race.combined.eset)$f.dscore.c<-f.dscore.c

lumA<-tumor6.race.combined.eset[,which(tumor6.race.combined.eset$subtype%in%c("LumA"))]

lumAB<-tumor6.race.combined.eset[,which(tumor6.race.combined.eset$subtype%in%c("LumA","LumB"))]

table(pData(lumA)$d_score_grp1)
# -2   0   2 
#54 66 32 
table(pData(lumAB)$d_score_grp1)
# -2   0   2 
# 76 118  65 

table(pData(lumA)$d_score_grp2)
# -4  -2   0   2   4 
#	50 53 25 21  3 

table(pData(lumAB)$d_score_grp2)
# -4  -2   0   2   4 
#55 70 56 54 24 

table(pData(lumA)$f.dscore.c)
#[-5,-1)  [-1,2)   [2,5) 
#    103      25      24 

table(pData(lumAB)$f.dscore.c)
#[-5,-1)  [-1,2)   [2,5) 
#    125      56      78  

####################################################
# do survival analyses in all dataset in 3 groups
####################################################

####################################################
# first do over all tumors 
####################################################

# group1 (CRYBB2, PSPH)
all.surv.death <- coxph(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ factor(tumor6.race.combined.eset$d_score_grp1))

summary(all.surv.death)
#factor(tumor6.race.combined.eset$d_score_grp1)0 0.4664    1.5942   0.2626 1.776  0.07575 . 
#factor(tumor6.race.combined.eset$d_score_grp1)2 0.7818    2.1853   0.2699 2.897  0.00377 **
#                                                exp(coef) exp(-coef) lower .95 upper .95
#factor(tumor6.race.combined.eset$d_score_grp1)0     1.594     0.6273    0.9528     2.667
#factor(tumor6.race.combined.eset$d_score_grp1)2     2.185     0.4576    1.2877     3.709
#Score (logrank) test = 8.74  on 2 df,   p=0.01264
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_FullGrp1NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ factor(pData(tumor6.race.combined.eset)$d_score_grp1)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()

# group2 (ACOX2,MUC1,TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ factor(tumor6.race.combined.eset$f.dscore.c))

summary(all.surv.death)
#factor(tumor6.race.combined.eset$f.dscore.c)[-1,2) 0.7427    2.1017   0.2899 2.562   0.0104 *  
#factor(tumor6.race.combined.eset$f.dscore.c)[2,5)  1.1223    3.0720   0.2494 4.500  6.8e-06 ***
#                                                   exp(coef) exp(-coef) lower .95 upper .95
#factor(tumor6.race.combined.eset$f.dscore.c)[-1,2)     2.102     0.4758     1.191     3.710
#factor(tumor6.race.combined.eset$f.dscore.c)[2,5)      3.072     0.3255     1.884     5.009
#Score (logrank) test = 22.21  on 2 df,   p=1.507e-05
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_FullGrp2NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ pData(tumor6.race.combined.eset)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()


# group3 - proliferation (TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ factor(tumor6.race.combined.eset$d_score_prol))

summary(all.surv.death)
#factor(tumor6.race.combined.eset$d_score_prol)0 1.0610    2.8892   0.3010 3.525 0.000424 ***
#factor(tumor6.race.combined.eset$d_score_prol)2 1.3494    3.8551   0.2939 4.591  4.4e-06 ***
#                                                exp(coef) exp(-coef) lower .95 upper .95
#factor(tumor6.race.combined.eset$d_score_prol)0     2.889     0.3461     1.602     5.212
#factor(tumor6.race.combined.eset$d_score_prol)2     3.855     0.2594     2.167     6.858
#Score (logrank) test = 23.81  on 2 df,   p=6.765e-06
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_FullGrpProlNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(tumor6.race.combined.eset$S_Time)),as.numeric(as.character(tumor6.race.combined.eset$S_Event))) ~ factor(pData(tumor6.race.combined.eset)$d_score_prol)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()



####################################################
# now do in Luminal AB tumors
####################################################

# group1 (CRYBB2, PSPH)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ factor(lumAB$d_score_grp1))

summary(all.surv.death)
#factor(lumAB$d_score_grp1)0 0.2670    1.3061   0.3877 0.689    0.491
#factor(lumAB$d_score_grp1)2 0.2505    1.2847   0.4372 0.573    0.567
#                                                exp(coef) exp(-coef) lower .95 upper .95
#factor(lumAB$d_score_grp1)0     1.306     0.7657    0.6109     2.792
#factor(lumAB$d_score_grp1)2     1.285     0.7784    0.5453     3.027
#Score (logrank) test = 0.52  on 2 df,   p=0.7711
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_LumABGrp1NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ factor(pData(lumAB)$d_score_grp1)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()

# group2 (ACOX2,MUC1,TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ factor(lumAB$f.dscore.c))

summary(all.surv.death)
#factor(lumAB$f.dscore.c)[-1,2) 0.5708    1.7698   0.4416 1.293  0.19611
#factor(lumAB$f.dscore.c)[2,5)  0.9665    2.6287   0.3654 2.645  0.00818 **
#                                  exp(coef) exp(-coef) lower .95 upper .95
#factor(lumAB$f.dscore.c)[-1,2)     1.770     0.5650    0.7448     4.205
#factor(lumAB$f.dscore.c)[2,5)      2.629     0.3804    1.2843     5.380
#Score (logrank) test = 7.47  on 2 df,   p=0.02387
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_LumABGrp2NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ pData(lumAB)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()


# group3 - proliferation (TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ factor(lumAB$d_score_prol))

summary(all.surv.death)
#factor(lumAB$d_score_prol)0 1.0589    2.8832   0.4881 2.170  0.03004 * 
#factor(lumAB$d_score_prol)2 1.5107    4.5300   0.4631 3.262  0.00111 **
#                      		  exp(coef) exp(-coef) lower .95 upper .95
#factor(lumAB$d_score_prol)0     2.883     0.3468     1.108     7.504
#factor(lumAB$d_score_prol)2     4.530     0.2208     1.828    11.228
#Score (logrank) test = 12.44  on 2 df,   p=0.001987
#Test set - N=450
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_LumABGrpProlNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$S_Time)),as.numeric(as.character(lumAB$S_Event))) ~ factor(pData(lumAB)$d_score_prol)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()

####################################################
# now do only Luminal A tumors; N =152
####################################################

# group1 (CRYBB2, PSPH)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ factor(lumA$d_score_grp1))

summary(all.surv.death)
#factor(lumA$d_score_grp1)0 0.09282   1.09726  0.67412 0.138    0.890
#factor(lumA$d_score_grp1)2 0.39680   1.48706  0.76520 0.519    0.604
#                       exp(coef) exp(-coef) lower .95 upper .95
#factor(lumA$d_score_grp1)0     1.097     0.9114    0.2927     4.113
#factor(lumA$d_score_grp1)2     1.487     0.6725    0.3319     6.663
#Score (logrank) test = 0.29  on 2 df,   p=0.8659
#Test set -   n= 152, number of events= 12 
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_LumAGrp1NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ factor(pData(lumA)$d_score_grp1)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()

# group2 (ACOX2,MUC1,TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ factor(lumA$f.dscore.c))

summary(all.surv.death)
#factor(lumA$f.dscore.c)[-1,2) 0.8015    2.2288   0.7122 1.125    0.260
#ffactor(lumA$f.dscore.c)[2,5)  0.7158    2.0458   0.7082 1.011    0.312
#                               exp(coef) exp(-coef) lower .95 upper .95
#factor(lumA$f.dscore.c)[-1,2)     2.229     0.4487    0.5519     9.001
#factor(lumA$f.dscore.c)[2,5)      2.046     0.4888    0.5106     8.197
#Score (logrank) test = 1.81  on 2 df,   p=0.4045

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_FullGrp2NoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ pData(lumA)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()


# group3 - proliferation (TYMS,SQLE)
all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ factor(lumA$d_score_prol))

summary(all.surv.death)
#factor(lumA$d_score_prol)0 0.3820    1.4653   0.6715 0.569    0.569
#factor(lumA$d_score_prol)2 0.9692    2.6359   0.7314 1.325    0.185
#           				 exp(coef) exp(-coef) lower .95 upper .95
#factor(lumA$d_score_prol)0     1.465     0.6825    0.3929     5.464
#factor(lumA$d_score_prol)2     2.636     0.3794    0.6286    11.053
#Score (logrank) test = 1.86  on 2 df,   p=0.3937
#  n= 152, number of events= 12 
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/UNC337NKI295_LumAGrpProlNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumA$S_Time)),as.numeric(as.character(lumA$S_Event))) ~ factor(pData(lumA)$d_score_prol)),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)
dev.off()

