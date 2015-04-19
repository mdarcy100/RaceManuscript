#MD 9/2013: make images 350*350; remove all the text from images (do LumAB image, full image)
#MD - uodate images for MS:  go to line 117, then to 566
# caldas analysis on 6 genes.
#read in larger dataset
# this code will create an aggregate score of race/survival associated genes:
# CYRYBB2 - up is bad 1*exprs HR = 1.4
# ACOX2 - up is good -1* exprs HR = 0.65
# SQLE - up is bad 1*exprs HR = 1.98
# MUC1 - up is good-1*exprs HR = 0.65
# TYMS - up is bad 1*exprs HR = 2.67
# PSPH - up is bad 1*exprs HR = 1.65


## MD 1/2014.  evaluate survival, tumor characteristics with full set: ACOX2, 
# evaluate both discrete score and score.
# make sure that the associations hold for all the genes (CRYBB2, ACOX1 etc)
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase")
   
    biocLite("hgug4112a.db") # this is what we need to load hgug4112a.db
    biocLite("limma")
biocLite("gage")
biocLite("survplot")  # this package is apparently not available for 3.0.2
library(Biobase)
library(gdata) #
library(hgug4112a.db)
library(gplots)
library(limma)
library(gage)
library(survival)
library(muhaz)
library(survplot)

data(egSymb)

setwd("/Users/mdarcy/Desktop/MTroester/")
source("AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("AgePaper/AgeManuscript_062011/pcaFuncs.R")
source("AgePaper/AgeManuscript_062011/age_paper_scripts.R")

setwd("/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper")

load(file="caldas_race_eset.RData")

dim(caldas.race.eset) # this contains the proliferation genes
#Features  Samples 
#       14     1584 


medians.caldas <- esApply(caldas.race.eset,1,median)
caldas.exprs.scale <- t(scale(t(exprs(caldas.race.eset)),center=medians.caldas,scale=FALSE))

## Create a new eset
caldas.scale.eset <-caldas.race.eset
exprs(caldas.scale.eset) <- caldas.exprs.scale
dim(caldas.scale.eset)
##################################################################
# quick check of the tumor characteristics available
##################################################################
colnames(pData(caldas.scale.eset))
# [1] "Train_Test"                 "Index"                     
# [3] "Name"                       "Order"                     
# [5] "PAM50"                      "grade"                     
# [7] "size"                       "Node"                      
# [9] "ER_IHC_status"              "Treatment"                 
#[11] "age_at_diagnosis"           "RFSE10" 

table(caldas.scale.eset$PAM50)
# Basal   Her2   LumA   LumB Normal 
#   339    223    401    397    224 
table(caldas.scale.eset$grade)
#  1   2   3 
# 126 618 840 
table(caldas.scale.eset$Node)
#  0   1   2 
# 708 585 291 

table(caldas.scale.eset$size)
# these are strange values - need to look up metric
hist(as.numeric(caldas.scale.eset$size))
mean(as.numeric(as.character(caldas.scale.eset$size)))
#[1] 26.62723

caldas.scale.eset$size_group<-ifelse(as.numeric(as.character(caldas.scale.eset$size))>=20,1,0)
table(caldas.scale.eset$size_group)
#   0    1 
# 481 1103 

mean(as.numeric(as.character(caldas.scale.eset$age_at_diagnosis)))
#[1] 60.90307

# create a few more variables for age and size
age_cat_c<-cut(as.numeric(as.character(caldas.scale.eset$age_at_diagnosis)),breaks=c(0,30,40,50,60,70,80,100),right=FALSE) #  
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
prolif_score<-apply(exprs(caldas.scale.eset)[c(5,7:14),],2,sum)
pData(caldas.scale.eset)$prolif_score<-prolif_score
colnames(pData(caldas.scale.eset))
rownames(exprs(caldas.race.eset))

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

colnames(pData(caldas.scale.eset))
#[43] "f.age.c"                    "prolif_score"              
#[45] "ACOX2"                      "CRYBB2"                    
#[47] "MUC1"                       "SQLE"                      
#[49] "PSPH"                       "TYMS" 
scale_gene <- matrix(c(-1,1,-1,1,1,1), ncol = 6)
d_score<-scale_gene%*%t(pData(caldas.scale.eset)[,45:50])

quantile(d_score)
#  0%  25%  50%  75% 100% 
#  -6   -2    0    2    6 
#<=4, >=3 are good cutpoints


pData(caldas.scale.eset)$d_score<-t(d_score)

#get a continuous measure of 'badness'
rownames(exprs(caldas.race.eset))
#[1] "CRYBB2" "MUC1"   "PSPH"   "SQLE"   "TYMS"   "ACOX2" 
# we need 1,-1,1,1,1,-1
scale_gene <- matrix(c(1,-1,1,1,1,-1), ncol = 6)

# curiosity - looking at distribution of continuous gene expression
score<-scale_gene%*%exprs(caldas.scale.eset)
quantile(score)
#         0%         25%         50%         75%        100% 
#-23.8620  -5.4830   0.2165   7.0050  28.0990 
pData(caldas.scale.eset)$score<-t(score)


#### NOTE: these images weren't use in manuscript

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Score.tiff',width=400,height=400,res=100,antialias="none")
hist(d_score)
dev.off()

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_ContinuousScore.tiff',width=400,height=400,res=100,antialias="none")
hist(score)
dev.off()

#### END NOTE

###################### continuous score
#### NOTE: these images weren't use in manuscript, checking consistency of results
aov(pData(caldas.scale.eset)$score~pData(caldas.scale.eset)$PAM50)
fit_full<- aov(score~factor(PAM50),data=pData(caldas.scale.eset))


summary(aov(score~factor(PAM50),data=pData(caldas.scale.eset)))
summary.lm(fit_full)

#does the mean score vary by subtype
aov(pData(caldas.scale.eset)$score~pData(caldas.scale.eset)$PAM50)
fit_full<- aov(score~factor(PAM50),data=pData(caldas.scale.eset))


summary(aov(d_score~factor(PAM50),data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           8.9996     0.3617   24.88   <2e-16 ***
#factor(PAM50)Her2    -7.0581     0.5742  -12.29   <2e-16 ***
#factor(PAM50)LumA   -14.4821     0.4913  -29.47   <2e-16 ***
#factor(PAM50)LumB    -8.0431     0.4925  -16.33   <2e-16 ***
#factor(PAM50)Normal -10.8945     0.5734  -19.00   <2e-16 ***

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_SubtypeScoreContinuous.tiff',width=350,height=350,res=100,antialias="none")
boxplot(score~factor(PAM50),data=pData(caldas.scale.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="continuous sum score",xlab="subtype", main="Cont. sum score",sub="F-stat= 228.6 on 4 & 1579 DF, p-value: < 2.2e-16", cex=0.7)
abline(h = 0, col = "black", lty="dotted")
dev.off()

#### END NOTE


###################### cdiscrete score
# first evaluate the association between 'points' and proliferation score.
table(pData(caldas.scale.eset)$d_score)
cor(pData(caldas.scale.eset)$prolif_score,pData(caldas.scale.eset)$d_score,use="complete.obs")
#          [,1]
#[1,] 0.5932023
points_fit<-lm(pData(caldas.scale.eset)$prolif_score~pData(caldas.scale.eset)$d_score) 

summary(points_fit)
#pData(caldas.scale.eset)$d_score  2.44783    0.08352  29.308   <2e-16 ***
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/Caldas1584_ProliferationPointsCorrelation.tiff',width=350,height=350,res=100,antialias="none")
plot(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$prolif_score)
abline(points_fit)
dev.off()


#does the mean score vary by subtype, by proliferation score
# points
aov(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$PAM50)
fit_full<- aov(d_score~factor(PAM50),data=pData(caldas.scale.eset))
summary(aov(d_score~factor(PAM50),data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           2.9086     0.1428  20.365   <2e-16 ***
#factor(PAM50)Her2    -2.1462     0.2267  -9.466   <2e-16 ***
#factor(PAM50)LumA    -5.2677     0.1940 -27.150   <2e-16 ***
#factor(PAM50)LumB    -2.6416     0.1945 -13.584   <2e-16 ***
#factor(PAM50)Normal  -4.2836     0.2264 -18.918   <2e-16 ***

#"F-stat= 205.4 on 4 & 1579 DF, p-value: < 2.2e-16"
#should make a variable for length/width
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_SubtypeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(PAM50),data=pData(caldas.scale.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$PAM50)

# now do proliferation score
aov(pData(caldas.scale.eset)$prolif_score~pData(caldas.scale.eset)$PAM50)
fit_full<- aov(prolif_score~factor(PAM50),data=pData(caldas.scale.eset))
summary(aov(prolif_score~factor(PAM50),data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          12.7061     0.4587  27.701  < 2e-16 ***
#factor(PAM50)Her2    -5.9428     0.7282  -8.161  6.7e-16 ***
#factor(PAM50)LumA   -24.2126     0.6231 -38.858  < 2e-16 ***
#factor(PAM50)LumB    -6.9407     0.6245 -11.113  < 2e-16 ***
#factor(PAM50)Normal -25.9008     0.7272 -35.617  < 2e-16 ***

#F-statistic:   596 on 4 and 1579 DF,  p-value: < 2.2e-16
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/ExtendedAnalysisImages/CALDAS1584_SubtypeProliferationScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(prolif_score~factor(PAM50),data=pData(caldas.scale.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="proliferation score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(caldas.scale.eset)$prolif_score,pData(caldas.scale.eset)$PAM50)
#group: Basal
#  var   n  mean   sd median trimmed  mad    min   max range  skew kurtosis   se
#1   1 339 12.71 8.55  13.69   13.06 8.94 -14.63 32.71 47.34 -0.39    -0.11 0.46
#------------------------------------------------------------ 
#group: Her2
#  var   n mean   sd median trimmed  mad    min  max range skew kurtosis   se
#1   1 223 6.76 7.65   6.53    6.66 7.41 -11.82 28.1 39.91 0.16    -0.25 0.51
#------------------------------------------------------------ 
#group: LumA
#  var   n   mean   sd median trimmed  mad    min  max range  skew kurtosis   se
#1   1 401 -11.51 8.54 -10.82  -11.18 8.31 -38.37 8.55 46.92 -0.37    -0.01 0.43
#------------------------------------------------------------ 
#group: LumB
#  var   n mean   sd median trimmed mad    min   max range skew kurtosis  se
#1   1 397 5.77 7.93   5.05    5.42 8.4 -10.47 26.93  37.4 0.35    -0.51 0.4
#------------------------------------------------------------ 
#group: Normal
#  var   n   mean   sd median trimmed  mad    min   max range skew kurtosis   se
#1   1 224 -13.19 9.67 -12.93   -13.2 10.3 -35.33 12.28  47.6 0.03    -0.46 0.65


########does the mean score vary by age
aov(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$f.age.c)
fit_full<- aov(d_score~f.age.c,data=pData(caldas.scale.eset))


summary(aov(d_score~f.age.c,data=pData(caldas.scale.eset)))
summary.lm(fit_full)
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       2.3636     0.9574   2.469 0.013663 * 
#f.age.c[30,40)   -0.4535     1.0149  -0.447 0.655018  
#f.age.c[40,50)   -2.1126     0.9792  -2.157 0.031119 *  
#f.age.c[50,60)   -2.1853     0.9715  -2.249 0.024632 *  
#f.age.c[60,70)   -2.3238     0.9690  -2.398 0.016592 *  
#f.age.c[70,80)   -3.2738     0.9730  -3.365 0.000785 ***
#f.age.c[80,100)  -2.6782     1.0149  -2.639 0.008396 ** 

#F-statistic: 11.56 on 6 and 1577 DF,  p-value: 1.077e-12
#should make a variable for length/width
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_AgeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~f.age.c,data=pData(caldas.scale.eset),names=c("<30","30-39","40-49","50-59","60-69","70-79","80+"), col=c("grey"),ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$f.age.c) #just one grouping variable 
#INDICES: [0,30)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 11 2.36 3.07      4    2.44 2.97  -2   6     8 -0.42    -1.51 0.93
#-------------------------------------------------------------------------------------------- 
#INDICES: [30,40)
#   var  n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 89 1.91 3.09      2    2.11 2.97  -6   6    12 -0.55    -0.52 0.33
#-------------------------------------------------------------------------------------------- 
#INDICES: [40,50)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 239 0.25 3.23      0    0.24 2.97  -6   6    12    0    -0.87 0.21
#-------------------------------------------------------------------------------------------- 
#INDICES: [50,60)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 370 0.18 3.23      0    0.14 2.97  -6   6    12 0.06    -0.82 0.17
#-------------------------------------------------------------------------------------------- 
#INDICES: [60,70)
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 452 0.04 3.26      0    0.03 2.97  -6   6    12 0.01    -0.86 0.15
#-------------------------------------------------------------------------------------------- 
#INDICES: [70,80)
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 334 -0.91 2.98     -2   -0.99 2.97  -6   6    12 0.21     -0.7 0.16
#-------------------------------------------------------------------------------------------- 
#INDICES: [80,100)
#   var  n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 89 -0.31 3.19      0   -0.38 2.97  -6   6    12 0.24    -0.76 0.34


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

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_SizeScore.tiff',width=350,height=350,res=100,antialias="none")
#plot(as.numeric(as.character(caldas.scale.eset$size))/10,pData(caldas.scale.eset)$d_score)
boxplot(d_score~size_group,data=pData(caldas.scale.eset),names=c("< 2cm",">= 2cm"), col=c("grey"),ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

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

tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_GradeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~grade,data=pData(caldas.scale.eset), col=c("grey"),names=c("1","2","3"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()

describeBy(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$grade) #just one grouping variable 
#INDICES: 1
#   var   n  mean  sd median trimmed  mad min max range skew kurtosis   se
#V1   1 126 -2.38 2.6     -2   -2.57 2.97  -6   6    12 0.63     0.14 0.23
#INDICES: 2
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 618 -1.06 2.85     -2   -1.18 2.97  -6   6    12 0.29    -0.49 0.11
#-------------------------------------------------------------------------------------------- 
#INDICES: 3
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 840 1.15 3.13      2    1.26 2.97  -6   6    12 -0.27    -0.71 0.11
#####################################################
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
describeBy(pData(caldas.scale.eset)$d_score,pData(caldas.scale.eset)$Node) #just one grouping variable 
#INDICES: 0
#   var   n  mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 708 -0.24 3.26      0    -0.3 2.97  -6   6    12 0.16    -0.84 0.12
#-------------------------------------------------------------------------------------------- 
#INDICES: 1
#   var   n mean   sd median trimmed  mad min max range skew kurtosis   se
#V1   1 585 0.17 3.29      0    0.15 2.97  -6   6    12 0.01    -0.86 0.14
#-------------------------------------------------------------------------------------------- 
#INDICES: 2
#   var   n mean   sd median trimmed  mad min max range  skew kurtosis   se
#V1   1 291 0.27 3.06      0    0.33 2.97  -6   6    12 -0.09    -0.73 0.18
tiff('/Users/mdarcy/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_NodeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~Node,data=pData(caldas.scale.eset), col=c("grey"),names=c("0","1","2"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


table(caldas.scale.eset$grade)
#  1   2   3 
# 126 618 840 
table(caldas.scale.eset$Node)
#  0   1   2 
# 708 585 291 






##########################################################################
#Survival analyses
######################################################## 
#Discrete (points) continuous (cox proportional) model


all.surv.death <- coxph(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$d_score)
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


factor(pData(caldas.scale.eset)$d_score)
# NOTE: we don't use this images
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_All_DiscreteContModel.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$d_score),
     lty = 1:3,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, N=1584")
     
     legend(0, 0.4, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(caldas.scale.eset$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(caldas.scale.eset$d_score)==0)),"; HR=1.42 (1.39,1.46)",sep=""),
     paste("Score = 6; N = ",length(which(factor(caldas.scale.eset$d_score)==6)),"; HR=2.03 (1.98,2.08)",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(4,0.1,"Discrete/Cont Model; exp(B)=1.06; p=3.7e-6",cex=0.7)
dev.off()
#END NOTE

######################################################## 
# NOTE: we don't use this model Continuous continuous model


all.surv.death <- coxph(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$score,surv=TRUE)
summary(all.surv.death)
#pData(caldas.scale.eset)$score 0.017161  1.017309 0.004868 3.525 0.000423 ***
#pData(caldas.scale.eset)$score     1.017      0.983     1.008     1.027


min(score)
################################### 
# testing plotting: need to work on these continuous plots
plot(survfit(fit.1, data.frame(Prison = 0, Dose = 62, Clinic = 1),conf.type="none")) 
lines(survfit(fit.1, data.frame(Prison = 0, Dose = 62, Clinic = 2)),col=2) 
plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$score,pData(caldas.scale.eset)$score=-23.862)
     lty = 1,col=c("blue"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, N=1584")
     
     legend(0, 0.4, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(caldas.scale.eset$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(caldas.scale.eset$d_score)==0)),"; HR=1.42 (1.39,1.46)",sep=""),
     paste("Score = 6; N = ",length(which(factor(caldas.scale.eset$d_score)==6)),"; HR=2.03 (1.98,2.08)",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     
     plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$score),
     lty = 1:3,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, N=1584")
     
     legend(0, 0.4, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(caldas.scale.eset$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(caldas.scale.eset$d_score)==0)),"; HR=1.42 (1.39,1.46)",sep=""),
     paste("Score = 6; N = ",length(which(factor(caldas.scale.eset$d_score)==6)),"; HR=2.03 (1.98,2.08)",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(4,0.1,"Discrete/Cont Model; exp(B)=1.06; p=3.7e-6",cex=0.7)


#given this spread
##         0%         25%         50%         75%        100% 
#-23.8620  -5.4830   0.2165   7.0050  28.0990 
#Calculate HR and 95% CI for Upper and lower values
#HR 6 points/-6 points = exp(12*0.0589) = [1] 2.027493
##
#Point estimate for high score:
 exp((28.09+23.86)*0.017161)
#[1] 2.438819


#Point estimate for middle score:
 exp((.2165+23.862)*0.017161)
#[1] 1.511664



L95 = exp((28.09+23.86)*0.017161-1.96*0.004868)
#[1] 2.41566
U95= exp((28.09+23.86)*0.017161+1.96*0.004868)
#[1] 2.4622
#HR 0 points/-6 points = exp(6*0.0589) =[1] 1.511664

L95 = exp((.2165+23.862)*0.017161-1.96*0.004868)
#[1] 1.49731
U95= exp((.2165+23.862)*0.017161+1.96*0.004868)
#[1] 1.526156


factor(pData(caldas.scale.eset)$d_score)



tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_All_DiscreteContModel.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$score),
     lty = 1:3,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, N=1584")
     
     legend(0, 0.4, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(caldas.scale.eset$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(caldas.scale.eset$d_score)==0)),"; HR=1.42 (1.39,1.46)",sep=""),
     paste("Score = 6; N = ",length(which(factor(caldas.scale.eset$d_score)==6)),"; HR=2.03 (1.98,2.08)",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(4,0.1,"Discrete/Cont Model; exp(B)=1.06; p=3.7e-6",cex=0.7)
dev.off()







#######################################################
# now look at survival by subtype
######################### #########################  #########################  
######################### ######################### 
#survival by score among both luminal A/B
################################ - we don't use this image
lumAB<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("LumA","LumB"))]
dim(lumAB)

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_LumABScore.tiff',width=400,height=400,res=100,antialias="none")
hist(lumAB$d_score)
dev.off()
quantile(lumAB$d_score)
#  0%  25%  50%  75% 100% 
#  -6   -4   -2    0    6 



all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ pData(lumAB)$d_score)
summary(all.surv.death)
#pData(lumAB)$d_score 0.03908   1.03985  0.02001 1.953   0.0508 .
#pData(lumAB)$d_score      1.04     0.9617    0.9999     1.081
#Score (logrank) test = 3.82  on 1 df,   p=0.05063

#also calculate the 95% CI 
#HR 6 points/-6 points = exp(12*0.03908) = [1] 1.60
L95 = exp(12*0.03908-1.96*0.02001)
#[1] 1.536859
U95= exp(12*0.03908+1.96*0.02001)
#[1] 1.662262
#HR 0 points/-6 points = exp(6*0.03908) = [1] 1.26


L95 = exp(6*0.03908-1.96*0.02001)
#[1] 1.215627
U95= exp(6*0.03908+1.96*0.02001)
#[1] 1.31482




tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_LumAB_ContinuousModel.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ pData(lumAB)$d_score),
     lty = 1:3,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, LumA, N=798")
     
     legend(0, 0.6, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(lumAB$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(lumAB$d_score)==0)),"; HR=1.26",sep=""),
     paste("Score = 6; N = ",length(which(factor(lumAB$d_score)==6)),"; HR=1.60",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(4,0.2,"Cont. Model, exp(B)=1.04, 95%CI=(1.00,1.08)",cex=0.7)
dev.off()



all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ pData(lumAB)$score)
summary(all.surv.death)


################################################################################################################################## - we don't use this image
lumANorm<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("LumA","Normal"))]
dim(lumANorm)
#Features  Samples 
#       6      625 
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_LumANormScore.tiff',width=400,height=400,res=100,antialias="none")
hist(lumANorm$d_score)
dev.off()
quantile(lumANorm$d_score)
#  0%  25%  50%  75% 100% 
#  -6   -4   -2    0    6 



all.surv.death <- coxph(Surv(as.numeric(as.character(lumANorm$RFS10yr)),as.numeric(as.character(lumANorm$RFSE10))) ~ pData(lumANorm)$d_score)
summary(all.surv.death)
#pData(lumANorm)$d_score 0.04970   1.05095  0.03022 1.644      0.1
#pData(lumANorm)$d_score     1.051     0.9515    0.9905     1.115
#Score (logrank) test = 2.71  on 1 df,   p=0.09992

#also calculate the 95% CI 
#HR 6 points/-6 points = exp(12*0.04970) = [1] 1.82
L95 = exp(12*0.04970-1.96*0.03022)
#[1] 1.71
U95= exp(12*0.04970+1.96*0.03022)
#[1] 1.93
#HR 0 points/-6 points = exp(6*0.04970) = [1] 1.35


L95 = exp(6*0.04970-1.96*0.03022)
#[1] 1.27
U95= exp(6*0.04970+1.96*0.03022)
#[1] 1.43




tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_LumANorm_ContinuousModel.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumANorm$RFS10yr)),as.numeric(as.character(lumANorm$RFSE10))) ~ pData(lumANorm)$d_score),
     lty = 1:6,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, LumA/Normal, N=625")
     
     legend(0, 0.6, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(lumANorm$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(lumANorm$d_score)==0)),"; HR=1.35",sep=""),
     paste("Score = 6; N = ",length(which(factor(lumANorm$d_score)==6)),"; HR=1.82",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(4,0.2,"Cont. Model, exp(B)=1.05, 95%CI=(0.99,1.12)",cex=0.7)
dev.off()



all.surv.death <- coxph(Surv(as.numeric(as.character(lumANorm$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ pData(lumAB)$score)
summary(all.surv.death)
##################################################################################################################################






lumA<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("LumA"))]

table(pData(lumA)$f.dscore.c)
#[-7,-2)  [-2,3)   [3,7) 
#    177     210      14 

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_LumAScore.tiff',width=400,height=400,res=100,antialias="none")
hist(lumA$d_score)
dev.off()
quantile(lumA$d_score)
#  0%  25%  50%  75% 100% 
#  -6   -4   -2    0    6 

lumA_dscore_cat_c<-cut(lumA$d_score,breaks=c(-7,-4,-2,0,3,7),right=FALSE) # corresponds to quantiles
f.lumA.dscore.c<-factor(lumA_dscore_cat_c)
table(f.lumA.dscore.c)
#[-7,-4) [-4,-2)  [-2,0)   [0,3)   [3,7) 
#     61     116     117      93      14 



all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ lumA$d_score)

summary(all.surv.death)
#pData(caldas.scale.eset)$f.dscore.c[-2,3)     1.658     0.6033     1.305     2.105
#pData(caldas.scale.eset)$f.dscore.c[3,7)      1.875     0.5333     1.427     2.463

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_All.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$f.dscore.c),
     lty = 1:3,col=c("darkblue","black","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, N=1584")
     
     
legend(0, 0.35, 
legend = c(paste("Bottom 25% (ref), N = ",length(which(factor(caldas.scale.eset$f.dscore.c)=='[-7,-2)')),sep=""), 
paste("Middle 50%, N = ",length(which(factor(caldas.scale.eset$f.dscore.c)=='[-2,3)')),"; HR=", round(summary(all.surv.death)$coef[3],2)," P=", round(summary(all.surv.death)$coef[9],2),sep=""),
paste("Top 25%, N = ",length(which(factor(caldas.scale.eset$f.dscore.c)=='[3,7)')),"; HR=", round(summary(all.surv.death)$coef[4],2)," P=", round(summary(all.surv.death)$coef[10],2),sep="")),

 fill = c("darkblue","black","darkred"),lty = c(1,2,3), bty = "n",cex=0.7)     

dev.off()




########################################################################################

all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ pData(lumA)$d_score)
summary(all.surv.death)
#pData(lumA)$d_score 0.04708   1.04820  0.03658 1.287    0.198
#pData(lumA)$d_score     1.048      0.954    0.9757     1.126


#also calculate the 95% CI 
#HR 6 points/-6 points = exp(12*0.04708) = [1] 1.76
L95 = exp(12*0.04708-1.96*0.03658)
#[1] 1.637651
U95= exp(12*0.04708+1.96*0.03658)
##[[1] 1.890151
#HR 0 points/-6 points = exp(6*0.04708) = [1] 1.33


L95 = exp(6*0.04708-1.96*0.03658)
#[1] 1.234645
U95= exp(6*0.04708+1.96*0.03658)
#[1] 1.425007




tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/SurvivalImages/SumScore/CALDAS1584_Survival_LumA_ContinuousModel.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ pData(lumA)$d_score),
     lty = 1:3,col=c("blue","darkblue","darkblue","black","darkred","darkred","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="CALDAS, LumA, N=401")
     
     legend(0, 0.6, legend = c(paste("Score = -6 (ref), N = ",length(which(factor(lumA$d_score)==-6)),sep=""), 
     paste("Score = 0; N = ",length(which(factor(lumA$d_score)==0)),"; HR=1.33",sep=""),
     paste("Score = 6; N = ",length(which(factor(lumA$d_score)==6)),"; HR=1.76",sep="")),
     fill = c("blue","black","red"),lty = c(1,2,3), bty = "n",cex=0.7)  
     
     text(7,0.2,"Continuous Model, exp(B)=1.05, p = 0.198",cex=0.7)
dev.off()



all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ pData(lumA)$score)
summary(all.surv.death)


########################################################################################
Basal<-caldas.scale.eset[,which(caldas.scale.eset$PAM50%in%c("Basal"))]

table(pData(Basal)$f.dscore.c)

all.surv.death <- coxph(Surv(as.numeric(as.character(Basal$RFS10yr)),as.numeric(as.character(Basal$RFSE10))) ~ pData(Basal)$d_score)
summary(all.surv.death)
#pData(Basal)$d_score -0.01466   0.98545  0.03412 -0.43    0.667
#pData(Basal)$d_score    0.9854      1.015    0.9217     1.054


#also calculate the 95% CI 
#HR 6 points/-6 points = exp(12*(-0.01466)) = #[1] 0.8386851
L95 = exp(12*-0.01466-1.96*0.03412)
#[[1] 0.7844321
U95= exp(12*-0.01466+1.96*0.03412)
#[1] 0.8966902
#HR 0 points/-6 points = exp(6*0.04708) = [1] 1.33


L95 =  exp(6*-0.01466-1.96*0.03412)
#[1] 1.234645
U95= exp(6*-0.01466-1.96*0.03412)
#[1] 1.425007






caldas.large[race.idx,2]
#[1] "CRYBB2" "MUC1"   "PSPH"   "SQLE"   "TYMS"   "ACOX2" 




##########################################################################################################
#categorized survival images
# MD - updated 8/21/2013.  put into categories not continuous points
# NOTE: images for manuscript 
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
##########################################################################################################
all.surv.death <- coxph(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ lumA$f.dscore.c)

summary(all.surv.death)
#pData(caldas.scale.eset)$f.dscore.c[-2,3)1.490     0.6713    1.0085     2.200
#pData(caldas.scale.eset)$f.dscore.c[3,7) 1.443     0.6930    0.5167     4.031

#NOTE: we don't use this image;  not enough power in top group
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/Caldas_LumA_NoText.tiff',width=400,height=400,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumA$RFS10yr)),as.numeric(as.character(lumA$RFSE10))) ~ pData(lumA)$f.dscore.c),
     lty = 1:3,col=c("darkblue","black","darkred"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="survival by total score",sub="LumA, N=401")
     
     
legend(0, 0.35, legend=c("low","med","high"),
#legend = c(paste("Lowest Risk (ref), N = ",length(which(factor(lumA$f.dscore.c)=='[-7,-2)')),sep=""), 
#paste("Middle Risk, N = ",length(which(factor(lumA$f.dscore.c)=='[-2,3)')),"; HR=", round(summary(all.surv.death)$coef[3],2)," P=", round(summary(all.surv.death)$coef[9],2),sep=""),
#paste("Top Risk, N = ",length(which(factor(lumA$f.dscore.c)=='[3,7)')),"; HR=", round(summary(all.surv.death)$coef[4],2)," P=", round(summary(all.surv.death)$coef[10],2),sep="")),

 fill = c("darkblue","black","darkred"),lty = c(1,2,3), bty = "n",cex=0.7)     

dev.off()


# now do lumA/B
all.surv.death <- coxph(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ lumAB$f.dscore.c)

summary(all.surv.death)
#pData(caldas.scale.eset)$f.dscore.c[-2,3)     1.593     0.6277    1.1907     2.132
#pData(caldas.scale.eset)$f.dscore.c[3,7)      1.501     0.6663    0.9644     2.336
#Score (logrank) test = 10.09  on 2 df,   p=0.006433
#CALDAS, Luminal, N=798
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/Caldas_LumABNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(lumAB$RFS10yr)),as.numeric(as.character(lumAB$RFSE10))) ~ pData(lumAB)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)

dev.off()



# survival in all subtype
all.surv.death <- coxph(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ caldas.scale.eset$f.dscore.c)

summary(all.surv.death)
#caldas.scale.eset$f.dscore.c[-2,3) 0.5054    1.6577   0.1220 4.143 3.43e-05 ***
#caldas.scale.eset$f.dscore.c[3,7)  0.6286    1.8749   0.1392 4.516 6.29e-06 ***

#pData(caldas.scale.eset)$f.dscore.c[-2,3)     1.658     0.6033     1.305     2.105
#pData(caldas.scale.eset)$f.dscore.c[3,7)      1.875     0.5333     1.427     2.463

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/Caldas_FullNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ pData(caldas.scale.eset)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)     
 dev.off()



######################Figure 5 images
aov(pData(caldas.scale.eset)$d_score~pData(caldas.scale.eset)$PAM50)
fit_full<- aov(d_score~factor(PAM50),data=pData(caldas.scale.eset))


summary(aov(d_score~factor(PAM50),data=pData(caldas.scale.eset)))
summary.lm(fit_full)

#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           2.9086     0.1428  20.365   <2e-16 ***
#factor(PAM50)Her2    -2.1462     0.2267  -9.466   <2e-16 ***
#factor(PAM50)LumA    -5.2677     0.1940 -27.150   <2e-16 ***
#factor(PAM50)LumB    -2.6416     0.1945 -13.584   <2e-16 ***
#factor(PAM50)Normal  -4.2836     0.2264 -18.918   <2e-16 ***

#"F-stat= 205.4 on 4 & 1579 DF, p-value: < 2.2e-16"
tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_SubtypeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(PAM50),data=pData(caldas.scale.eset), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


########################################################################################################## AGE STRATIFIED ANALYSES.  CALDAS is a lot older than NKI (45) + UNC337 (57) 
young<-caldas.scale.eset[,which(as.numeric(as.character(caldas.scale.eset$age_at_diagnosis))<60)]
old<-caldas.scale.eset[,which(as.numeric(as.character(caldas.scale.eset$age_at_diagnosis))>=60)]
#look in people who are <= 60; this has a different age distribution than the other data
mean(as.numeric(as.character(pData(caldas.scale.eset)$age_at_diagnosis)))
#60.9
#max(as.numeric(as.character(pData(caldas.scale.eset)$age_at_diagnosis)))
#[1] 96.29

#stratify by age and look at those >= 60 and those < 60
 mean(as.numeric(as.character(old$age_at_diagnosis)))
#[1] 70.47789
 mean(as.numeric(as.character(young$age_at_diagnosis)))
#[1] 49.08647


# survival in all subtype in all young

all.surv.death <- coxph(Surv(as.numeric(as.character(young$RFS10yr)),as.numeric(as.character(young$RFSE10))) ~ young$f.dscore.c)

summary(all.surv.death)
#young$f.dscore.c[-2,3) 1.1764    3.2427   0.3023 3.891 9.97e-05 ***
#young$f.dscore.c[3,7)  1.4267    4.1649   0.3139 4.545 5.50e-06 ***
#young$f.dscore.c[-2,3)     3.243     0.3084     1.793     5.864
#young$f.dscore.c[3,7)      4.165     0.2401     2.251     7.706

# survival in all subtype in all young


tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/CaldasYOUNG_FullNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(young$RFS10yr)),as.numeric(as.character(young$RFSE10))) ~ pData(young)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)     
 dev.off()

# survival in all subtype in all old

all.surv.death <- coxph(Surv(as.numeric(as.character(old$RFS10yr)),as.numeric(as.character(old$RFSE10))) ~ old$f.dscore.c)

summary(all.surv.death)
#old$f.dscore.c[-2,3) 0.3812    1.4640   0.1359 2.805  0.00503 **
#old$f.dscore.c[3,7)  0.4680    1.5968   0.1680 2.786  0.00533 **
#old$f.dscore.c[-2,3)     1.464     0.6831     1.122     1.911
#old$f.dscore.c[3,7)      1.597     0.6263     1.149     2.219

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/CaldasOLD_FullNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(old$RFS10yr)),as.numeric(as.character(old$RFSE10))) ~ pData(old)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)     
 dev.off()


##################################################################
# Now analyze only in Luminal A/B in young or old
young_lumAB<-young[,which(young$PAM50%in%c("LumA","LumB"))]
old_lumAB<-old[,which(old$PAM50%in%c("LumA","LumB"))]



all.surv.death <- coxph(Surv(as.numeric(as.character(young_lumAB$RFS10yr)),as.numeric(as.character(young_lumAB$RFSE10))) ~ young_lumAB$f.dscore.c)

#  n= 254, number of events= 47 
summary(all.surv.death)
#young_lumAB$f.dscore.c[-2,3) 1.0637    2.8970   0.4403 2.416   0.0157 *
#young_lumAB$f.dscore.c[3,7)  0.8916    2.4391   0.6456 1.381   0.1673  
#young_lumAB$f.dscore.c[-2,3)     2.897     0.3452    1.2223     6.866
#young_lumAB$f.dscore.c[3,7)      2.439     0.4100    0.6881     8.646

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/CaldasYOUNG_LumABNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(young_lumAB$RFS10yr)),as.numeric(as.character(young_lumAB$RFSE10))) ~ pData(young_lumAB)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)     
 dev.off()

######## ########

all.surv.death <- coxph(Surv(as.numeric(as.character(old_lumAB$RFS10yr)),as.numeric(as.character(old_lumAB$RFSE10))) ~ old_lumAB$f.dscore.c)

#  n= 544, number of events= 221 
summary(all.surv.death)
#old_lumAB$f.dscore.c[-2,3) 0.4134    1.5119   0.1593 2.596  0.00944 **
#old_lumAB$f.dscore.c[3,7)  0.2921    1.3393   0.2414 1.210  0.22614   
#old_lumAB$f.dscore.c[-2,3)     1.512     0.6614    1.1065     2.066
#old_lumAB$f.dscore.c[3,7)      1.339     0.7467    0.8345     2.149

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure4/350Images/CaldasOLD_LumABNoText.tiff',width=350,height=350,res=100,antialias="none")
plot(survfit(Surv(as.numeric(as.character(old_lumAB$RFS10yr)),as.numeric(as.character(old_lumAB$RFSE10))) ~ pData(old_lumAB)$f.dscore.c),
     lty = 1:3,col=c("blue","black","red"), mark.time = TRUE)     
 dev.off()


###################### look at distribution of points in young, old
aov(pData(young)$d_score~pData(young)$PAM50)
fit_full<- aov(d_score~factor(PAM50),data=pData(young))


summary(aov(d_score~factor(PAM50),data=pData(young)))
summary.lm(fit_full)

#does the mean score vary by subtype
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           3.1080     0.1753  17.729  < 2e-16 ***
#factor(PAM50)Her2    -2.0729     0.2969  -6.982 6.76e-12 ***
#factor(PAM50)LumA    -5.2366     0.2784 -18.811  < 2e-16 ***
#factor(PAM50)LumB    -2.8975     0.2969  -9.759  < 2e-16 ***
#factor(PAM50)Normal  -4.5455     0.2861 -15.886  < 2e-16 ***

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_YOUNG_SubtypeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(PAM50),data=pData(young), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


###################### look at distribution of points in old
aov(pData(old)$d_score~pData(old)$PAM50)
fit_full<- aov(d_score~factor(PAM50),data=pData(old))


summary(aov(d_score~factor(PAM50),data=pData(old)))
summary.lm(fit_full)

#does the mean score vary by subtype
#                    Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           2.5714     0.2389  10.763  < 2e-16 ***
#factor(PAM50)Her2    -2.0944     0.3508  -5.970 3.45e-09 ***
#factor(PAM50)LumA    -5.0542     0.2909 -17.373  < 2e-16 ***
#factor(PAM50)LumB    -2.2817     0.2872  -7.944 6.04e-15 ***
#factor(PAM50)Normal  -3.8631     0.3633 -10.633  < 2e-16 ***

tiff('/Users/mdarcy100/Desktop/MTroester/AncestryMicroarray/Paper/Figure5/350Images/CALDAS1584_OLD_SubtypeScore.tiff',width=350,height=350,res=100,antialias="none")
boxplot(d_score~factor(PAM50),data=pData(old), col=c("red","pink","lightblue","darkblue","green"), ylab="risk score")
abline(h = 0, col = "black", lty="dotted")
dev.off()


#adjustment for age ("age_at_diagnosis"), node  "Node" , grade ("grade"), size (TsizeG)

full.adj.surv.death <- coxph(Surv(as.numeric(as.character(caldas.scale.eset$RFS10yr)),as.numeric(as.character(caldas.scale.eset$RFSE10))) ~ caldas.scale.eset$f.dscore.c+as.numeric(as.character(caldas.scale.eset$age_at_diagnosis))+factor(caldas.scale.eset$Node)+factor(caldas.scale.eset$grade)+factor(caldas.scale.eset$TsizeG))

summary(full.adj.surv.death)
                                                                 coef exp(coef) se(coef)     z Pr(>|z|)    
#caldas.scale.eset$f.dscore.c[-2,3)                           0.410146  1.507037 0.124912 3.283 0.001025 ** 
#caldas.scale.eset$f.dscore.c[3,7)                            0.511283  1.667430 0.146928 3.480 0.000502 ***
#as.numeric(as.character(caldas.scale.eset$age_at_diagnosis)) 0.024686  1.024993 0.003544 6.966 3.26e-12 ***
#factor(caldas.scale.eset$Node)1                              0.433959  1.543355 0.101148 4.290 1.78e-05 ***
#factor(caldas.scale.eset$Node)2                              1.037750  2.822858 0.110167 9.420  < 2e-16 ***
#factor(caldas.scale.eset$grade)2                             0.317957  1.374317 0.218032 1.458 0.144757    
#factor(caldas.scale.eset$grade)3                             0.630142  1.877877 0.216866 2.906 0.003665 ** 
#factor(caldas.scale.eset$TsizeG)1                            0.258277  1.294698 0.092727 2.785 0.005347 ** 
#factor(caldas.scale.eset$TsizeG)2                            0.617838  1.854914 0.216185 2.858 0.004264 ** 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#                                                             exp(coef) exp(-coef) lower .95 upper .95
#caldas.scale.eset$f.dscore.c[-2,3)                               1.507     0.6636    1.1798     1.925
#caldas.scale.eset$f.dscore.c[3,7)                                1.667     0.5997    1.2502     2.224
#as.numeric(as.character(caldas.scale.eset$age_at_diagnosis))     1.025     0.9756    1.0179     1.032
#factor(caldas.scale.eset$Node)1                                  1.543     0.6479    1.2658     1.882
#factor(caldas.scale.eset$Node)2                                  2.823     0.3543    2.2746     3.503
#factor(caldas.scale.eset$grade)2                                 1.374     0.7276    0.8964     2.107
#factor(caldas.scale.eset$grade)3                                 1.878     0.5325    1.2276     2.873
#factor(caldas.scale.eset$TsizeG)1                                1.295     0.7724    1.0795     1.553
#factor(caldas.scale.eset$TsizeG)2                                1.855     0.5391    1.2142     2.834





young.adj.surv.death <- coxph(Surv(as.numeric(as.character(young$RFS10yr)),as.numeric(as.character(young$RFSE10))) ~ young$f.dscore.c+as.numeric(as.character(young$age_at_diagnosis))+factor(young$Node)+factor(young$grade)+factor(young$TsizeG))

summary(young.adj.surv.death)

                                                      coef exp(coef)  se(coef)      z Pr(>|z|)    
#young$f.dscore.c[-2,3)                            0.787274  2.197399  0.307248  2.562  0.01040 *  
#young$f.dscore.c[3,7)                             0.727036  2.068939  0.326841  2.224  0.02612 *  
#as.numeric(as.character(young$age_at_diagnosis)) -0.007175  0.992851  0.008922 -0.804  0.42129    
#factor(young$Node)1                               0.319330  1.376205  0.180002  1.774  0.07606 .  
#factor(young$Node)2                               1.170602  3.223934  0.189849  6.166 7.01e-10 ***
#factor(young$grade)2                              0.750197  2.117418  0.524743  1.430  0.15282    
#factor(young$grade)3                              1.461868  4.314012  0.515647  2.835  0.00458 ** 
#factor(young$TsizeG)1                             0.097179  1.102057  0.153539  0.633  0.52678    
#factor(young$TsizeG)2                            -0.028849  0.971563  0.405937 -0.071  0.94334    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

                                                 exp(coef) exp(-coef) lower .95 upper .95
#young$f.dscore.c[-2,3)                              2.1974     0.4551    1.2033     4.013
#young$f.dscore.c[3,7)                               2.0689     0.4833    1.0903     3.926
#as.numeric(as.character(young$age_at_diagnosis))    0.9929     1.0072    0.9756     1.010
#factor(young$Node)1                                 1.3762     0.7266    0.9671     1.958
#factor(young$Node)2                                 3.2239     0.3102    2.2222     4.677
#factor(young$grade)2                                2.1174     0.4723    0.7571     5.922
#factor(young$grade)3                                4.3140     0.2318    1.5702    11.852
#factor(young$TsizeG)1                               1.1021     0.9074    0.8157     1.489
#factor(young$TsizeG)2                               0.9716     1.0293    0.4385     2.153

old.adj.surv.death <- coxph(Surv(as.numeric(as.character(old$RFS10yr)),as.numeric(as.character(old$RFSE10))) ~ old$f.dscore.c+as.numeric(as.character(old$age_at_diagnosis))+factor(old$Node)+factor(old$grade)+factor(old$TsizeG))

summary(old.adj.surv.death)


#                                                   coef exp(coef) se(coef)     z Pr(>|z|)    
#old$f.dscore.c[-2,3)                           0.345257  1.412353 0.140730 2.453 0.014154 *  
#old$f.dscore.c[3,7)                            0.445138  1.560706 0.175843 2.531 0.011359 *  
#as.numeric(as.character(old$age_at_diagnosis)) 0.052649  1.054060 0.007667 6.867 6.57e-12 ***
#factor(old$Node)1                              0.415482  1.515101 0.124252 3.344 0.000826 ***
#factor(old$Node)2                              0.954056  2.596219 0.137232 6.952 3.60e-12 ***
#factor(old$grade)2                             0.153902  1.166377 0.241657 0.637 0.524214    
#factor(old$grade)3                             0.334188  1.396806 0.244305 1.368 0.171338    
#factor(old$TsizeG)1                            0.269072  1.308750 0.116756 2.305 0.021191 *  
#factor(old$TsizeG)2                            0.871318  2.390060 0.256841 3.392 0.000693 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

                                               exp(coef) exp(-coef) lower .95 upper .95
#old$f.dscore.c[-2,3)                               1.412     0.7080    1.0719     1.861
#old$f.dscore.c[3,7)                                1.561     0.6407    1.1057     2.203
#as.numeric(as.character(old$age_at_diagnosis))     1.054     0.9487    1.0383     1.070
#factor(old$Node)1                                  1.515     0.6600    1.1876     1.933
#factor(old$Node)2                                  2.596     0.3852    1.9839     3.397
#factor(old$grade)2                                 1.166     0.8574    0.7263     1.873
#factor(old$grade)3                                 1.397     0.7159    0.8653     2.255
#factor(old$TsizeG)1                                1.309     0.7641    1.0411     1.645
#factor(old$TsizeG)2                                2.390     0.4184    1.4447     3.954

library(psych)
corr.test(t(exprs((caldas.scale.eset)))
#Correlation matrix 
#       CRYBB2  MUC1  PSPH  SQLE  TYMS ACOX2
#CRYBB2   1.00 -0.08  0.03  0.05  0.04 -0.04
#MUC1    -0.08  1.00 -0.10 -0.27 -0.35  0.27
#PSPH     0.03 -0.10  1.00  0.19  0.24 -0.12
#SQLE     0.05 -0.27  0.19  1.00  0.42 -0.18
#TYMS     0.04 -0.35  0.24  0.42  1.00 -0.18
#ACOX2   -0.04  0.27 -0.12 -0.18 -0.18  1.00

#Probability values (Entries above the diagonal are adjusted for multiple tests.) 
#       CRYBB2 MUC1 PSPH SQLE TYMS ACOX2
#CRYBB2   0.00 0.01 0.35 0.23 0.35  0.35
#MUC1     0.00 0.00 0.00 0.00 0.00  0.00
#PSPH     0.18 0.00 0.00 0.00 0.00  0.00
#SQLE     0.06 0.00 0.00 0.00 0.00  0.00
#TYMS     0.12 0.00 0.00 0.00 0.00  0.00
#ACOX2    0.13 0.00 0.00 0.00 0.00  0.00


