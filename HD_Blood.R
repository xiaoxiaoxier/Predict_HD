source("https://bioconductor.org/biocLite.R")
biocLite("NanoStringDiff")
library(NanoStringDiff)
library(xlsx)
rm(list=ls())
assay=read.xlsx("~/neurogen/predictHD/data/_for_dbGap_/PREDICT-HD.dbGap.20160608/miRNA.NanoString.20160608.xlsx",1,stringsAsFactors=F)
# get covariances of all subjects
#------------------------------------------
covs=read.xlsx("~/neurogen/predictHD/doc/20151104_scherzer_all_pax_unblinded_phenotypes.xlsx",1,stringsAsFactors=F)
covs=unique(covs[,-1])
rownames(covs)=covs[,1]
covs=covs[,-1]
covs$group=ifelse(covs$cap_d_group=="cont","control","HD")
covs$is_diagnosed[is.na(covs$is_diagnosed)]="Not_diagnosed"
#read in gene data from text file
#-----------------------------------------
expr=read.table("~/neurogen/predictHD/data/nanoString/miRNA.nanostring.raw.201606.txt",header = T,sep="\t",stringsAsFactors = F,fill = T,comment.char = "")
rownames(expr)=expr$Gene.Name
table(covs$Class.Name)


#confirm There must be have six positive control genes order by different concentrations 
#from high to low, since NanoString nCounter analyzer provide six different samll 
#positive controls with six different concen- trations
test=subset(expr,Class.Name=="Positive")


expr1=subset(expr,select=-c(1:3))


#----------------pick one from two replicates, and remove UR
subjectsID=assay$Subject.ID[match(colnames(expr1),assay$raw.data.filename)]
subjectsID1=unique(subjectsID[grep("UR",subjectsID,invert = T)])
colnames(expr1)=subjectsID
rrrrrr=expr1[,subjectsID1]
EXPR=data.frame(Code_Class=expr$Class.Name,miRNA=expr$Gene.Name,access=expr$Accession..,rrrrrr)
Type=covs$group[match(subjectsID1,rownames(covs))]
colnames(EXPR)[4:220]=Type


trainSet=read.xlsx("~/neurogen/predictHD/doc/20150529_scherzer_randomization_request_with_demography.xls",1,stringsAsFactors=F)
sample2subject=read.xlsx("~/neurogen/predictHD/doc/PREDICT-HD.Samples.Summary.20150520.xlsx",1,stringsAsFactors=F)
trainSetSubjectID1=sample2subject$SubjectID[match(trainSet$SampleID,sample2subject$SampleID)]

rrrrrr_train=rrrrrr[,trainSetSubjectID1]
EXPR_train=data.frame(Code_Class=expr$Class.Name,miRNA=expr$Gene.Name,access=expr$Accession..,rrrrrr_train)
Type=covs$group[match(trainSetSubjectID1,rownames(covs))]
colnames(EXPR_train)[4:103]=Type

rrrrrr_test=rrrrrr[,!subjectsID1%in%trainSetSubjectID1]
EXPR_test=data.frame(Code_Class=expr$Class.Name,miRNA=expr$Gene.Name,access=expr$Accession..,rrrrrr_test)
Type=covs$group[match(subjectsID1[!subjectsID1%in%trainSetSubjectID1],rownames(covs))]
colnames(EXPR_test)[4:120]=Type

designs=data.frame(group=colnames(EXPR_train)[4:103])
train_endogenous=as.matrix(EXPR_train[1:800,4:103])
train_housekeeping=as.matrix(EXPR_train[801:805,4:103])
train_negative=as.matrix(EXPR_train[806:811,4:103])
train_positive=as.matrix(EXPR_train[812:817,4:103])

NanoStringData1=createNanoStringSet(train_endogenous,train_positive,train_negative,train_housekeeping,designs)
pheno=pData(NanoStringData1)
group=pheno$group
design.full=model.matrix(~0+group)
contrast=c(-1,1)
NanoStringData1=estNormalizationFactors(NanoStringData1)
result_train=glm.LRT(NanoStringData1,design.full,contrast=contrast)
train_DE=rownames(result_train$table)[result_train$table$pvalue<0.005]



designs=data.frame(group=colnames(EXPR_test)[4:120])
test_endogenous=as.matrix(EXPR_test[1:800,4:120])
test_housekeeping=as.matrix(EXPR_test[801:805,4:120])
test_negative=as.matrix(EXPR_test[806:811,4:120])
test_positive=as.matrix(EXPR_test[812:817,4:120])

NanoStringData1=createNanoStringSet(test_endogenous,test_positive,test_negative,test_housekeeping,designs)
pheno=pData(NanoStringData1)
group=pheno$group
design.full=model.matrix(~0+group)
contrast=c(-1,1)
NanoStringData1=estNormalizationFactors(NanoStringData1)
result_test=glm.LRT(NanoStringData1,design.full,contrast=contrast)
test_DE=rownames(result_test$table)[result_test$table$pvalue<0.005]

#train_DE=rownames(result_train$table)[result_train$table$qvalue<0.05]
#test_DE=rownames(result_test$table)[result_test$table$qvalue<0.05]
#intersect(train_DE,test_DE)
