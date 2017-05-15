##HD biomarker analysis
##default path /Users/zhuowang
#negative contral >10% positive scaling [0.3,3] global scaling [0.1,10]
#replicates 16 pairs without UR
#before 800 240 samples after QC and norm 727 miRNAs and 194 samples

#log10(expr+0.1)
#start feature=10 final feature=8 `hsa-miR-146b-5p` `hsa-miR-521`  `hsa-miR-581`  `hsa-miR-562`  `hsa-miR-323a-3p`  `hsa-miR-191-5p` `hsa-miR-139-3p``hsa-miR-329`
#cross validation DE from 1 to 20 100 times
source("http://bioconductor.org/biocLite.R")
rm(list=ls())
library(xlsx)
library(sva)
library(MASS)
library(ggplot2)
library(ROCR)
library(pROC)
library(VennDiagram)
#load data read in patient data from excel file
#----------------------------------------------
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
#------------------------------------
expr=read.table("~/neurogen/predictHD/data/nanoString/miRNA.nanostring.raw.201606.txt",header = T,sep="\t",stringsAsFactors = F,fill = T,comment.char = "")
rownames(expr)=expr$Gene.Name
table(expr$Class.Name)

##----------------------------------------------
##quanlity control
#using negative control miRNAs to filter outlier miRNAs
negative=subset(expr,Class.Name=="Negative",select = -c(1:3))
scale=colSums(negative)/nrow(negative)
filter=apply(expr[,-(1:3)],1,function(x) length(x[x>scale])>floor(ncol(expr)/10))
expr=expr[filter,]
table(expr$Class.Name)
lll=t(t(expr[,-(1:3)])-scale)
expr[,-(1:3)]=lll
for (i in 1:nrow(expr))
{for (j in 4:(ncol(expr)))
{
  if(expr[i,j]<0)
    expr[i,j]=0
}}


#------------------------------------------------------------
#using positive control miRNAs and Endogenous miRNA to filter samples
#1.normalize to positive control
positive=subset(expr,Class.Name=="Positive",select=-c(1:3))
control=positive[,which(covs[assay$Subject.ID[match(colnames(positive),assay$raw.data.filename)],'group']=='control')]
disease=positive[,which(covs[assay$Subject.ID[match(colnames(positive),assay$raw.data.filename)],'group']=='HD')]

positive=colSums(positive)
positive_scaling_factor=mean(positive)/positive
positive_scaling_factor2diag=diag(positive_scaling_factor)
colnames(positive_scaling_factor2diag)=names(positive_scaling_factor)
rownames(positive_scaling_factor2diag)=names(positive_scaling_factor)
expr[,-c(1:3)]=as.matrix(expr[,-c(1:3)])%*% positive_scaling_factor2diag
#remove samples with positive_scaling_factor outside a range [0.3,3]
toberemoved1=colnames(expr[,-c(1:3)])[positive_scaling_factor<0.3 | positive_scaling_factor>3]

#2.globally normalize to Endogenous miRNAs
expr=subset(expr,Class.Name=="Endogenous1",select=-c(1:3))
total=colSums(expr)
globaltotal_scaling_factor=mean(total)/total
globaltotal_scaling_factor2diag=diag(globaltotal_scaling_factor)
colnames(globaltotal_scaling_factor2diag)=names(globaltotal_scaling_factor)
rownames(globaltotal_scaling_factor2diag)=names(globaltotal_scaling_factor)
expr=as.matrix(expr)%*% globaltotal_scaling_factor2diag
#remove samples with reference_scaling_factor outside a range [0.1,10]
toberemoved2=colnames(expr)[globaltotal_scaling_factor<0.1 | globaltotal_scaling_factor>10]

#-----------------remove replicates
subjectsID=assay$Subject.ID[match(colnames(expr),assay$raw.data.filename)]
rep=names(table(subjectsID[!grepl("UR",subjectsID)]))[table(subjectsID[!grepl("UR",subjectsID)])>1]

#-------------------------filter out samples(N=240-->194)

expr0=expr
N1=ncol(expr)

#remove samples with normalization flag
toberemoved=unique(c(toberemoved1,toberemoved2))
expr=expr[,!(colnames(expr)%in%toberemoved)]
expr=expr[,grep("HD_",colnames(expr),invert = T)]

#----------------pick one from two replicates, and remove UR
subjectsID=assay$Subject.ID[match(colnames(expr),assay$raw.data.filename)]
expr=expr[,match(unique(subjectsID[grep("UR",subjectsID,invert = T)]),subjectsID)]
N2=ncol(expr)
message(paste("# Before filtering: N=",N1,"\n# After filtering: N=",N2))
#-------------construct EXP, COV for modeling
#-------------get the batch, plate, and RIN info
batch=sub("batch","",assay$batch[match(colnames(expr),assay$raw.data.filename)])
plate=sub("-.*","",assay$NanoString.ID[match(colnames(expr),assay$raw.data.filename)])
#use SubjectID from sample table as colnames
colnames(expr)=assay$Subject.ID[match(colnames(expr),assay$raw.data.filename)]

EXPR=as.matrix(log10(expr+0.1))
COVS=cbind(covs[colnames(expr),],batch=batch,plate=plate)
col2facotr=c('visit','site','gender','race','ethnicity','dcl','is_diagnosed','cap_d_group','baseline_cap_group','batch','plate','group')
COVS[,col2facotr]=lapply(COVS[,col2facotr],factor)
col2num=c("RIN","visit_age","cap_d_score","baseline_cap","motor_total","education_max")
COVS[,col2num]=lapply(COVS[,col2num],as.numeric)

data=cbind(t(EXPR),COVS)
table(data[,which(colnames(data)=="group")])
#-------------dividing into training and text sets
#--------training set
trainSet=read.xlsx("~/neurogen/predictHD/doc/20150529_scherzer_randomization_request_with_demography.xls",1,stringsAsFactors=F)
sample2subject=read.xlsx("~/neurogen/predictHD/doc/PREDICT-HD.Samples.Summary.20150520.xlsx",1,stringsAsFactors=F)
trainSetSubjectID1=sample2subject$SubjectID[match(trainSet$SampleID,sample2subject$SampleID)]
trainSetSubjectID=intersect(intersect(trainSetSubjectID1,rownames(covs)),colnames(expr))
trainData=EXPR[,trainSetSubjectID] ##742 rows 95 columns
trainPheno=COVS[trainSetSubjectID,] ##95 rows 19 columns

testData = EXPR[,!(colnames(expr) %in% trainSetSubjectID)];
testPheno = COVS[!(colnames(expr) %in% trainSetSubjectID),] #99

### ------------------------------------------------
### Use sva to call differentially expressed (DE) genes and logistic regression
### ------------------------------------------------
# sva adjusted
trainMod = model.matrix(~group+gender+visit_age,data=trainPheno)
trainMod0 = model.matrix(~gender+visit_age,data=trainPheno)
pValuesSv=f.pvalue(trainData,trainMod,trainMod0)
svobj=sva(trainData,trainMod,trainMod0)
fsvaobj=fsva(trainData,trainMod,svobj,newdat = trainData)
new=fsvaobj$db
new_test=fsva(trainData,trainMod,svobj,newdat = testData)$new

#DEgenes_both_ttest_and_sva<-names(sort(pValuesSv[pValuesSv<=0.05]))

DEgenes_both_ttest_and_sva<-names(pValuesSv[order(pValuesSv)[1:10]])

#trainData.ExpCov.combined=cbind(t(trainData),trainPheno)
trainData.ExpCov.combined=cbind(t(new),trainPheno)
HDHC.miRNALogistic<-glm(paste0("group~",paste0("`",DEgenes_both_ttest_and_sva,"`",collapse="+")),
                        family=binomial(logit),data=trainData.ExpCov.combined)
HDHC.miRNALogistic.AIC=stepAIC(HDHC.miRNALogistic,direction="backward")

p<-predict(HDHC.miRNALogistic.AIC,newdata=trainData.ExpCov.combined,type=c("response"))
pp=roc(trainData.ExpCov.combined$group,p,plot =TRUE, main="ROC on training set", print.thres = FALSE, print.auc =TRUE,col="blue",auc=TRUE)
coords(pp, "best")
combined=cbind(t(new_test),testPheno)

q<-predict(HDHC.miRNALogistic.AIC,newdata=combined,type=c("response"))
qq=roc(combined$group,q,plot =TRUE, main="ROC on test set", print.thres = FALSE, print.auc =TRUE,col="1",auc=TRUE)
coords(qq, "best")


final_feature=c("hsa-miR-146b-5p","hsa-miR-329","hsa-miR-139-3p","hsa-miR-581","hsa-miR-562",
                "hsa-miR-323a-3p","hsa-miR-191-5p","hsa-miR-521")
ppp=round(pValuesSv[final_feature],3)
Y=t(trainData) 
hh=Y[,final_feature]
hhc=hh[which(COVS[rownames(hh),'group']=='control'),]
hhd=hh[which(COVS[rownames(hh),'group']=='HD'),]
o=rbind(data.frame(hhc,type="control"),data.frame(hhd,type="HD"))
k1=ggplot(o,aes_string("type",colnames(o)[1],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[1],"/p=",ppp[1]))+ylab("log10(expression+0.1)")
k2=ggplot(o,aes_string("type",colnames(o)[2],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[2],"/p=",ppp[2]))+ylab("log10(expression+0.1)")
k3=ggplot(o,aes_string("type",colnames(o)[3],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[3],"/p=",ppp[3]))+ylab("log10(expression+0.1)")
k4=ggplot(o,aes_string("type",colnames(o)[4],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[4],"/p=",ppp[4]))+ylab("log10(expression+0.1)")
k5=ggplot(o,aes_string("type",colnames(o)[5],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[5],"/p=",ppp[5]))+ylab("log10(expression+0.1)")
k6=ggplot(o,aes_string("type",colnames(o)[6],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[6],"/p=",ppp[6]))+ylab("log10(expression+0.1)")
k7=ggplot(o,aes_string("type",colnames(o)[7],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[7],"/p=",ppp[7]))+ylab("log10(expression+0.1)")
k8=ggplot(o,aes_string("type",colnames(o)[8],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(final_feature[8],"/p=",ppp[8]))+ylab("log10(expression+0.1)")
################################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 4)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(k1, vp = vplayout(1, 1))
print(k2, vp = vplayout(1, 2))
print(k3, vp = vplayout(1, 3))
print(k4, vp = vplayout(1, 4))
print(k5, vp = vplayout(2, 1))
print(k6, vp = vplayout(2, 2))
print(k7, vp = vplayout(2, 3))
print(k8, vp = vplayout(2, 4))

#cross validation in training set 75% train 25% validation

trainData_raw=new
trainPheno_raw=trainPheno
testData_raw=new_test
testPheno_raw=testPheno


value=vector()
type=vector()
DE_number=vector()
auc_optimal<-vector()
auc_train<-vector()
auc_test<-vector()

for(j in 1:20)
{
  for(i in 1:100)
  {
    id=sample(1:95,70)
    trainData=trainData_raw[,id]
    trainPheno=trainPheno_raw[id,]
    testData=trainData_raw[,!c(1:95)%in%id]
    testPheno=trainPheno_raw[!c(1:95)%in%id,]
    
    combined=cbind(t(testData),testPheno)
    trainData.ExpCov.combined=cbind(t(trainData),trainPheno)
    
    trainMod = model.matrix(~group,data=trainPheno)
    trainMod0 = model.matrix(~1,data=trainPheno)
    pValuesSv=f.pvalue(trainData,trainMod,trainMod0)
    DEgenes_both_ttest_and_sva<-names(pValuesSv[order(pValuesSv)[1:j]])
    
    HDHC.miRNALogistic<-glm(paste0("group~",paste0("`",DEgenes_both_ttest_and_sva,"`",collapse="+")),
                            family=binomial(logit),data=trainData.ExpCov.combined)
    HDHC.miRNALogistic.AIC=stepAIC(HDHC.miRNALogistic,direction="backward") 
    
    p<-predict(HDHC.miRNALogistic.AIC,newdata=trainData.ExpCov.combined,type=c("response"))
    pp=roc(trainData.ExpCov.combined$group,p,plot =TRUE, main="ROC for Model Assessment", print.thres = FALSE, print.auc =TRUE,col="blue",auc=TRUE)
    
    
    q<-predict(HDHC.miRNALogistic.AIC,newdata=combined,type=c("response"))
    qq=roc(combined$group,q,plot =TRUE, main="ROC", print.thres = FALSE, print.auc =TRUE,col="1",auc=TRUE)
    
    real_test_combine=cbind(t(testData_raw),testPheno_raw)
    real_test=predict(HDHC.miRNALogistic.AIC,newdata = real_test_combine,type=c("response"))
    
    
    auc_train[i]=auc(trainData.ExpCov.combined$group,p)
    auc_optimal[i]=auc(combined$group,q)
    auc_test[i]=auc(real_test_combine$group,real_test)
  }
  value=c(value,auc_train,auc_optimal,auc_test)
  type=c(type,rep("auc_train",100),rep("auc_optimal",100),rep("auc_test",100))
  DE_number=c(DE_number,rep(paste(j),300))
}

result=data.frame(DE_number,value,type)
DE_Number=factor(result$DE_number,levels=1:20)
result=data.frame(DE_Number,value,type)
ggplot(result,aes(DE_Number,value,col=type))+geom_boxplot()


