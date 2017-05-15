#whole pipeline: see predictHD.pdf (gliff.com)
#step 1: Quality control and normalization
#step 2: Divide into training and test set
#step 3: Effect removing
####repeat 100 times step 4 to step 7
#step 4: In adjusted training set: Raw feature selection
#step 5: Model building
#step 6: Model selection
#step 7: Model assessment
#step 8: Features' frequency
#step 9: Final model 
## Predictive diagostics in adjusted test set

library(stringr)
library(xlsx)
library(sva)
library(pROC)
library(ROCR)
rm(list=ls())
##load data
#------------------------------------------
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
table(expr$Class.Name)

#-------------------------------------------
#step 1: Quality control and normalization
#------------------------------------------=
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

#-----------------------------------------------
#step 2: dividing into training and text sets
#-----------------------------------------------

trainSet=read.xlsx("~/neurogen/predictHD/doc/20150529_scherzer_randomization_request_with_demography.xls",1,stringsAsFactors=F)
sample2subject=read.xlsx("~/neurogen/predictHD/doc/PREDICT-HD.Samples.Summary.20150520.xlsx",1,stringsAsFactors=F)
trainSetSubjectID1=sample2subject$SubjectID[match(trainSet$SampleID,sample2subject$SampleID)]
trainSetSubjectID=intersect(intersect(trainSetSubjectID1,rownames(covs)),colnames(expr))
trainData=EXPR[,trainSetSubjectID] ##727 rows 95 columns
trainPheno=COVS[trainSetSubjectID,] ##95 rows 19 columns

testData = EXPR[,!(colnames(expr) %in% trainSetSubjectID)];
testPheno = COVS[!(colnames(expr) %in% trainSetSubjectID),] #99

#---------------------------------------------
#step 3: Effect removing (sva adjusted)
#---------------------------------------------

trainMod = model.matrix(~group+gender+visit_age,data=trainPheno)
trainMod0 = model.matrix(~gender+visit_age,data=trainPheno)
svobj=sva(trainData,trainMod,trainMod0)
fsvaobj=fsva(trainData,trainMod,svobj,newdat = trainData)
new=fsvaobj$db
new_test=fsva(trainData,trainMod,svobj,newdat = testData)$new

#-------------------------------------------
#step 4-----step 7: repeat 100 times
#-------------------------------------------

trainData_raw=new
trainPheno_raw=trainPheno
testData_raw=new_test
testPheno_raw=testPheno
raw_train=cbind(t(trainData_raw),trainPheno_raw)
raw_test=cbind(t(testData_raw),testPheno_raw)

auc_test=vector()
auc_train=vector()
j=vector()
repeat_time=100
for(k in 1:repeat_time)
{
  set.seed(k)
  #-----------random split into 75% and 25%------------------
  id=sample(1:95,70)
  trainData=trainData_raw[,id]
  trainPheno=trainPheno_raw[id,]
  testData=trainData_raw[,!c(1:95)%in%id]
  testPheno=trainPheno_raw[!c(1:95)%in%id,]
  combined=cbind(t(testData),testPheno)
  trainData.ExpCov.combined=cbind(t(trainData),trainPheno)
  #------------step 4: raw feature selection----------------
  trainMod = model.matrix(~group,data=trainPheno)
  trainMod0 = model.matrix(~1,data=trainPheno)
  pValuesSv=f.pvalue(trainData,trainMod,trainMod0)
  DEgenes_both_ttest_and_sva<-names(pValuesSv[order(pValuesSv)[1:8]])
  #------------step 5&6: Model building and AIC selection-----------------------
  HDHC.miRNALogistic<-glm(paste0("group~",paste0("`",DEgenes_both_ttest_and_sva,"`",collapse="+")),
                          family=binomial(logit),data=trainData.ExpCov.combined)
  HDHC.miRNALogistic.AIC=stepAIC(HDHC.miRNALogistic,direction="backward") 
  #------------step 7: Model assessment------------------
  p<-predict(HDHC.miRNALogistic.AIC,newdata=trainData.ExpCov.combined,type=c("response"))
  q<-predict(HDHC.miRNALogistic,newdata=combined,type=c("response"))
  auc_train[k]=roc(trainData.ExpCov.combined$group,p)$auc
  auc_test[k]=roc(combined$group,q)$auc
  #-----------step 8: Features' frequency---------------
  if (auc_test[k]>0.5)
  {
    jj=strsplit(as.character(HDHC.miRNALogistic.AIC$formula)[3],"[``]")[[1]]
    jj=jj[str_detect(jj,"a")]
    j=cbind(t(jj),j)
  }
}
#-----------step 9: Final features
result_fresquence=as.matrix(table(j))
DE=names(result_fresquence[which(result_fresquence>=floor(repeat_time/5)),])

if(!is.null(DE))
{
  par(mfrow=c(1,2),las="1")
  HDHC.miRNALogistic<-glm(paste0("group~",paste0("`",DE,"`",collapse="+")),
                          family=binomial(logit),data=raw_train)
  
  p<-predict(HDHC.miRNALogistic,newdata=raw_train,type=c("response"))
  pp=roc(raw_train$group,p,plot =TRUE, main="ROC on training set", print.thres = FALSE, print.auc =TRUE,col="blue",auc=TRUE)
  coords(pp, "best")
  q<-predict(HDHC.miRNALogistic,newdata=raw_test,type=c("response"))
  qq=roc(raw_test$group,q,plot =TRUE, main="ROC on test set", print.thres = FALSE, print.auc =TRUE,col="1",auc=TRUE)
  coords(qq, "best")
}

#threshold specificity sensitivity  Feq>repeat_time/3
#0.2306421   0.3396226   0.8043478 
#threshold specificity sensitivity  Feq>repeat_time/4
#0.3181322   0.6415094   0.4565217  
#threshold specificity sensitivity  Feq>repeat_time/5
#0.9818422   0.9056604   0.1956522 

##-------------Frequency plot------------------------
par(mfrow=c(1,1),las="2")
Frequency_95_samples=result_fresquence[order(result_fresquence,decreasing = TRUE)[1:15],]
kkkk=barplot(Frequency_95_samples,cex.names=0.6,main="top 15 (100 times repetition)",ylab="Frequency",ylim=c(0,80))
text(kkkk, Frequency_95_samples, format(Frequency_95_samples), col = "blue")
#####--------------export the ROC curve on training set and test set------------
#####--------------export the frequency of top 10 miRNAs------------
