#blood sample nanostring data
###########################################################
rm(list=ls())
library(xlsx)
library(sva)
library(VennDiagram)
library(grid)
library(ggplot2)
library(fitdistrplus)
library(pheatmap)
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
#table(covs$cap_d_group)
#cont high  med 
#111   100  15 

#read in gene data from text file
#------------------------------------
expr=read.table("~/neurogen/predictHD/data/nanoString/miRNA.nanostring.raw.201606.txt",header = T,sep="\t",stringsAsFactors = F,fill = T,comment.char = "")
rownames(expr)=expr$Gene.Name
table(expr$Class.Name)
#Endogenous1 Housekeeping     Negative     Positive      SpikeIn 
#800            5            6            6            5 

##----------------------------------------------
##quanlity control
#using negative control miRNAs to filter outlier miRNAs
negative=subset(expr,Class.Name=="Negative",select = -c(1:3))
scale=colSums(negative)/nrow(negative)
filter=apply(expr[,-(1:3)],1,function(x) length(x[x>scale])>floor(ncol(expr)/10))
expr=expr[filter,]
table(expr$Class.Name)

#Endogenous1 Housekeeping     Negative     Positive      SpikeIn 
#727            5            4            6            5 
lll=t(t(expr[,-(1:3)])-scale)
expr[,-(1:3)]=lll
for (i in 1:nrow(expr))
{for (j in 4:(ncol(expr)))
{
  if(expr[i,j]<0)
    expr[i,j]=0
}}
expr[,-c(1:3)]=expr[,-c(1:3)]+1

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

#######Figure
pdf("replicates.pdf", paper='us') 
par(mfrow=c(3,3))
for(i in rep)
{
  e=data.matrix(expr[, grep(i,subjectsID)])+0.1
  par(pty="s"); #generates square plotting region
  #plot(e, log='xy', main=i, xlab=names(e)[1], ylab=names(e)[2], ylim=range(e), xlim=range(e), pch=1, cex=.8, asp=1, xaxt="n", yaxt="n")
  plot(e, log='xy', main=i, xlab="Normalized miRNA count", ylab="Normalized miRNA count", ylim=range(e), xlim=range(e), pch=1, cex=.8, asp=1, xaxt="n", yaxt="n")
  
  abline(a=0,b=1, lty=2, col='red') #add straight line to plot
  pcc=cor(as.numeric(e[,1]), as.numeric(e[,2])) #compute correlation
  legend("topleft", paste0("Pearson's r = ", round(pcc,7)), bty='n')
  #print legend, pcc vales as Pearons's r, rounded to 7 decimal places
  #bty='n' suppresses box around plot
  
  x <- floor(log10(range(e)))
  pow <- seq(x[1], x[2]+1)
  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  #note global assignment within function using <-
  axis(1, 10^pow)
  axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
  axis(2, 10^pow)
  axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)    
}

# for UR all-to-all only
e=log10(data.matrix(expr[, grep("UR",subjectsID)]))
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(e, lower.panel = panel.smooth, upper.panel = panel.cor,asp=1, main="Pearson's correlations for UR samples (log10)")
dev.off()

#-------------------------filter out samples(N=240-->194)
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
EXPR=as.matrix(log10(expr))
COVS=cbind(covs[colnames(expr),],batch=batch,plate=plate)
col2facotr=c('visit','site','gender','race','ethnicity','dcl','is_diagnosed','cap_d_group','baseline_cap_group','batch','plate','group')
COVS[,col2facotr]=lapply(COVS[,col2facotr],factor)
col2num=c("RIN","visit_age","cap_d_score","baseline_cap","motor_total","education_max")
COVS[,col2num]=lapply(COVS[,col2num],as.numeric)
data=cbind(t(EXPR),COVS)
data_raw=cbind(t(as.matrix(expr)),COVS)
table(data[,which(colnames(data)=="group")])

############################################################
#whole dataset to perform DE analysis
############################################################
#data  194*746
#EXPR  727*194
#COVS  194*19
wholeMod = model.matrix(~group+gender+visit_age,data=COVS)
wholeMod0 = model.matrix(~gender+visit_age,data=COVS)
whole_pValuesSv=f.pvalue(EXPR,wholeMod,wholeMod0)
whole_qValues_BH = p.adjust(whole_pValuesSv,method="BH")
whole_qValues_fdr = p.adjust(whole_pValuesSv,method="fdr")
whole_DE_sva<-names(whole_pValuesSv[order(whole_pValuesSv)[1:length(which(whole_pValuesSv<0.05))]])
#whole_DE_sva=c("hsa-miR-135a-5p","hsa-miR-24-3p","hsa-miR-3185","hsa-miR-577")
pp=round(whole_pValuesSv[whole_DE_sva],3)

whole_data=cbind(t(EXPR),COVS)
whole_HD<-subset(whole_data,group %in%  "HD")[,1:727]
whole_control<-subset(whole_data,group %in% "control")[,1:727]
log2fc=log2(apply(whole_HD,2,mean)/apply(whole_control,2,mean))
res=data.frame(log2FoldChange=log2fc,pvalue=whole_pValuesSv,fdr=whole_qValues_fdr,BH=whole_qValues_BH)
write.csv(res,"res.csv")
####Volcano plot#############
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot(p<0.05 |lf|(HD/control)>0.5)", xlab="log2FoldChange(HD/control)",xlim=c(-3,3)))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
library(calibrate)
with(subset(res, pvalue<.05 & abs(log2FoldChange)>0.5), textxy(log2FoldChange, -log10(pvalue), labs=rownames(subset(res, pvalue<.05 & abs(log2FoldChange)>0.5)), cex=.6))
lines(y=c(0,2),x=c(1,1),col='red')
lines(y=c(1.3013,1.3013),x=c(-3,3),col='red')


pdf("miRNA.nanostring.volcano.pdf", width=5, height=5, useDingbats=F)
plot(-log10(pvalue) ~ log2FoldChange, res, pch=20, col=ifelse(pvalue<0.05,'red','black'), xlab="log2(fold change of HD/HC)", main="Volcano plot for DE miRNAs in blood")
with(subset(res, abs(log2FoldChange)>0.5&pvalue<0.05) %>% mutate(names=gsub("hsa-","",rownames(.))), text(log2FoldChange, -log10(pvalue), names, cex=.5, pos=4))
abline(v=c(0.5,-0.5),h=-log10(0.05), lty=2, col='gray')
dev.off()
#################################################
###heatmap#############


#test=subset(res,pvalue<0.05)
test=subset(res,pvalue<0.05& abs(log2FoldChange)>0.5)
test=test[order(test$log2FoldChange,decreasing=F),]
ddd=data_raw[order(data_raw[,744],decreasing=T),]
ttest=t(log10(ddd[,rownames(test)]+0.1))
ttest=t(scale(t(ttest),scale=T,center=F))

mat_data_sorted = t(apply(ttest, 1, function(x) c(sort(x[1:102],decreasing = T),sort(x[103:194],decreasing = T))))

annotation_col =data.frame(CellType=ddd$group)
rownames(annotation_col) <- colnames(ttest)
colnames(mat_data_sorted) = colnames(ttest)
pheatmap(mat_data_sorted,annotation_col = annotation_col,cluster_rows = T,cluster_cols = F,show_colnames=F,main ="Blood Samples/ fold change >2")
#pheatmap(mat_data_sorted,annotation_col = annotation_col,cluster_rows = T,cluster_cols = F,show_colnames=F,main ="Blood Samples/ P value<0.05")

#### ------------------------------------------------
#### barplot for 20 DE Use linear regression to call DE genes 
#### ------------------------------------------------

result_group=data.frame();
Y=t(EXPR) 

pdf("miRNA.nanostring.barplot.topDE.pdf", width=2, height=3)
for(i in 1:ncol(Y)){
  gname = colnames(Y)[i];
  df=cbind(y=Y[,i],COVS)
  
  #mod = lm(as.formula(paste("y ~", paste(colnames(COVS),collapse = " + "))), df) # complete mode
  mod = lm(as.formula("y ~ group + gender + visit_age"), df) # simplified
  
  df$y=10^df$y;
  message(paste(i,gname))
  
  ## use a barplot (incl. error bar) starting from y=0
  df.summary <- data.frame(
    group=levels(df$group),
    mean=tapply(df$y, df$group, mean),
    n=tapply(df$y, df$group, length),
    sd=tapply(df$y, df$group, sd)
  )
  df.summary$sem <- df.summary$sd/sqrt(df.summary$n)
  
  r=summary.lm(mod)$coefficients['groupHD',]
  r=c(r, as.vector(t(df.summary[,2:5])), as.numeric(log2(df.summary$mean[2]/df.summary$mean[1])))
  names(r)[5:length(r)]=c(as.vector(sapply(rownames(df.summary), function(x) paste(x, colnames(df.summary)[2:5], sep="."))),'log2_foldchange');
  
  if(ncol(result_group)==0) {
    result_group = as.data.frame(r)
  } else {
    result_group = cbind(result_group, r)
  }
  #changed 3 instances of group1 to groupHD
  p=r[4]
  if(p<0.05 || gname==""){  # plot box plot for DE genes
    g=ggplot(df.summary, aes(x = group, y = mean)) +  
      geom_bar(width=.5, position = position_dodge(), stat="identity", fill=c('blue2','red2')) + 
      geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=0.25) + 
      ylab("Normalized expression") + 
      ggtitle(paste0(gname,'\n(p=',round(p,3),"; log2(HD/HC)=",round(as.numeric(log2(df.summary$mean[2]/df.summary$mean[1])),2))) + 
      theme_bw() + 
      #geom_text(data=NULL, x=Inf, y=Inf, label=paste0('p-value: ',round(p,4)), hjust=-1) + 
      theme(panel.grid.major = element_blank(), plot.title = element_text(size = rel(.6)))
    print(g)
  }
}
dev.off()

####################################################################################
##########two Pearson’s r distribution plot(0 with no 0) ############################
whole_data=cbind(t(EXPR),COVS)

whole_HD<-subset(whole_data,group %in%  "HD")
whole_control<-subset(whole_data,group %in% "control")
whole_HD_matrix=whole_HD[,(1:nrow(EXPR))]
whole_control_matrix=whole_control[,(1:nrow(EXPR))]
df_whole=t(rbind(whole_control_matrix,whole_HD_matrix)) #control 92 HD 102
# df_res=t(apply(df_whole, 1, function(x) {hc=x[1:92]; hd=x[93:194]; hc=hc[hc>-1]; hd=hd[hd>-1]; 
#               pvalues=t.test(hc, hd)$p.value; n_HC=length(hc); n_HD=length(hd); 
#               log2FC=log2(mean(hd)/mean(hc));c(n_HD, n_HC, log2FC, pvalues);}))
# BH=p.adjust(df_res[,4],method = "BH")
# fdr=p.adjust(df_res[,4],method="fdr")
# bonferroni=p.adjust(df_res[,4],method = "bonferroni")
# df_res_all=data.frame(df_res,BH,fdr,bonferroni)

HD_cap_pearson=cor(whole_HD_matrix,whole_HD$cap_d_score,method = "pearson")

###########################################################
df=as.numeric(HD_cap_pearson)
fitn=fitdist(df,'norm')
pdf("Blood_pearson_distribution.pdf", width=8.5, height=6)
hist(df, breaks=30, prob=TRUE, xlab='Pearson coefficient', main='Distribution of Pearson coefficient',xlim=c(-0.3,0.3))
lines(density(df, bw=0.05))
m=round(as.numeric(fitn$estimate[1]),digits=3)
sd=round(as.numeric(fitn$estimate[2]),digits=3)
lines(density(rnorm(n=2000000, mean=m, sd=sd),bw=0.05), col='blue',lty=2)
p=round(qnorm(.05, mean=m, sd=sd, lower.tail = F), digits=3)
q=round(qnorm(.95, mean=m, sd=sd, lower.tail = F), digits=3)
lines(y=c(0,0.3),x=c(p,p),col='red')
lines(y=c(0,0.3),x=c(q,q),col='red')
text(p+0.02,0.9,paste0("P(X>",p,") = 0.05"))
text(q,0.9,paste0("P(X<",q,") = 0.05"))
legend("topright", c("empirical density curve", paste0("fitted normal distribution \n(mean=",m,", sd=",sd,")")), col=c('black','blue'), lty=c(1,2), bty='n')
dev.off()
Blood_CAP=c(rownames(HD_cap_pearson)[HD_cap_pearson>p],rownames(HD_cap_pearson)[HD_cap_pearson<q])
venn.diagram(x=list(whole_feature=whole_DE_sva,CAP_feature=Blood_CAP),col=c(1,2),main = "Overlap miRNA",filename = "overlap",
             sub = "hsa-miR-135a-5p,hsa-miR-24-3p,hsa-miR-3185",width = 4000)

####################
colnames(whole_HD_matrix)[348]
miR_367_3p_expr=data.frame(expr=whole_HD$`hsa-miR-367-3p`,CAP_score=whole_HD$cap_d_score,HD_group=whole_HD$cap_d_group)
#tt=
k1=ggplot(miR_367_3p_expr,aes(CAP_score,expr,col=HD_group))+geom_point()+labs(title="Scatterplot for miR-367-3p (r = -0.3)")+ylab("log10(expr+0.1)")
k2=ggplot(miR_367_3p_expr,aes(x=HD_group,y=expr,col=HD_group))+geom_boxplot(position = position_dodge(1),width=0.5)+labs(title="Scatterplot for miR-367-3p (r=-0.3,p=0.06)")+ylab("log10(expr+0.1)")
miR_367_subset=subset(miR_367_3p_expr,expr>-1)
k3=ggplot(miR_367_subset,aes(CAP_score,expr,col=HD_group))+geom_point()+labs(title="Scatterplot for miR-367-3p (r = -0.3)")+ylab("log10(expr+0.1)")
k4=ggplot(miR_367_subset,aes(HD_group,expr,col=HD_group))+geom_boxplot(position = "identity")+labs(title="Scatterplot for miR-367-3p (r=-0.3,p=0.21)")+ylab("log10(expr+0.1)")

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(k1, vp = vplayout(1, 1))
print(k2, vp = vplayout(1, 2))
print(k3, vp = vplayout(2, 1))
print(k4, vp = vplayout(2, 2))

uu=subset(miR_367_subset,HD_group %in% "med")
vv=subset(miR_367_subset,HD_group %in% "high")
t.test(uu$expr,vv$expr)$p.value #0.2101473

uu=subset(miR_367_3p_expr,HD_group %in% "med")
vv=subset(miR_367_3p_expr,HD_group %in% "high")
t.test(uu$expr,vv$expr)$p.value #0.06

# ####################################################################################
# CAG=whole_HD$cap_d_score
# HD_cap_pearson_no0=apply(whole_HD_matrix,2,function(x) {cor(x[x>-1],CAG[x>-1])})
# df_no_zero=as.numeric(apply(whole_HD_matrix,2,function(x) {cor(x[x>-1],CAG[x>-1])}))
# df=df_no_zero
# fitn=fitdist(df,'norm')
# pdf("Blood_pearson_distribution_no0.pdf", width=8.5, height=6)
# hist(df, breaks=100, prob=TRUE, xlab='Pearson coefficient', main='Distribution of Pearson coefficient',xlim=c(-1,1))
# lines(density(df, bw=0.05))
# m=round(as.numeric(fitn$estimate[1]),digits=3)
# sd=round(as.numeric(fitn$estimate[2]),digits=3)
# lines(density(rnorm(n=2000000, mean=m, sd=sd),bw=0.05), col='blue',lty=2)
# p=round(qnorm(.05, mean=m, sd=sd, lower.tail = F), digits=3)
# q=round(qnorm(.95, mean=m, sd=sd, lower.tail = F), digits=3)
# lines(y=c(0,0.5),x=c(p,p),col='red')
# lines(y=c(0,0.5),x=c(q,q),col='red')
# text(p+0.4,0.7,paste0("P(X>",p,") = 0.05"))
# text(q-0.4,0.7,paste0("P(X<",q,") = 0.05"))
# legend("topright", c("empirical density curve", paste0("fitted normal distribution \n(mean=",m,", sd=",sd,")")), col=c('black','blue'), lty=c(1,2), bty='n')
# dev.off()
# 
# Blood_CAP_no0=c(names(HD_cap_pearson_no0)[HD_cap_pearson_no0>p],names(HD_cap_pearson_no0)[HD_cap_pearson_no0<q])
# venn.diagram(x=list(whole_feature=whole_DE_sva,CAP_feature=Blood_CAP_no0),col=c(1,2),main = "Overlap miRNA",filename = "overlap",sub = "hsa-miR-3185",width = 4000)
# 
# plot(whole_HD$`hsa-miR-577`, whole_HD$cap_d_score, main="Scatterplot for miR-577", xlab="log10(expr+0.1)", ylab="CAP score", pch=19)
# plot(whole_HD$`hsa-miR-3185`, whole_HD$cap_d_score, main="Scatterplot for miR-3185", xlab="log10(expr+0.1)", ylab="CAP score", pch=19)
# 
# plot(whole_HD$`hsa-miR-577`[which(whole_HD$`hsa-miR-577`>-1)], whole_HD$cap_d_score[which(whole_HD$`hsa-miR-577`>-1)], main="Scatterplot for miR-577 (r = 0.91)", xlab="log10(expr+0.1)", ylab="CAP score", pch=19)
# HD_cap_pearson_no0["hsa-miR-3185"]
# plot(whole_HD$`hsa-miR-3185`[which(whole_HD$`hsa-miR-3185`>-1)], whole_HD$cap_d_score[which(whole_HD$`hsa-miR-3185`>-1)], main="Scatterplot for miR-3185 (r = 0.28)", xlab="log10(expr+0.1)", ylab="CAP score", pch=19)


##########two Pearson’s r distribution plot(0 with no 0) ############################
####################################################################################




################figure 1 Volcano Plot##################################
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot(p<0.05 |lf|(HD/control)>1)", xlab="log2FoldChange(HD/control)",xlim=c(-3.1,3.1)))
# Add colored points: red if pvalue<0.05, orange of log2FC>0.5, green if both)
#with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
library(calibrate)
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=rownames(subset(res, pvalue<.05 & abs(log2FoldChange)>1)), cex=.5))
lines(y=c(0,2),x=c(1,1),col='red')
lines(y=c(0,2),x=c(-1,-1),col='red')
lines(y=c(1.3013,1.3013),x=c(-3,3),col='red')
#------------------------------------------------
#log10(0.05)=-1.3013
#############################figure 1 Volcano Plot#####################


#whole_DE_sva=c("hsa-miR-1252","hsa-miR-433","hsa-miR-553","hsa-miR-138-5p","hsa-miR-2053")
pp=round(whole_pValuesSv[order(whole_pValuesSv)[1:length(which(whole_pValuesSv<0.05))]],3)
pp=round(whole_pValuesSv[whole_DE_sva],3)
#Y=t(EXPR) 
Y=t(expr)
hh=Y[,whole_DE_sva]
hhc=hh[which(COVS[rownames(hh),'group']=='control'),]
hhd=hh[which(COVS[rownames(hh),'group']=='HD'),]
o=rbind(data.frame(hhc,type="control"),data.frame(hhd,type="HD"))

o=rbind(data.frame(miRNA=colnames(hhc),mean=colMeans(hhc),se=apply(hhc,2,sd),significant=pp,type="control",ymin=(mean-se)/sqrt(nrow(hhc))),
        data.frame(miRNA=colnames(hhd),mean=colMeans(hhd),se=apply(hhd,2,sd),significant=pp,type="HD"))
o=rbind(data.frame(miRNA=colnames(hhc),mean=colMeans(hhc),significant=pp,type="control",ymin=colMeans(hhc)-apply(hhc,2,sd)/sqrt(nrow(hhc)),ymax=colMeans(hhc)+apply(hhc,2,sd)/sqrt(nrow(hhc))),
        data.frame(miRNA=colnames(hhd),mean=colMeans(hhd),significant=pp,type="HD",ymin=colMeans(hhd)-apply(hhd,2,sd)/sqrt(nrow(hhd)),ymax=colMeans(hhd)+apply(hhd,2,sd)/sqrt(nrow(hhd))))

o$miRNA<- factor(o$miRNA, levels=c("hsa-miR-1252","hsa-miR-433","hsa-miR-553","hsa-miR-138-5p","hsa-miR-2053"), ordered=TRUE)
g1=ggplot(o, aes(x=miRNA, y=mean, fill=type)) +geom_bar(position=position_dodge(.8), stat="identity", width=.8)+
      geom_errorbar(aes(ymin=ymin, ymax=ymax),width=.5, # Width of the error bars
                size=.5,position=position_dodge(.8))+geom_text(aes(label=significant, y=ifelse(mean>1,ymax,ymin)), colour="black", position=position_dodge(.8), vjust=0, size=3)+ylab("Normalized expression")+labs(title="Barplot for DE miRNAs in blood")


k1=ggplot(o,aes_string("type",colnames(o)[1],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[1],"/p=",pp[1]))+ylab("Normalized expression")
k2=ggplot(o,aes_string("type",colnames(o)[2],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[2],"/p=",pp[2]))+ylab("Normalized expression")
k3=ggplot(o,aes_string("type",colnames(o)[3],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[3],"/p=",pp[3]))+ylab("Normalized expression")
k4=ggplot(o,aes_string("type",colnames(o)[4],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[4],"/p=",pp[4]))+ylab("Normalized expression")
k5=ggplot(o,aes_string("type",colnames(o)[5],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[5],"/p=",pp[5]))+ylab("Normalized expression")
k6=ggplot(o,aes_string("type",colnames(o)[6],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[6],"/p=",pp[6]))+ylab("Normalized expression")
k7=ggplot(o,aes_string("type",colnames(o)[7],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[7],"/p=",pp[7]))+ylab("Normalized expression")
k8=ggplot(o,aes_string("type",colnames(o)[8],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[8],"/p=",pp[8]))+ylab("Normalized expression")
k9=ggplot(o,aes_string("type",colnames(o)[9],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[9],"/p=",pp[9]))+ylab("Normalized expression")
k10=ggplot(o,aes_string("type",colnames(o)[10],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[10],"/p=",pp[10]))+ylab("Normalized expression")
#################################################

k1=ggplot(o,aes_string("type",colnames(o)[11],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[11],"/p=",pp[11]))+ylab("Normalized expression")
k2=ggplot(o,aes_string("type",colnames(o)[12],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[12],"/p=",pp[12]))+ylab("Normalized expression")
k3=ggplot(o,aes_string("type",colnames(o)[13],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[13],"/p=",pp[13]))+ylab("Normalized expression")
k4=ggplot(o,aes_string("type",colnames(o)[14],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[14],"/p=",pp[14]))+ylab("Normalized expression")
k5=ggplot(o,aes_string("type",colnames(o)[15],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[15],"/p=",pp[15]))+ylab("Normalized expression")
k6=ggplot(o,aes_string("type",colnames(o)[16],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[16],"/p=",pp[16]))+ylab("Normalized expression")
k7=ggplot(o,aes_string("type",colnames(o)[17],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[17],"/p=",pp[17]))+ylab("Normalized expression")
k8=ggplot(o,aes_string("type",colnames(o)[18],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[18],"/p=",pp[18]))+ylab("Normalized expression")
k9=ggplot(o,aes_string("type",colnames(o)[19],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[19],"/p=",pp[19]))+ylab("Normalized expression")
k10=ggplot(o,aes_string("type",colnames(o)[20],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0(whole_DE_sva[20],"/p=",pp[20]))+ylab("Normalized expression")

####################################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 5)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(k1, vp = vplayout(1, 1))
print(k2, vp = vplayout(1, 2))
print(k3, vp = vplayout(1, 3))
print(k4, vp = vplayout(1, 4))
print(k5, vp = vplayout(1, 5))
print(k6, vp = vplayout(2, 1))
print(k7, vp = vplayout(2, 2))
print(k8, vp = vplayout(2, 3))
print(k9, vp = vplayout(2, 4))
print(k10, vp = vplayout(2, 5))
#######################################################
#miR-577 and miR-3185 correlated with age gender and batch
whole_577_3185=data.frame("miR-577"=factor(whole_data$`hsa-miR-577`==-1,labels = c("expressed","un-expressed")),"miR-3185"=factor(whole_data$`hsa-miR-3185`==-1,labels = c("expressed","un-expressed")),
                          age=whole_data$visit_age,gender=as.factor(whole_data$gender),batch=as.factor(whole_data$batch))
ggplot(whole_577_3185,aes(miR.577,age,fill=miR.577))+geom_boxplot()
ggplot(whole_577_3185,aes(miR.3185,age,fill=miR.3185))+geom_boxplot()

table(whole_577_3185$miR.577)
table(whole_577_3185$miR.3185)
dim(subset(whole_577_3185, miR.577=="expressed" & gender=="m"))
dim(subset(whole_577_3185, miR.577=="expressed" & gender=="f"))
dim(subset(whole_577_3185, miR.577=="un-expressed" & gender=="m"))
dim(subset(whole_577_3185, miR.577=="un-expressed" & gender=="f"))

dim(subset(whole_577_3185, miR.577=="expressed" & batch=="1"))
dim(subset(whole_577_3185, miR.577=="expressed" & batch=="2"))
dim(subset(whole_577_3185, miR.577=="un-expressed" & batch=="1"))
dim(subset(whole_577_3185, miR.577=="un-expressed" & batch=="2"))


dim(subset(whole_577_3185, miR.3185=="expressed" & gender=="m"))
dim(subset(whole_577_3185, miR.3185=="expressed" & gender=="f"))
dim(subset(whole_577_3185, miR.3185=="un-expressed" & gender=="m"))
dim(subset(whole_577_3185, miR.3185=="un-expressed" & gender=="f"))

dim(subset(whole_577_3185, miR.3185=="expressed" & batch=="1"))
dim(subset(whole_577_3185, miR.3185=="expressed" & batch=="2"))
dim(subset(whole_577_3185, miR.3185=="un-expressed" & batch=="1"))
dim(subset(whole_577_3185, miR.3185=="un-expressed" & batch=="2"))


COVS_control=subset(COVS,group=="control")
COVS_HD=subset(COVS,group=="HD")
mean(COVS_control$visit_age)
sd(COVS_control$visit_age)
mean(COVS_HD$visit_age)
sd(COVS_HD$visit_age)
table(COVS_control$gender)
table(COVS_HD$gender)
mean(COVS_HD$cap_d_score)
sd(COVS_HD$cap_d_score)
mean(COVS_control$RIN)
sd(COVS_control$RIN)
mean(COVS_HD$RIN)
sd(COVS_HD$RIN)
range(COVS_control$education_max)
range(COVS_HD$education_max)





#install.packages('gsheet')
library(gsheet)
brain_normcount=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1ivkFQB-8Ltt8LTOz9t7d0hq2YH9k4PPcM69O2Uc8SdQ/edit#gid=48752243"))
rownames(brain_normcount)=brain_normcount$miRNA
dim(brain_normcount)
brain_normcount=brain_normcount[,-1]
dim(brain_normcount)
brain_20DE=brain_normcount[whole_DE_sva,]
brain_13DE=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1ivkFQB-8Ltt8LTOz9t7d0hq2YH9k4PPcM69O2Uc8SdQ/edit#gid=490463446"))

brain_result=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1ivkFQB-8Ltt8LTOz9t7d0hq2YH9k4PPcM69O2Uc8SdQ/edit#gid=348039719"))
rownames(brain_result)=brain_result$miRNA

DE_13=c("hsa-miR-1252-5p","hsa-miR-323a-3p","hsa-miR-1270","hsa-miR-154-5p","hsa-miR-135a-5p",
        "hsa-miR-1275","hsa-miR-138-5p","hsa-miR-146b-5p","hsa-miR-409-3p","hsa-miR-31-5p","hsa-miR-130a-3p",
        "hsa-miR-24-3p","hsa-miR-574-5p")
pp=round(brain_result[DE_13,]$adj.P.Val,3)
o=brain_13DE
k1=ggplot(o,aes_string("type",colnames(o)[1],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[1]))
k2=ggplot(o,aes_string("type",colnames(o)[2],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[2]))
k3=ggplot(o,aes_string("type",colnames(o)[3],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[3]))
k4=ggplot(o,aes_string("type",colnames(o)[4],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[4]))
k5=ggplot(o,aes_string("type",colnames(o)[5],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[5]))
k6=ggplot(o,aes_string("type",colnames(o)[6],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[6]))
k7=ggplot(o,aes_string("type",colnames(o)[7],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[7]))
k8=ggplot(o,aes_string("type",colnames(o)[8],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[8]))

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

k1=ggplot(o,aes_string("type",colnames(o)[9],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[9]))
k2=ggplot(o,aes_string("type",colnames(o)[10],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[10]))
k3=ggplot(o,aes_string("type",colnames(o)[11],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[11]))
k4=ggplot(o,aes_string("type",colnames(o)[12],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[12]))
k5=ggplot(o,aes_string("type",colnames(o)[13],fill="type"))+geom_boxplot()+theme(axis.title.x = element_blank(),legend.position = "none")+labs(title=paste0("pvalue=",pp[13]))


GO_genes=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1B_evBmFWSFqTVDyo3xWevxSBozcma8IqC5YklHL3M4g/edit#gid=812975757"))
write.csv(unique(GO_genes$genes.for.6),"cc.csv")
HTT_network=c("HIP1","CLTCL1","ITPR1","GRIN1","DLG4","TP53","CREBBP","REST","CLTC","AP2A2","HTT")

co_miRNA=intersect(rownames(EXPR),rownames(whole_data))
co_sample=intersect(colnames(EXPR),colnames(whole_data))
ff=EXPR[co_miRNA,co_sample]
gg=whole_data[co_miRNA,co_sample]
aaaaa=cor(ff,gg)
pheatmap(aaaaa,cluster_rows = FALSE,cluster_cols = FALSE)

