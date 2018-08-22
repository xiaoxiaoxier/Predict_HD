#> intersect(rownames(miRNA),rownames(HD_cap_pearson)[HD_cap_pearson>p])
#[1] "hsa-miR-125a-3p" "hsa-miR-133a"    "hsa-miR-326"     "hsa-miR-342-3p" 
#[5] "hsa-miR-410"     "hsa-miR-592"     "hsa-miR-885-5p"  "hsa-miR-1246"   
#> intersect(rownames(miRNA),rownames(HD_cap_pearson)[HD_cap_pearson<q])
#[1] "hsa-miR-339-3p" "hsa-miR-375"    "hsa-miR-627"    "hsa-miR-1180"   "hsa-miR-2054" 
rm(list=ls())
library(pheatmap)
library(gsheet)
library(fitdistrplus)
library(tidyr)
mRNA=read.table("ROSMAP.expression.postSVA.xls.renameCols",header = T,sep="\t")
rownames(mRNA) = mRNA[,1]
mRNA=mRNA[,-1]
mRNA=mRNA[,1:194]
for (i in 1:length(colnames(mRNA)))
{
  colnames(mRNA)[i]=unlist(strsplit(colnames(mRNA)[i], "X"))[2]
}


miRNA=read.table("ROSMAP_arraymiRNA.gct",sep="\t",header = T,stringsAsFactors=FALSE, skip=2)
rownames(miRNA) = miRNA[,1]
miRNA=miRNA[,-1]
miRNA=miRNA[,-1]

DE_13=c("hsa-miR-125a-3p","hsa-miR-133a","hsa-miR-326","hsa-miR-342-3p","hsa-miR-410",
        "hsa-miR-592","hsa-miR-885-5p","hsa-miR-1246","hsa-miR-339-3p",
        "hsa-miR-375","hsa-miR-627","hsa-miR-1180","hsa-miR-2054")
miRNA=miRNA[DE_13,]
for (i in 1:length(colnames(miRNA)))
{
  colnames(miRNA)[i]=unlist(strsplit(colnames(miRNA)[i], "[_]"))[2]
}


miRNA_cov=read.csv("ROSMAP_arraymiRNA_covariates.csv",header = T)
miRNA_cov$Sample_ID=as.character(miRNA_cov$Sample_ID)
for (i in 1:length(miRNA_cov$Sample_ID))
{
  miRNA_cov$Sample_ID[i]=unlist(strsplit(miRNA_cov$Sample_ID[i], "[_]"))[1]
}
sdf=unique(miRNA_cov$Sample_ID)
miRNA_cov=data.frame(projid=miRNA_cov$projid[match(sdf,miRNA_cov$Sample_ID)],Sample_ID=sdf)                                      
overlapsample=intersect(miRNA_cov$Sample_ID[match(colnames(mRNA),miRNA_cov$projid)],colnames(miRNA))
miRNA_overlap=miRNA[,overlapsample]
over_proj=miRNA_cov$projid[match(overlapsample,miRNA_cov$Sample_ID)]
mRNA_overlap=mRNA[,match(over_proj,colnames(mRNA))]

A=miRNA_overlap
A=apply(A, 2, as.numeric)
B=mRNA_overlap
ggggg=cor(t(A),t(B))
rownames(ggggg)=DE_13
colnames(ggggg)=rownames(mRNA)

#genes=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=0"))
#rownames(genes)=genes$ENSG

kkk=read.table("gencode.v24lift37.annotation.gtf.gene.txt",1,stringsAsFactors=FALSE)
kkk=separate(kkk,gene_id,into = c("gene_id","c"))

i=1
 df=as.numeric(ggggg[i,])
  fitn=fitdist(df,'norm')
  m=round(as.numeric(fitn$estimate[1]),digits=3)
  sd=round(as.numeric(fitn$estimate[2]),digits=3)
  p=round(qnorm(.05, mean=m, sd=sd, lower.tail = F), digits=3)
  q=round(qnorm(.95, mean=m, sd=sd, lower.tail = F), digits=3)
  
miR_410=c(names(ggggg[i,])[ggggg[i,]>p],names(ggggg[i,])[ggggg[i,]<q])
write.csv(miR_410,"jjj.csv")

miR_2054=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=1940840788"))
write.csv(kkk$gene_name[match(miR_433$cor_mRNA,kkk$gene_id)],"zzz.csv")

miR_125a_3p=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=642055842"))
miR_133a=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=504251739"))
miR_326=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=1549607819"))
miR_342_3p=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=956207299"))
miR_410=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=994470347"))
miR_592=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=783736962"))
miR_885_5p=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=0"))
miR_1246=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=1543933053"))
miR_339_3p=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=401904638"))
miR_375=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=2036528692"))
miR_627=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=566149628"))
miR_1180=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=427801054"))
miR_2054=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/1dF3_yMOGVU-mQOg0t7RFcBsw3xdo0wDJ1n9s7iS2PnA/edit#gid=1940840788"))
##############################
miR_138_cor=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=611314930"))
miR_433=data.frame(gsheet2tbl("https://docs.google.com/spreadsheets/d/19XWD6bLzzf-KZDcXBUR0BYG-Aya7dIcAsMmWNWgBIA8/edit#gid=696337344"))

################################

miR=miR_2054
mRNA_coexpressed = unique(miR$cor[!is.na(miR$cor)])
miRNA_target = unique(miR$target[!is.na(miR$target)])
DE_mRNA_paper = miR$DE
intersect(mRNA_coexpressed,miRNA_target)
intersect(intersect(mRNA_coexpressed,miRNA_target),DE_mRNA_paper)
