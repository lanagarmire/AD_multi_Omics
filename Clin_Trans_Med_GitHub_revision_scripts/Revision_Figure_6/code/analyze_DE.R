args<-commandArgs(TRUE);root<-args[1]
suppressPackageStartupMessages({library(limma);library(mgcv)})
e<-read.csv(file.path(root,'inputs/EFIGA.csv'),row.names=1,check.names=FALSE)
a<-read.csv(file.path(root,'inputs/ADNI.csv'),row.names=1,check.names=FALSE)
patient_file<-if(file.exists(file.path(root,'results/linear_monocle/patients.csv')))file.path(root,'results/linear_monocle/patients.csv') else file.path(root,'inputs/patients.csv')
d<-read.csv(patient_file,row.names=1)
al<-read.csv(file.path(root,'inputs/ADNI_labels.csv'),row.names=1)
stopifnot(nrow(e)==149,nrow(a)==401,identical(colnames(e),colnames(a)),setequal(rownames(e),rownames(d)))
d<-d[rownames(e),];al<-al[rownames(a),,drop=FALSE]
de<-function(x,label,base,comparison,prefix){
 keep<-label%in%c(base,comparison);g<-factor(label[keep],levels=c(base,comparison))
 design<-model.matrix(~g);fit<-eBayes(lmFit(t(as.matrix(x[keep,])),design))
 res<-topTable(fit,coef=2,number=Inf,sort.by='none',adjust.method='BH');res$HMDB<-rownames(res)
 res$baseline<-base;res$comparison<-comparison
 pooled<-sqrt(((sum(g==base)-1)*apply(x[keep,][g==base,],2,var)+(sum(g==comparison)-1)*apply(x[keep,][g==comparison,],2,var))/(sum(keep)-2))
 res$standardized_effect<-res$logFC/pooled
 write.csv(res,file.path(root,'results',paste0(prefix,'.csv')),row.names=FALSE)
}
de(e,substr(d$label,1,4),'EMCI','LMCI','EFIGA_broad_DE')
de(a,substr(al$label,1,4),'EMCI','LMCI','ADNI_broad_DE')
de(e,d$label,'EMCI1','EMCI2','EFIGA_EMCI_within_DE')
de(e,d$label,'LMCI1','LMCI2','EFIGA_LMCI_within_DE')
de(a,al$label,'EMCI1','EMCI2','ADNI_EMCI_within_DE')
de(a,al$label,'LMCI1','LMCI2','ADNI_LMCI_within_DE')
capture.output(sessionInfo(),file=file.path(root,'results/DE_R_sessionInfo.txt'))
cat('Six limma DE contrasts completed.\n')
