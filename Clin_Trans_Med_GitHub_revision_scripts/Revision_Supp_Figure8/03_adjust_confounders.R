#!/usr/bin/env Rscript
library(limma)
args <- commandArgs(trailingOnly=FALSE)
here <- dirname(normalizePath(sub('^--file=', '', args[grepl('^--file=', args)])))
src <- normalizePath(file.path(here,'inputs'))
out <- file.path(here,'results/additional_confounders');dir.create(out,recursive=TRUE,showWarnings=FALSE)
a <- read.csv(file.path(here,'results/additional_covariates.csv'),row.names=1)
stopifnot(!anyDuplicated(rownames(a)))
models <- list(baseline=character(),APOE4='APOE4_carrier',education='education',
 APOE4_education=c('APOE4_carrier','education'))
summaries <- list();diagnostics <- list();inventory <- list()
for(stage in c('EMCI','LMCI')) for(omics in c('MRNA','META')) {
 prefix<-paste(stage,omics,sep='_');lower<-tolower(stage)
 meta<-read.csv(file.path(here,'results',paste0(prefix,'_model_metadata.csv')),row.names=1)
 meta$group<-factor(meta$group,levels=paste0(lower,'subtype',1:2))
 meta$PTGENDER<-factor(meta$PTGENDER);meta$PTRACCAT<-factor(meta$PTRACCAT)
 meta$education<-a[rownames(meta),'PTEDUCAT']
 covars<-setdiff(names(meta),c('SNF','group','APOE4_source','APOE4_carrier','education'))
 extra<-unique(unlist(models))
 inventory[[prefix]]<-data.frame(analysis=prefix,covariate=extra,available_n=sapply(extra,function(v)sum(!is.na(meta[[v]]))),missing_n=sapply(extra,function(v)sum(is.na(meta[[v]]))))
 keep<-complete.cases(meta[,c('group',covars,extra)])
 write.csv(data.frame(PID=rownames(meta),included_common_sample=keep),file.path(out,paste0(prefix,'_sample_audit.csv')),row.names=FALSE)
 meta<-droplevels(meta[keep,,drop=FALSE])
 write.csv(meta,file.path(out,paste0(prefix,'_metadata.csv')))
 filename<-if(omics=='MRNA')paste0('mrna20k_',lower,'.csv')else paste0('all_',lower,'_metab_hmdb_lilikoi.csv')
 d<-read.csv(file.path(src,filename));ids<-if(omics=='MRNA')d$index else d$PID
 y<-t(as.matrix(d[,-c(1,2)]));colnames(y)<-ids;y<-y[,rownames(meta),drop=FALSE]
 cutoff<-if(omics=='META').3 else if(stage=='EMCI').35 else .4
 fits<-list()
 for(model in names(models)) {
  f<-reformulate(c('group',covars,models[[model]]));design<-model.matrix(f,meta)
  stopifnot(identical(rownames(design),colnames(y)))
  rank<-qr(design)$rank;condition<-kappa(design)
  status<-if(rank<ncol(design))'not_estimable_rank_deficient'else if(nrow(design)<=ncol(design))'not_estimable_insufficient_df'else 'ok'
  diagnostics[[paste(prefix,model)]]<-data.frame(analysis=prefix,model=model,n=nrow(design),parameters=ncol(design),rank=rank,condition_number_raw=condition,status=status,formula=paste(deparse(f),collapse=' '))
  if(status!='ok')next
  fit<-eBayes(lmFit(y,design));coef<-paste0('group',lower,'subtype2')
  tt<-topTable(fit,coef=coef,number=Inf,sort.by='none');tt$feature<-rownames(tt)
  tt$DE<-tt$adj.P.Val<.05&abs(tt$logFC)>cutoff
  write.csv(tt,file.path(out,paste0(prefix,'_',model,'_all_features.csv')),row.names=FALSE)
  write.csv(tt[tt$DE,],file.path(out,paste0(prefix,'_',model,'_DE.csv')),row.names=FALSE)
  write.csv(design,file.path(out,paste0(prefix,'_',model,'_design.csv')))
  fits[[model]]<-tt
 }
 for(model in names(fits)) {
  b<-fits$baseline;t<-fits[[model]];stopifnot(identical(b$feature,t$feature))
  summaries[[paste(prefix,model)]]<-data.frame(stage=stage,omics=omics,model=model,n=ncol(y),n_subtype1=sum(meta$group==levels(meta$group)[1]),n_subtype2=sum(meta$group==levels(meta$group)[2]),FDR05=sum(t$adj.P.Val<.05),DE=sum(t$DE),baseline_DE=sum(b$DE),retained=sum(b$DE&t$DE),retained_rate_percent=100*sum(b$DE&t$DE)/sum(b$DE),lost=sum(b$DE&!t$DE),gained=sum(!b$DE&t$DE),effect_correlation=cor(b$logFC,t$logFC),direction_agreement=mean(sign(b$logFC)==sign(t$logFC)))
  comp<-data.frame(feature=b$feature,effect_baseline=b$logFC,effect_adjusted=t$logFC,FDR_baseline=b$adj.P.Val,FDR_adjusted=t$adj.P.Val,DE_baseline=b$DE,DE_adjusted=t$DE)
  write.csv(comp,file.path(out,paste0(prefix,'_',model,'_comparison.csv')),row.names=FALSE)
 }
}
write.csv(do.call(rbind,summaries),file.path(out,'confounder_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,inventory),file.path(out,'covariate_inventory.csv'),row.names=FALSE)
inputs<-c(file.path(here,'results/additional_covariates.csv'),file.path(here,'03_adjust_confounders.R'))
write.csv(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_manifest.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(do.call(rbind,summaries))
