args<-commandArgs(TRUE);root<-normalizePath(args[1]);.libPaths(c(file.path(root,'code/R_library'),.libPaths()))
suppressPackageStartupMessages({library(MetaboAnalystR);library(qs2)})
# Compatibility only: current MetaboAnalyst databases use qs2 serialization.
# Preserve the original MetaboAnalyst local Globaltest/topology implementation.
ns<-asNamespace('MetaboAnalystR')
loader<-function(filenm,sub.dir=NULL){
 if(!file.exists(filenm)) download.file(paste0('https://www.metaboanalyst.ca/resources/libs/',if(!is.null(sub.dir))paste0(sub.dir,'/') else '',filenm),filenm,method='libcurl')
 tryCatch(qs::qread(filenm),error=function(err) qs2::qs_read(filenm))
}
assignInNamespace('.get.my.lib',loader,ns='MetaboAnalystR')
prep<-get('.prepare.qea.score',ns);b<-body(prep);patched<-FALSE
for(i in seq_along(b)) if(is.call(b[[i]])&&identical(b[[i]][[1]],as.name('if'))&&grepl('!.on.public.web',paste(deparse(b[[i]][[2]]),collapse=' '),fixed=TRUE)){b[[i]][[2]]<-FALSE;patched<-TRUE}
stopifnot(patched);body(prep)<-b
# Preserve matrix dimensions for single-compound matches.
prep<-eval(parse(text=gsub('mSetObj$dataSet$norm[, hit.inx]', 'mSetObj$dataSet$norm[, hit.inx, drop=FALSE]',paste(deparse(prep),collapse='\n'),fixed=TRUE)),envir=ns)

e<-read.csv(file.path(root,'inputs/EFIGA.csv'),row.names=1,check.names=FALSE)
d<-read.csv(file.path(root,'inputs/patients.csv'),row.names=1);d<-d[rownames(e),]
for(group in c('EMCI','LMCI')){
 out<-file.path(root,'results/metaboanalyst',group);dir.create(out,recursive=TRUE,showWarnings=FALSE);setwd(out)
 de<-read.csv(file.path(root,paste0('results/EFIGA_',group,'_within_DE.csv')))
 hits<-de$HMDB[de$P.Value<.05 & abs(de$logFC)>.3];write.csv(de[de$P.Value<.05 & abs(de$logFC)>.3,], 'DE_metabolites_raw_p005_effect03.csv',row.names=FALSE)
 ix<-d$label%in%paste0(group,1:2)
 x<-data.frame(Sample=rownames(e)[ix],Group=d$label[ix],e[ix,hits,drop=FALSE],check.names=FALSE)
 write.csv(x,'MetaboAnalyst_DE_abundances.csv',row.names=FALSE,quote=FALSE)
 stopifnot(length(hits)>0)
 mSet<-InitDataObjects('conc','pathqea',FALSE)
 mSet<-Read.TextData(mSet,file.path(out,'MetaboAnalyst_DE_abundances.csv'),'rowu','disc')
 mSet<-SanityCheckData(mSet)
 mSet<-ReplaceMin(mSet)
 mSet<-CrossReferencing(mSet,'hmdb');mSet<-CreateMappingResultTable(mSet)
 mSet<-PreparePrenormData(mSet)
 mSet<-Normalization(mSet,'NULL','NULL','NULL',ratio=FALSE,ratioNum=20)
 mSet<-SetKEGG.PathLib(mSet,'hsa','current');mSet<-SetMetabolomeFilter(mSet,FALSE)
 current.kegglib<-loader('hsa.qs','kegg/metpa')
 mSet<-prep(mSet,'rbc','gt')
 dat.in<-qs::qread('dat.in.qs')
 gt.obj<-globaltest::gt(dat.in$cls,dat.in$data,subsets=dat.in$subsets)
 dat.in$my.res<-globaltest::result(gt.obj)[,c(5,1),drop=FALSE]
 qs::qsave(dat.in,'dat.in.qs')
 mSet<-get('.save.qea.score',ns)(mSet)
 saveRDS(mSet,'mSet.rds')
 write.csv(mSet$analSet$qea.mat,'pathway_results_ids.csv')
 pathways <- mSet$analSet$qea.mat
 pathway_names <- setNames(names(current.kegglib$path.ids), unname(current.kegglib$path.ids))
 stopifnot(!is.null(pathway_names), all(rownames(pathways) %in% names(pathway_names)))
 rownames(pathways) <- unname(pathway_names[rownames(pathways)])
 write.csv(pathways,'pathway_results.csv')
 capture.output(sessionInfo(),file='sessionInfo.txt')
}
