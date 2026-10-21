#!/usr/bin/env Rscript
# Transcriptomic baseline: finalized Revision_Supp_Fig4CDEF_volcano.ipynb,
# after cell-type adjustment for BOTH stages. Metabolomics: original META notebooks.
library(limma)
args <- commandArgs(trailingOnly=FALSE)
here <- dirname(normalizePath(sub('^--file=', '', args[grepl('^--file=', args)])))
src <- normalizePath(file.path(here, 'inputs'))
out <- file.path(here, 'results')
dir.create(out, showWarnings=FALSE)
clinical <- read.csv(file.path(src, 'cn_emci_lmci_ad_clin_pseudotime.csv'), row.names=1)
stopifnot(!anyDuplicated(rownames(clinical)))
visits_path <- normalizePath(file.path(here, 'inputs/ADNIMERGE.csv'))
visits <- read.csv(visits_path, stringsAsFactors=FALSE)
apoe <- lapply(split(visits$APOE4, visits$PTID), function(v) unique(v[!is.na(v)]))
stopifnot(all(lengths(apoe) <= 1), all(unlist(apoe) %in% 0:2))
apoe <- vapply(apoe, function(v) if(length(v)) v else NA_real_, numeric(1))
clinical$APOE4_source <- apoe[rownames(clinical)]
known <- !is.na(clinical$APOE4) & !is.na(clinical$APOE4_source)
stopifnot(all(clinical$APOE4[known] == clinical$APOE4_source[known]))
clinical$APOE4_carrier <- ifelse(is.na(clinical$APOE4_source), NA, as.integer(clinical$APOE4_source >= 1))
cells <- c('T.cells.CD8','T.cells.CD4.naive','T.cells.CD4.memory.activated',
 'T.cells.CD4.memory.resting','Neutrophils','Eosinophils','Monocytes',
 'NK.cells.resting','Macrophages.M2')
summaries <- list(); inputs <- c(visits_path, file.path(src, 'cn_emci_lmci_ad_clin_pseudotime.csv'))
for(stage in c('EMCI','LMCI')) for(omics in c('MRNA','META')) {
  lower <- tolower(stage); prefix <- paste(stage, omics, sep='_')
  filename <- if(omics=='MRNA') paste0('mrna20k_',lower,'.csv') else paste0('all_',lower,'_metab_hmdb_lilikoi.csv')
  path <- file.path(src, filename); inputs <- c(inputs, path)
  dat <- read.csv(path, check.names=TRUE)
  ids <- if(omics=='MRNA') dat$index else dat$PID
  stopifnot(!anyDuplicated(ids))
  y <- t(as.matrix(dat[, -c(1,2)])); colnames(y) <- ids
  stopifnot(is.numeric(y), !anyDuplicated(rownames(y)), all(is.finite(y)))
  meta <- clinical[match(ids, rownames(clinical)),,drop=FALSE]; rownames(meta) <- ids
  levels <- paste0(lower,'subtype',1:2)
  meta$group <- factor(meta$SNF, levels=levels)
  meta$PTGENDER <- factor(as.character(meta$PTGENDER))
  meta$PTRACCAT <- factor(as.character(meta$PTRACCAT))
  covars <- c('AGE','PTGENDER','PTRACCAT')
  if(omics=='MRNA') {
    path <- file.path(src,paste0('mrna_',lower,'_deconv.tsv')); inputs <- c(inputs,path)
    dc <- read.delim(path, row.names=1)
    stopifnot(!anyDuplicated(rownames(dc)))
    # Original merge limits both transcriptomic models to deconvolution subjects.
    eligible <- ids %in% rownames(dc)
    {
      for(cell in cells) meta[[cell]] <- dc[ids,cell]
      covars <- c(cells,covars)
    }
  } else eligible <- rep(TRUE,length(ids))
  keep <- eligible & complete.cases(meta[,c('group',covars,'APOE4_carrier')])
  audit <- data.frame(PID=ids, SNF=meta$SNF, APOE4=meta$APOE4_source,
    APOE4_carrier=meta$APOE4_carrier, deconvolution_eligible=eligible, included=keep)
  write.csv(audit,file.path(out,paste0(prefix,'_participant_audit.csv')),row.names=FALSE)
  meta <- droplevels(meta[keep,,drop=FALSE]); y <- y[,keep,drop=FALSE]
  stopifnot(nlevels(meta$group)==2, length(unique(meta$APOE4_carrier))==2)
  if(omics=='MRNA') for(cell in cells) meta[[cell]] <- as.numeric(scale(meta[[cell]],center=FALSE))
  write.csv(meta[,c('SNF','group',covars,'APOE4_source','APOE4_carrier')],file.path(out,paste0(prefix,'_model_metadata.csv')))
  cutoff <- if(omics=='META') 0.3 else if(stage=='EMCI') 0.35 else 0.4
  fits <- list(); formulas <- c()
  for(model in c('baseline','APOE4_adjusted')) {
    formula <- reformulate(c('group',covars,if(model=='APOE4_adjusted') 'APOE4_carrier'))
    design <- model.matrix(formula, meta)
    stopifnot(identical(rownames(design),colnames(y)), qr(design)$rank==ncol(design),nrow(design)>ncol(design))
    coef <- paste0('group',levels[2]); stopifnot(coef %in% colnames(design))
    fit <- eBayes(lmFit(y,design))
    tt <- topTable(fit,coef=coef,number=Inf,sort.by='none',adjust.method='BH')
    tt$feature <- rownames(tt)
    tt$DE <- tt$adj.P.Val < 0.05 & abs(tt$logFC)>cutoff
    tt$direction <- ifelse(tt$DE,ifelse(tt$logFC>0,'UP','DOWN'),'NO')
    write.csv(tt,file.path(out,paste0(prefix,'_',model,'_all_features.csv')),row.names=FALSE)
    write.csv(tt[tt$DE,],file.path(out,paste0(prefix,'_',model,'_DE.csv')),row.names=FALSE)
    write.csv(design,file.path(out,paste0(prefix,'_',model,'_design.csv')))
    fits[[model]] <- tt; formulas[model] <- paste(deparse(formula),collapse=' ')
    png(file.path(out,paste0(prefix,'_',model,'_volcano.png')),width=1500,height=1200,res=180)
    plot(tt$logFC,-log10(pmax(tt$adj.P.Val,.Machine$double.xmin)),pch=16,cex=.5,
      col=ifelse(tt$direction=='UP','red',ifelse(tt$direction=='DOWN','blue','grey')),
      xlab='Subtype 2 - subtype 1 (input-scale coefficient)',ylab='-log10(BH FDR)',
      main=paste(prefix,model,'DE =',sum(tt$DE)))
    abline(v=c(-cutoff,cutoff),h=-log10(.05),lty=2); dev.off()
  }
  b <- fits$baseline; a <- fits$APOE4_adjusted
  stopifnot(identical(b$feature,a$feature))
  comparison <- data.frame(feature=b$feature, effect_baseline=b$logFC,effect_APOE4=a$logFC,
    effect_change=a$logFC-b$logFC,FDR_baseline=b$adj.P.Val,FDR_APOE4=a$adj.P.Val,
    DE_baseline=b$DE,DE_APOE4=a$DE,same_direction=sign(b$logFC)==sign(a$logFC))
  write.csv(comparison,file.path(out,paste0(prefix,'_comparison.csv')),row.names=FALSE)
  summaries[[prefix]] <- data.frame(stage=stage,omics=omics,n=ncol(y),n_subtype1=sum(meta$group==levels[1]),
    n_subtype2=sum(meta$group==levels[2]),features=nrow(y),cutoff=cutoff,
    baseline_formula=formulas['baseline'],adjusted_formula=formulas['APOE4_adjusted'],
    baseline_FDR05=sum(b$adj.P.Val<.05),adjusted_FDR05=sum(a$adj.P.Val<.05),
    baseline_DE=sum(b$DE),adjusted_DE=sum(a$DE),retained_DE=sum(b$DE & a$DE),
    lost_DE=sum(b$DE & !a$DE),gained_DE=sum(!b$DE & a$DE),
    effect_correlation=cor(b$logFC,a$logFC),direction_agreement=mean(comparison$same_direction))
}
summary <- do.call(rbind,summaries)
write.csv(summary,file.path(out,'sensitivity_summary.csv'),row.names=FALSE)
inputs <- unique(c(inputs, file.path(here,'01_baseline_APOE.R')))
write.csv(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(out,'input_manifest.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary)
