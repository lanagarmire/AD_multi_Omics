# Original Figure6E: group-defined subsets, GAM y ~ s(pt,k=5), REML.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
library(mgcv)
patients <- read.csv(file.path(root,'inputs/patients.csv'),row.names=1,check.names=FALSE)
groups <- c('EMCI2','LMCI1','LMCI2')
for(panel in c('shared54')) {
 x <- read.csv(file.path(root,'inputs',paste0(panel,'.csv')),row.names=1,check.names=FALSE)
 x <- x[rownames(patients),,drop=FALSE]
 rows <- list(); rownum <- 0
 for(group in groups) {
  keep <- patients$label %in% c('EMCI1',group) & is.finite(patients$pseudotime)
  pt <- patients$pseudotime[keep]
  pt <- (pt-min(pt))/diff(range(pt))
  for(feature in colnames(x)) {
   fit <- gam(y ~ s(pt,k=5),data=data.frame(y=x[keep,feature],pt=pt),method='REML')
   ss <- summary(fit)
   rownum <- rownum+1
   rows[[rownum]] <- data.frame(feature=feature,comparison=paste0('EMCI1->',group),R2_adjusted=ss$r.sq,deviance_explained=ss$dev.expl,p_raw=ss$s.table[1,'p-value'],n_patients=sum(keep),n_reference=sum(patients$label[keep]=='EMCI1'),n_comparison=sum(patients$label[keep]==group),edf=ss$s.table[1,'edf'])
  }
  cat(panel,group,'complete\n');flush.console()
 }
 result <- do.call(rbind,rows)
 # Match original: one BH family across all metabolite/subset pairs per panel.
 result$q_BH <- p.adjust(result$p_raw,method='BH')
 write.csv(result,file.path(root,'results',paste0(panel,'_original_style_pair_GAM.csv')),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(root,'results/pair_GAM_sessionInfo.txt'))
