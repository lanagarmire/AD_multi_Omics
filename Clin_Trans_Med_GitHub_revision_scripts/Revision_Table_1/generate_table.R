Sys.setlocale('LC_CTYPE', 'en_US.utf8')
args <- commandArgs(trailingOnly=FALSE)
here <- dirname(normalizePath(sub('^--file=', '', args[grepl('^--file=', args)])))
review <- dirname(here)
clinical_path <- file.path(here, 'inputs/cn_emci_lmci_ad_clin_pseudotime.csv')
merge_path <- file.path(here, 'inputs/ADNIMERGE.csv')
subtype_labels <- c(emcisubtype1='EMCI-1', emcisubtype2='EMCI-2', lmcisubtype1='LMCI-1', lmcisubtype2='LMCI-2')
labels <- c(CN='CN', AD='AD', subtype_labels)
clinical <- read.csv(clinical_path, stringsAsFactors=FALSE)
paired_path <- file.path(here, 'inputs/681p_172f.csv')
paired <- read.csv(paired_path, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(nrow(paired)==681, !anyDuplicated(paired$PID))
stopifnot(identical(as.integer(table(factor(paired$Label, levels=c('CN','EMCI','LMCI','AD')))), c(239L,201L,200L,41L)))
clinical <- clinical[clinical$PID %in% paired$PID, ]
stopifnot(setequal(clinical$PID, paired$PID))
write.csv(paired[,c('RID','PID','Label')], file.path(here,'paired_cohort_manifest.csv'), row.names=FALSE)
visits <- read.csv(merge_path, stringsAsFactors=FALSE, check.names=FALSE)
baseline <- visits[visits$VISCODE=='bl' & visits$PTID %in% clinical$PID,
                   c('PTID','AGE','PTGENDER','PTRACCAT','PTEDUCAT','APOE4')]
stopifnot(nrow(clinical)==681, !anyDuplicated(clinical$PID), !anyDuplicated(baseline$PTID))
d <- baseline[match(clinical$PID, baseline$PTID), ]
d$source <- 'Baseline ADNIMERGE'
stopifnot(identical(d$PTID, clinical$PID))
d$Group <- unname(labels[clinical$SNF])
stopifnot(identical(as.integer(table(factor(d$Group, levels=labels))), c(239L,41L,108L,93L,85L,115L)))
for (v in c('AGE','PTEDUCAT','APOE4')) d[[v]] <- as.numeric(d[[v]])
for (v in c('PTGENDER','PTRACCAT')) {
  d[[v]][trimws(d[[v]]) %in% c('', 'Unknown', 'NA')] <- NA
}
stopifnot(all(is.na(d$APOE4) | d$APOE4 %in% 0:2))
d$Carrier <- ifelse(is.na(d$APOE4), NA, ifelse(d$APOE4>0, 'Carrier', 'Noncarrier'))
# Cross-check numerical covariates against the finalized analysis clinical table.
audit <- data.frame(variable=c('AGE','PTEDUCAT','APOE4'),
                    n_discrepant=sapply(c('AGE','PTEDUCAT','APOE4'), function(v)
                      sum(abs(d[[v]]-as.numeric(clinical[[v]]))>1e-8, na.rm=TRUE)))
write.csv(audit, file.path(here,'covariate_source_audit.csv'), row.names=FALSE)
write.csv(d, file.path(here,'participant_characteristics_audit.csv'), row.names=FALSE, na='')
groups <- c('CN', 'All MCI', 'AD', unname(subtype_labels))
subsets <- lapply(groups, function(g) {
  if(g=='All MCI') d[d$Group %in% unname(subtype_labels),] else d[d$Group==g,]
})
counts <- sapply(subsets,nrow)
columns <- c('Characteristic', sprintf('%s (n=%d)',groups,counts), 'P: EMCI-1 vs EMCI-2', 'P: LMCI-1 vs LMCI-2')
rows <- list(); tests <- list()
pformat <- function(p) if(is.na(p)) '\u2014' else if(p<.001) '<0.001' else sprintf('%.3f',p)
test_variable <- function(v, categorical=FALSE) {
  sapply(c('EMCI','LMCI'), function(stage) {
    x <- d[startsWith(d$Group,stage) & !is.na(d[[v]]),]
    if(categorical) {
      tab <- table(x[[v]],x$Group)
      result <- if(nrow(tab)<2) list(p.value=NA_real_) else fisher.test(tab, workspace=2e7)
      method <- 'Two-sided Fisher exact test'
    } else {
      result <- wilcox.test(x[[v]] ~ x$Group, exact=FALSE)
      method <- 'Two-sided Wilcoxon rank-sum test'
    }
    tests[[length(tests)+1]] <<- data.frame(variable=v,stage=stage,n=nrow(x),method=method,p=result$p.value)
    result$p.value
  })
}
add <- function(name, values, p=c(NA,NA)) {
  rows[[length(rows)+1]] <<- c(name, values, sapply(p,pformat))
}
for(v in c('AGE','PTEDUCAT')) {
  title <- if(v=='AGE') 'Age, years, mean (SD)' else 'Education, years, mean (SD)'
  add(title, sapply(subsets,function(s) sprintf('%.1f (%.1f)',mean(s[[v]],na.rm=TRUE),sd(s[[v]],na.rm=TRUE))),test_variable(v))
  missing <- sapply(subsets,function(s) sum(is.na(s[[v]])))
  if(any(missing)) add(paste0(if(v=='AGE') 'Age' else 'Education', ', missing, n'), missing)
}
for(v in c('PTGENDER','PTRACCAT','Carrier')) {
  title <- switch(v, PTGENDER='Sex, n (%)', PTRACCAT='Race, n (%)', Carrier='APOE \u03b54 status, n (%)')
  add(title, rep('',length(groups)), test_variable(v,TRUE))
  categories <- if(v=='PTGENDER') c('Female','Male') else if(v=='Carrier') c('Carrier','Noncarrier') else sort(unique(na.omit(d[[v]])))
  stopifnot(all(na.omit(d[[v]]) %in% categories))
  for(category in categories) {
    add(paste0('  ',category),sapply(subsets,function(s) {
      n <- sum(s[[v]]==category,na.rm=TRUE); den <- sum(!is.na(s[[v]]))
      sprintf('%d (%.1f%%)',n,100*n/den)
    }))
  }
  missing <- sapply(subsets,function(s) sum(is.na(s[[v]])))
  if(any(missing)) add('  Missing, n',missing)
}
out <- as.data.frame(do.call(rbind,rows),stringsAsFactors=FALSE)
names(out) <- columns
write.csv(out,file.path(here,'ADNI_paired_omics_characteristics.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,tests),file.path(here,'comparison_tests.csv'),row.names=FALSE)
notes <- 'Values are mean (SD) or n (%). Percentages use participants with nonmissing values for the relevant characteristic. APOE \u03b54 carriers have one or two \u03b54 alleles. P values compare subtypes within EMCI and LMCI using two-sided Wilcoxon rank-sum tests for age and education and Fisher exact tests for sex, race, and APOE \u03b54 carrier status. P values are unadjusted. Characteristics are from baseline ADNIMERGE records; diagnoses and subtype assignments are from the finalized clinical analysis table. The cohort is restricted to the 681 participants in the paired gene-expression and metabolomics input. No missing values are imputed. All MCI contains the four MCI subtypes and is an aggregate column, not an additional group.'
table_title <- 'Table 1. Baseline characteristics of the ADNI paired-omics cohort by diagnosis and MCI subtype'
md <- c(paste('#',table_title),'',
         paste('|',paste(columns,collapse=' | '),'|'),
         paste('|',paste(rep('---',length(columns)),collapse=' | '),'|'),
         apply(out,1,function(r) paste('|',paste(r,collapse=' | '),'|')),'',notes)
writeLines(md,file.path(here,'ADNI_paired_omics_characteristics.md'),useBytes=TRUE)
html <- c('<!doctype html><meta charset="utf-8"><style>body{font-family:Arial;margin:30px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px 10px;border-bottom:1px solid #ddd;text-align:right}td:first-child,th:first-child{text-align:left}th{background:#eef2f5}p{max-width:1100px;font-size:12px}</style>',
 paste0('<h2>',table_title,'</h2><table><thead><tr>'),
 paste0('<th>',columns,'</th>'),'</tr></thead><tbody>',
 apply(out,1,function(r) paste0('<tr>',paste0('<td>',gsub('<','&lt;',r,fixed=TRUE),'</td>',collapse=''),'</tr>')),
 '</tbody></table>',paste0('<p>',notes,'</p>'))
writeLines(html,file.path(here,'ADNI_paired_omics_characteristics.html'),useBytes=TRUE)
print(out,row.names=FALSE)
print(audit,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(here,'R_session.txt'))
