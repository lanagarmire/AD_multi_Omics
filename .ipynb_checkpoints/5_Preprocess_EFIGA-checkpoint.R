library(dplyr)
## Read in the data and annotation
metabolo<-read.csv("EFIGA/both_combined_merged_no_outliers_short.csv")
anno<- read.table("EFIGA/hilpos_annotations_level1-5.txt", header = TRUE, sep = "\t")




## extract metabolomics data and make LAB NY NUM become numeric
metabolo <- metabolo %>%
  mutate(LAB..NYNUM = as.numeric(sub("^..", "", LAB..NYNUM)))

meta_dat <- metabolo %>%
  select(matches("^hilpos_|^c18neg_.*\\d$"))
meta_dat <- meta_dat[, -1]


## Add a new column called matching_column_index, according to mz_time column, to match annotation rows to metabolomics data column names. The index indicates where is the column at this time point.
anno$matching_column_index <- NA

for (i in 1:nrow(anno)) {

  column_suffix <- gsub(".*_", "", anno$mz_time[i])

  matching_columns <- grep(paste0("_", column_suffix, "$"), colnames(meta_dat), value = TRUE)
  
  if (length(matching_columns) > 0) {

    matching_index <- which(colnames(meta_dat) == matching_columns[1])
    anno$matching_column_index[i] <- matching_index
  }
}
              
              
## Make annotation mz.theor column numeric, delete rows with NA. Add a column called difference to calculate the difference between actual mz and mz.theor. Then extract for confidence Level 1,3,5, what's the mean difference
## For level 5, rows_above_mean vector contains rows with difference larger than mean in Level 5. For example, If Level 5 has a mean difference 135, all rows in column difference where value>135 will be collected in rows_above_mean.
## Using the mathing index column before, we delete metabolites in rows_above_mean. So we deleted motabolites in Level 5 where has lagrer than mean difference in mz and mz theory.
anno$mz.theor<-as.numeric(anno$mz.theor)
anno <- anno[complete.cases(anno$mz.theor), ]

table(anno$confidence)

anno$difference <- abs(anno$mz - anno$mz.theor)
confidence_anno <- anno[anno$confidence %in% c("Level 1", "Level 3", "Level 5"), ]

average_diff <- aggregate(difference ~ confidence, data = confidence_anno, mean)
average_diff
level_5_mean <-average_diff$difference[average_diff$confidence == "Level 5"]
rows_above_mean <-which(anno$confidence == "Level 5" & anno$difference > level_5_mean)
dim(meta_dat)
columns_to_remove <- unique(anno$matching_column_index[rows_above_mean])
columns_to_remove <- na.omit(columns_to_remove)

if (length(columns_to_remove) > 0) {
  # Ensure column indices are within valid range
  columns_to_remove <- columns_to_remove[columns_to_remove <= ncol(meta_dat)]
  
  if (length(columns_to_remove) > 0) {
    meta_dat <- meta_dat[, -columns_to_remove]
    cat("How many columns being removed:\n")
    cat(sprintf("Number of columns removed: %d\n", length(columns_to_remove)))
  } else {
    cat("No valid columns were removed.\n")
  }
} else {
  cat("No columns were removed.\n")
}
dim(meta_dat)

              
              

## Compute correction
#calculate average for each column 
global_avg <- colMeans(meta_dat, na.rm = TRUE)
meta_dat <- sweep(meta_dat, 2, global_avg, FUN = "/")


## Compute CV (Coefficient of Variation) for all columns
cv.all <- data.frame()
for (col_name in names(meta_dat)) {
  col <- meta_dat[[col_name]] 
  
  if (is.numeric(col)) {  
    cv_value <- sd(col, na.rm = TRUE) / mean(col, na.rm = TRUE)
  } else {
    cv_value <- NA  # Set CV to NA for non-numeric columns
  }

  cv_result <- data.frame(Column = col_name, CV = cv_value)

  cv.all <- rbind(cv.all, cv_result)
}


              
              
## Remove metabolites with CV > 0.2 
meta_dat <- meta_dat[, !(names(meta_dat) %in% cv.all$Column[cv.all$CV > 0.2])]
              

#meta_dat <- meta_dat[1:10, ]
## Check mets with non-numeric entries
nancols <- colSums(is.na(meta_dat))
rm_cols <- 1/nrow(meta_dat) * nancols > 0.40
cols2rm <- colnames(meta_dat)[rm_cols]

if (sum(rm_cols) > 0) {
  cat(sprintf('The following %d columns were flagged for removal due to missingness:\n', sum(rm_cols)))
  for (i in 1:sum(rm_cols)) {
    cat(sprintf('>%s\n', cols2rm[i]))
  }
} else {
  cat(sprintf('No columns were removed due to missingness.\n'))
}

meta_dat <- meta_dat[, !colnames(meta_dat) %in% cols2rm]


boxplot(meta_dat)
print(dim(meta_dat))

#metabolo<- metabolo[1:10, ]              
## Check subjects with non-numeric entries
nansubs <- rowSums(is.na(meta_dat))
zerosubs <- rowSums(abs(meta_dat) < 0.00001, na.rm = TRUE)
rm_subs <- 1/ncol(meta_dat) * (nansubs + zerosubs) > 0.40


meta_dat <- cbind(projid = metabolo$LAB..NYNUM, meta_dat)


projids2discard <- meta_dat$projid[rm_subs]


if (sum(rm_subs) > 0) {
  cat(sprintf('The following %d subjects were removed due to missingness:\n', length(projids2discard)))
  for (i in 1:sum(rm_subs)) {
    cat(sprintf('>%d\n', projids2discard[i]))
  }
} else {
  cat(sprintf('No subjects were removed due to missingness.\n'))
}

meta_dat <- meta_dat[!rm_subs, ]
              
              
## Consolidate biological replicates

#find all duplicates,if find duplicate, replace it with mean value between the rows.
if (length(unique(meta_dat$projid)) == nrow(meta_dat)) {
  cat(sprintf('There are no replicates.\n'))
} else {
  u.projid <- sort(unique(meta_dat$projid))
  counts <- as.numeric(table(meta_dat$projid))
  dup <- unique(meta_dat$projid[meta_dat$projid %in% u.projid[counts > 1]])
  
  cat(sprintf('Replicates were consolidated for the following PROJIDs:\n'))
  for (i in 1:length(dup)) {
    cat(sprintf('>%d\n', dup[i]))
    projid.index <- which(meta_dat$projid == dup[i])
    dat.slice <- meta_dat[projid.index, ]
    meta_dat[projid.index, -1] <- sapply(colMeans(dat.slice[, -1], na.rm = TRUE), rep, length(projid.index))
    meta_dat <- meta_dat[-projid.index[2:length(projid.index)], ]
  }
}

## Replace 0 or NaN LODs with 0.5*column min
for (i in 1:nrow(meta_dat)) {
  for (j in 1:ncol(meta_dat)) {
    if ((meta_dat[i, j] == 0) || (is.na(meta_dat[i, j]))) {
      meta_dat[i, j] <- 0.5*min(meta_dat[meta_dat[ , j] != 0, j], na.rm = TRUE)
      print(meta_dat[i, j])
    }
  }
}
               



id<-meta_dat[, 1, drop = FALSE]
meta_dat<- meta_dat[, -1]
## Center based on sample means
meta_dat <- log2(meta_dat)
meta_dat <- scale(meta_dat, scale = TRUE)

## Replace anything greater than >3 sd away from center
meta_dat[meta_dat < -3] <- -3
meta_dat[meta_dat > 3] <- 3

boxplot(meta_dat)
dim(meta_dat)

meta_dat <- cbind(id, meta_dat)


write.csv(meta_dat, 'outputs/LQ_preprocessed_uplc.csv', row.names = FALSE)










rosmap_metab <- meta_dat
rosmap_metab<-as.data.frame(rosmap_metab)
anno<- read.table("EFIGA/hilpos_annotations_level1-5.txt", header = TRUE, sep = "\t")

colnames_rosmap_metab <- colnames(rosmap_metab)
suffix_rosmap <- sapply(strsplit(colnames_rosmap_metab, "_"), function(x) paste(x[2], x[3], sep = "_"))
mz_to_hmdb <- setNames(anno$HMDB_ID, anno$mz_time)
new_colnames <- sapply(seq_along(suffix_rosmap), function(i) {
  suffix <- suffix_rosmap[i]
  if (suffix %in% names(mz_to_hmdb)) {
    mz_to_hmdb[[suffix]]
  } else {
    colnames_rosmap_metab[i] 
  }
})
colnames(rosmap_metab) <- new_colnames



## Load pre-processed ADNI metabolomics data
adni_emci_metab <- read.csv("extdata/all_emci_metab_hmdb_lilikoi.csv", row.names = 1)  # EMCI metabolomics
adni_lmci_metab <- read.csv("extdata/all_lmci_metab_hmdb_lilikoi.csv", row.names = 1)  # LMCI metabolomics
adni_emci_metab <- adni_emci_metab[,2:ncol(adni_emci_metab)]
adni_lmci_metab <- adni_lmci_metab[,2:ncol(adni_lmci_metab)]
head(adni_emci_metab)
head(adni_lmci_metab)


## Label and merge ADNI metabolomics data
subtype_info<-read.csv("extdata/cn_emci_lmci_ad_clin_pseudotime_snf3.csv",header=TRUE)
table(subtype_info$SNF)
name_emci1<-subtype_info$PID[which(subtype_info$SNF3=="emcisubtype1")]
name_emci2<-subtype_info$PID[which(subtype_info$SNF3=="emcisubtype2")]
name_lmci1<-subtype_info$PID[which(subtype_info$SNF3=="lmcisubtype1")]
name_lmci2<-subtype_info$PID[which(subtype_info$SNF3=="lmcisubtype2")]

adni_emci_metab$subtype=rep("EMCI1",nrow(adni_emci_metab))
adni_emci_metab$subtype[which(rownames(adni_emci_metab) %in% name_emci2)]="EMCI2"

adni_lmci_metab$subtype=rep("LMCI1",nrow(adni_lmci_metab))
adni_lmci_metab$subtype[which(rownames(adni_lmci_metab) %in% name_lmci2)]="LMCI2"

adni_metab_merge = rbind(adni_emci_metab,adni_lmci_metab)
table(adni_metab_merge$subtype)




## Load Differential Analysis Result ADNI metabolomics data
New_adj_DEmetab_LMCI <- read.csv("extdata/New_adj_DEmetab_LMCI.csv", row.names = 1)
New_adj_DEmetab_EMCI <- read.csv("extdata/New_adj_DEmetab_EMCI.csv", row.names = 1)
new_markers = union(row.names(New_adj_DEmetab_LMCI), row.names(New_adj_DEmetab_EMCI))


## Find overlapping metabolites between ROSMAP and ADNI
intersect_metab=intersect(new_markers,colnames(rosmap_metab)) 
intersect_adni=subset(adni_metab_merge, select=intersect_metab)
intersect_rosmap=subset(rosmap_metab, select=intersect_metab)
length(intersect_metab)
rosmap_metab_mci=rosmap_metab
intersect_rosmap_mci=subset(rosmap_metab_mci, select=intersect_metab)


## Map EFIGA data to ADNI data  
target_vector=normalize.quantiles.determine.target(as.matrix(intersect_adni))
norm_rosmap_mci=normalize.quantiles.use.target(as.matrix(intersect_rosmap_mci),target=target_vector)
colnames(norm_rosmap_mci)<-colnames(intersect_rosmap_mci)
rownames(norm_rosmap_mci)<-rownames(intersect_rosmap_mci)
norm_rosmap_mci<-data.frame(norm_rosmap_mci)


write.csv(intersect_adni, "extdata/intersect_adni_metab.csv")
write.csv(norm_rosmap_mci, "extdata/intersect_rosmap_mci_metab.csv")
