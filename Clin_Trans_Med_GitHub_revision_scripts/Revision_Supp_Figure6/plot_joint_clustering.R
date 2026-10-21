# Plot the fixed final joint-clustering contingency table; no parameter search.
args <- commandArgs(trailingOnly=FALSE)
here <- dirname(normalizePath(sub('^--file=', '', args[grepl('^--file=', args)])))
tab_reordered <- as.matrix(read.csv(file.path(here,'inputs/contingency_table.csv'), row.names=1, check.names=FALSE))
stopifnot(identical(dim(tab_reordered), c(4L,4L)), sum(tab_reordered)==401)
stopifnot(identical(rownames(tab_reordered),c('EMCI1','EMCI2','LMCI1','LMCI2')),
          identical(colnames(tab_reordered),c('EMCI1','EMCI2','LMCI1','LMCI2')))
stopifnot(sum(diag(tab_reordered))==287)
library(ggplot2)
library(reshape2)
library(dplyr)

df <- melt(tab_reordered)
colnames(df) <- c("Unsupervised", "Inferred", "Count")

# Add original cluster id to unsupervised labels
df$Unsupervised <- recode(df$Unsupervised,
  "EMCI1" = "EMCI1\n(Cluster 4)",
  "EMCI2" = "EMCI2\n(Cluster 1)",
  "LMCI1" = "LMCI1\n(Cluster 3)",
  "LMCI2" = "LMCI2\n(Cluster 2)"
)

unsup_confusion_4grp = ggplot(df, aes(x = Inferred, y = Unsupervised, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Count), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Count") +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Unsupervised vs Inferred Subtypes",
    x     = "Inferred Subtypes (EMCI1/2, LMCI1/2)",
    y     = "Unsupervised 4 Clusters (SNF + Spectral)"
  ) +
  theme_minimal() +
  theme(
    axis.text      = element_text(size = 11),
    axis.title     = element_text(size = 12),
    plot.title     = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid     = element_blank(),
    legend.position = "right"
  )
unsup_confusion_4grp

png(file.path(here,'joint_clustering.png'),width = 7,height=7,res=300,unit='in')
unsup_confusion_4grp
dev.off()
