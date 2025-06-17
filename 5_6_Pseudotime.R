library(monocle3)
library(Matrix)
library(uwot)
library(SingleCellExperiment)
library(ggplot2)

EFIGA <- read.csv("extdata/labels_matching_uc_refined.csv")

meta_cols <- c("...1", "ID", "labels_refined", "labels_matching_UC_big")
metab_matrix <- EFIGA[, !names(EFIGA) %in% meta_cols]
rownames(EFIGA) <- EFIGA$ID
rownames(metab_matrix) <- EFIGA$ID

cell_metadata <- EFIGA[, c("labels_refined"), drop = FALSE]

expr_matrix <- t(as.matrix(metab_matrix))
mode(expr_matrix) <- "numeric"

zero_cols <- which(colSums(abs(expr_matrix), na.rm = TRUE) == 0)
if (length(zero_cols) > 0) {
  message("Removing ", length(zero_cols), " all-zero samples.")
  expr_matrix <- expr_matrix[, -zero_cols]
  cell_metadata <- cell_metadata[-zero_cols, , drop = FALSE]
}

cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = data.frame(
    gene_short_name = rownames(expr_matrix),
    row.names = rownames(expr_matrix)
  )
)



#Pseudotime_Refined
reducedDims(cds)$UMAP <- umap_res
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds)

png("output/pseudotime.png", width = 1200, height = 1000)
plot_cells(cds,
           cell_size = 4,
           show_trajectory_graph = TRUE,
           reduction_method = "UMAP",
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

dev.off()


#UMAP_Refined
png("output/Umap_ellipse_filled.png", width = 1200, height = 1000)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = label, color = label)) +
  stat_ellipse(
    aes(group = label),
    type       = "norm",
    level      = 0.85,
    geom       = "polygon",
    alpha      = 0.2,
    color      = NA
  ) +
  stat_ellipse(
    aes(group = label),
    type       = "norm",
    level      = 0.85,
    geom       = "path",
    linetype   = "dashed",
    size       = 1
  ) +
  geom_point(
    size  = 8,
    alpha = 1
  ) +
  scale_fill_manual(
    name = "Stage",
    values = c(
      "EMCI1" = "#F8766D",
      "EMCI2" = "#7CAE00",
      "LMCI1" = "#00BFC4",
      "LMCI2" = "#C77CFF"
    )
  ) +
  scale_color_manual(
    name = "Stage",
    values = c(
      "EMCI1" = "#F8766D",
      "EMCI2" = "#7CAE00",
      "LMCI1" = "#00BFC4",
      "LMCI2" = "#C77CFF"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position  = "right"
  ) +
  guides(
    fill  = guide_legend(override.aes = list(alpha = 0.5)),
    color = guide_legend(override.aes = list(size = 5))
  ) +
  xlab("UMAP1") +
  ylab("UMAP2")
dev.off()






#EMCI_LMCI
library(monocle3)
library(Matrix)
library(uwot)
library(SingleCellExperiment)
library(ggplot2)

EFIGA <- read.csv("extdata/labels_matching_uc_refined.csv")

meta_cols <- c("...1", "ID", "labels_refined", "labels_matching_UC_big")
metab_matrix <- EFIGA[, !names(EFIGA) %in% meta_cols]
rownames(EFIGA) <- EFIGA$ID
rownames(metab_matrix) <- EFIGA$ID

cell_metadata <- EFIGA[, c("labels_matching_UC_big"), drop = FALSE]

expr_matrix <- t(as.matrix(metab_matrix))
mode(expr_matrix) <- "numeric"

zero_cols <- which(colSums(abs(expr_matrix), na.rm = TRUE) == 0)
if (length(zero_cols) > 0) {
  message("Removing ", length(zero_cols), " all-zero samples.")
  expr_matrix <- expr_matrix[, -zero_cols]
  cell_metadata <- cell_metadata[-zero_cols, , drop = FALSE]
}

cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = data.frame(
    gene_short_name = rownames(expr_matrix),
    row.names = rownames(expr_matrix)
  )
)

set.seed(123)
umap_res <- uwot::umap(t(expr_matrix))
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$label <- cell_metadata$labels_matching_UC_big

#UMAP
png("output/Umap_ellipse_filled_big.png", width = 1200, height = 1000)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = label, color = label)) +
  # 1) Filled ellipse per label
  stat_ellipse(
    aes(group = label),
    type       = "norm",
    level      = 0.85,
    geom       = "polygon",
    alpha      = 0.2,
    color      = NA
  ) +
  # 2) Dashed outline around each ellipse
  stat_ellipse(
    aes(group = label),
    type       = "norm",
    level      = 0.85,
    geom       = "path",
    linetype   = "dashed",
    size       = 1
  ) +
  # 3) Points on top, colored by label
  geom_point(
    size  = 8,
    alpha = 1
  ) +
  scale_fill_manual(
    name = "Stage",
    values = c(
      "EMCI" = "#E41A1C", "LMCI" = "#377EB8"
    )
  ) +
  scale_color_manual(
    name = "Stage",
    values = c(
      "EMCI" = "#E41A1C", "LMCI" = "#377EB8"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position  = "right"
  ) +
  guides(
    fill  = guide_legend(override.aes = list(alpha = 0.5)),
    color = guide_legend(override.aes = list(size = 5))
  ) +
  xlab("UMAP1") +
  ylab("UMAP2")
dev.off()





#Pseudotime_EMCI_LMCI
reducedDims(cds)$UMAP <- umap_res
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds)

png("output/pseudotime_big.png", width = 1200, height = 1000)
plot_cells(cds,
           cell_size = 4,
           show_trajectory_graph = TRUE,
           reduction_method = "UMAP",
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE)

dev.off()