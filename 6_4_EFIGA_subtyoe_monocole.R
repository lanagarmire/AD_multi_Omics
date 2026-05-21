library(monocle3)
library(Matrix)
library(uwot)
library(SingleCellExperiment)
library(ggplot2)

EFIGA <- read.csv("/nfs/turbo/umms-lgarmire/home/yhdu/Bowei_NAS/ADNI_Tidy/Code/Figure56_Validation/Labels/l1.csv")
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

set.seed(123)
umap_res <- uwot::umap(t(expr_matrix))
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$label <- cell_metadata$labels_refined


png('checkroot.png',width=5,height=5,unit='in',res=300)
plot_cells(
  cds,
  color_cells_by = "labels_refined",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 3
)
dev.off()

reducedDims(cds)$UMAP <- umap_res
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE, learn_graph_control = list(
  ncenter= 60,prune_graph=TRUE,minimal_branch_len=15))
cds <- order_cells(cds,root_pr_nodes = "Y_14")

## 1. Extract principal graph node coordinates
coords <- as.data.frame(t(cds@principal_graph_aux[["UMAP"]]$dp_mst))
colnames(coords)[1:2] <- c("UMAP1", "UMAP2")
coords$node <- rownames(coords)
coords$node_number <- gsub("Y_", "", coords$node)

## 2. Extract graph
p <- principal_graph(cds)[["UMAP"]]

## 3. Extract selected root node from order_cells()
root_node <- cds@principal_graph_aux[["UMAP"]]$root_pr_nodes
root_node
# should return something like "Y_18"

## 4. Calculate node distance from selected root
node_dist <- igraph::distances(
  p,
  v = root_node,
  to = igraph::V(p)$name,
  weights = igraph::E(p)$weight
)

node_order_df <- data.frame(
  node = colnames(node_dist),
  distance_from_root = as.numeric(node_dist[1, ])
)

node_order_df <- node_order_df[order(node_order_df$distance_from_root), ]
node_order_df$order_from_root <- seq_len(nrow(node_order_df))

coords <- merge(coords, node_order_df, by = "node", all.x = TRUE)

## 5. Mark the selected root node
coords$is_root <- coords$node %in% root_node

## 6. Plot pseudotime and overlay numbered principal nodes
g <- plot_cells(
  cds,
  cell_size = 2,
  show_trajectory_graph = TRUE,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

g +
  geom_point(
    data = coords,
    aes(x = UMAP1, y = UMAP2),
    inherit.aes = FALSE,
    shape = 21,
    size = 5,
    stroke = 1.2,
    fill = "white",
    color = "black"
  ) +
  geom_point(
    data = subset(coords, is_root),
    aes(x = UMAP1, y = UMAP2),
    inherit.aes = FALSE,
    shape = 21,
    size = 6,
    stroke = 1.5,
    fill = "white",
    color = "black"
  ) +
  geom_text(
    data = coords,
    aes(x = UMAP1, y = UMAP2, label = order_from_root),
    inherit.aes = FALSE,
    size = 3,
    color = "black"
  )




g<- plot_cells(cds,
               cell_size = 2,
               show_trajectory_graph = TRUE,
               reduction_method = "UMAP",
               color_cells_by = "pseudotime",
               label_groups_by_cluster = FALSE,
               label_leaves = TRUE,
               label_branch_points = TRUE)
pseudotime_plot = g+
  geom_point(
    data = coords,
    aes(x=UMAP1, y = UMAP2),
    inherit.aes = FALSE,
    shape=21,
    size = 5,
    stroke =1.2,
    fill = "white",
    color = "black"
  )+
  geom_text(
    data = coords,
    aes ( x= UMAP1 ,y= UMAP2, label = order_from_root),
    inherit.aes = FALSE ,
    size =3,
    color = "black"
  )
pseudotime_plot

png('/nfs/turbo/umms-lgarmire/home/yhdu/Bowei_NAS/EFIGA/Review/AD_RnT/pseudotime_plot.png',
    width=7.5,height=6,res=300,unit='in')
pseudotime_plot
dev.off()





#trajectory on UMAP
edges_df  <- igraph::as_data_frame(p, what = "edges")

seg_data <- merge(edges_df, coords[, c("node", "UMAP1", "UMAP2")],
                  by.x = "from", by.y = "node")
names(seg_data)[names(seg_data) == "UMAP1"] <- "x"
names(seg_data)[names(seg_data) == "UMAP2"] <- "y"

seg_data <- merge(seg_data, coords[, c("node", "UMAP1", "UMAP2")],
                  by.x = "to", by.y = "node")
names(seg_data)[names(seg_data) == "UMAP1"] <- "xend"
names(seg_data)[names(seg_data) == "UMAP2"] <- "yend"

# Now use seg_data in the plot
final_umap_with_node_trj = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = label, color = label)) +
  stat_ellipse(aes(group = label), type = "norm", level = 0.95,
               geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(group = label), type = "norm", level = 0.95,
               geom = "path", linetype = "dashed", linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  scale_fill_manual(name = "Stage",
                    values = c("EMCI1"="#F8766D","EMCI2"="#7CAE00",
                               "LMCI1"="#00BFC4","LMCI2"="#C77CFF")) +
  scale_color_manual(name = "Stage",
                     values = c("EMCI1"="#F8766D","EMCI2"="#7CAE00",
                                "LMCI1"="#00BFC4","LMCI2"="#C77CFF")) +
  geom_segment(data = seg_data,
               aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, color = "black", linewidth = 0.8) +
  geom_point(data = coords, aes(x = UMAP1, y = UMAP2),
             inherit.aes = FALSE, shape = 21, size = 5, stroke = 1.2,
             fill = "white", color = "black") +
  geom_text(data = coords, aes(x = UMAP1, y = UMAP2, label = order_from_root),
            inherit.aes = FALSE, size = 3, color = "black") +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "right") +
  guides(fill  = guide_legend(override.aes = list(alpha = 0.5)),
         color = guide_legend(override.aes = list(size = 4))) +
  xlab("UMAP1") + ylab("UMAP2")


final_umap_with_node_trj

png('/nfs/turbo/umms-lgarmire/home/yhdu/Bowei_NAS/EFIGA/Review/AD_RnT/final_umap_with_node_trj.png',
    width=8,height=7,res=300,unit='in')
final_umap_with_node_trj
dev.off()



library(ggplot2)

# Extract pseudotime
pt <- monocle3::pseudotime(cds)

# Build dataframe
df <- data.frame(
  sample = names(pt),
  pseudotime = as.numeric(pt),
  labels_refined = colData(cds)$labels_refined
)

# Optional ordering by median pseudotime
library(dplyr)

label_order <- df %>%
  group_by(labels_refined) %>%
  summarise(median_pt = median(pseudotime, na.rm = TRUE)) %>%
  arrange(median_pt) %>%
  pull(labels_refined)

df$labels_refined <- factor(df$labels_refined,
                            levels = label_order)

# Boxplot
plot_time_by_group = ggplot(df, aes(x = labels_refined,
               y = pseudotime,
               fill = labels_refined)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15,
              size = 1,
              alpha = 0.5) +
  theme_bw(base_size = 14) +
  labs(
    x = "Refined labels",
    y = "Pseudotime",
    title = "Pseudotime distribution across refined labels"
  ) +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    legend.position = "none"
  )

png('plot_time_by_group.png',res=300,width=5,height=5,unit='in')
plot_time_by_group
dev.off()

##############################
# ── Section 8: Branch-level regression (group-defined branches) ────────────────
library(dplyr)

R2_THRESHOLD <- 0.3

pt_vec     <- monocle3::pseudotime(cds)
grp_vec    <- colData(cds)$labels_refined   # EMCI1/EMCI2/LMCI1/LMCI2
finite_idx <- is.finite(pt_vec)

pt_all   <- pt_vec[finite_idx]
grp_all  <- grp_vec[finite_idx]
expr_all <- expr_matrix[, finite_idx]

# EMCI1 = root; each other group = a branch endpoint
branches <- list(
  "EMCI1->EMCI2" = c("EMCI1", "EMCI2"),
  "EMCI1->LMCI1" = c("EMCI1", "LMCI1"),
  "EMCI1->LMCI2" = c("EMCI1", "LMCI2")
)

library(mgcv)

run_branch_regression <- function(branch_name, groups) {
  idx         <- grp_all %in% groups
  pt_branch   <- pt_all[idx]
  expr_branch <- expr_all[, idx, drop = FALSE]
  
  # Normalize pseudotime to [0,1] per branch (matches ClinicTraj behavior)
  pt_scaled <- (pt_branch - min(pt_branch)) / (max(pt_branch) - min(pt_branch))
  
  lapply(rownames(expr_branch), function(metab) {
    y   <- as.numeric(expr_branch[metab, ])
    df  <- data.frame(y = y, pt = pt_scaled)
    
    # GAM with thin-plate spline — non-linear, equivalent to GPR
    fit  <- mgcv::gam(y ~ s(pt, k = 5), data = df, method = "REML")
    s    <- summary(fit)
    
    # Compute R² manually from deviance (gam summary R² = dev.expl)
    data.frame(
      branch     = branch_name,
      metabolite = metab,
      R2         = s$r.sq,          # adjusted R² from GAM
      dev_expl   = s$dev.expl,      # deviance explained (analogous to R²)
      pvalue     = s$s.table[1, "p-value"],
      n_samples  = sum(idx),
      stringsAsFactors = FALSE
    )
  }) |> do.call(what = rbind)
}

branch_regression_df <- mapply(
  run_branch_regression,
  names(branches), branches,
  SIMPLIFY = FALSE
) |> do.call(what = rbind)

branch_regression_df$padj <- p.adjust(branch_regression_df$pvalue, method = "BH")
branch_regression_df <- branch_regression_df[
  order(branch_regression_df$branch, -branch_regression_df$R2), ]

sig_branch <- branch_regression_df[branch_regression_df$R2 >= R2_THRESHOLD, ]
message(nrow(sig_branch), " branch-metabolite pairs with R² ≥ ", R2_THRESHOLD)

# ── Bubble plot: R² across branches per metabolite ─────────────────
# Keep only metabolites significant in at least one branch
sig_metabs <- unique(sig_branch$metabolite)
plot_df    <- branch_regression_df[branch_regression_df$metabolite %in% sig_metabs, ]
plot_df$sig_label <- ifelse(plot_df$padj < 0.05, "*", "")

bubble_plot <- ggplot(plot_df, aes(x = branch, y = metabolite, size = R2, color = R2)) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = sig_label), size = 5, vjust = 0.3, color = "black") +
  scale_size_continuous(range = c(1, 8), name = expression(R^2)) +
  scale_color_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                        midpoint = 0, name = "Slope") +
  labs(x = "Branch", y = "Metabolite",
       title = "Branch-level metabolite ~ pseudotime associations") +
  theme_classic() +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1),
        axis.text.y  = element_text(size = 8))

bubble_plot
# ── R² heatmap across branches ─────────────────────────────────────
library(pheatmap)

r2_wide <- reshape(
  branch_regression_df[branch_regression_df$metabolite %in% sig_metabs,
                        c("metabolite", "branch", "R2")],
  idvar     = "metabolite",
  timevar   = "branch",
  direction = "wide"
)
rownames(r2_wide) <- r2_wide$metabolite
r2_wide$metabolite <- NULL
colnames(r2_wide)  <- gsub("R2\\.", "", colnames(r2_wide))
r2_mat <- as.matrix(r2_wide)

heat_plot <- pheatmap(
  r2_mat,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = colorRampPalette(c("white", "#fee8c8", "#e34a33"))(50),
  display_numbers = round(r2_mat, 2),
  fontsize_number = 7,
  fontsize_row  = 8,
  main          = expression(paste("R"^2, " by branch")),
  silent        = TRUE
)
heat_plot


##### group proxy  + GAM

library(mgcv)
library(igraph)
library(pheatmap)
library(mgcv)
library(pheatmap)

R2_THRESHOLD <- 0.3

pt_vec     <- monocle3::pseudotime(cds)
finite_idx <- is.finite(pt_vec)
pt_all     <- pt_vec[finite_idx]
expr_all   <- expr_matrix[, finite_idx, drop = FALSE]
grp_all    <- colData(cds)$labels_refined[finite_idx]

# Manually defined branches by group condition
branches <- list(
    "EMCI1->EMCI2" = c("EMCI1", "EMCI2"),
  "EMCI1->LMCI1" = c("EMCI1", "LMCI1"),
  "EMCI1->LMCI2" = c("EMCI1", "LMCI2")
)

#branches <- list(
#  "EMCI1->EMCI2" = c("EMCI1", "EMCI2"),
#  "EMCI2->LMCI1" = c("EMCI2", "LMCI1"),
#  "LMCI1->LMCI2" = c("LMCI1", "LMCI2")
#)

run_branch_gam <- function(branch_name, groups) {
  idx <- grp_all %in% groups
  if (sum(idx) < 20) return(NULL)
  
  pt_scaled <- (pt_all[idx] - min(pt_all[idx])) / diff(range(pt_all[idx]))
  
  do.call(rbind, lapply(rownames(expr_all), function(metab) {
    fit <- tryCatch(mgcv::gam(y ~ s(pt, k = 5), method = "REML",
                              data = data.frame(y = expr_all[metab, idx], pt = pt_scaled)),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)
    data.frame(branch = branch_name, metabolite = metab,
               R2 = s$r.sq, dev_expl = s$dev.expl,
               pvalue = s$s.table[1, "p-value"], n_cells = sum(idx))
  }))
}

branch_regression_df <- do.call(rbind, Filter(Negate(is.null),
                                              mapply(run_branch_gam, names(branches), branches, SIMPLIFY = FALSE)
))
branch_regression_df$padj <- p.adjust(branch_regression_df$pvalue, method = "BH")

sig_branch <- branch_regression_df[branch_regression_df$R2 >= R2_THRESHOLD, ]
message(nrow(sig_branch), " branch-metabolite pairs with R² >= ", R2_THRESHOLD)

# ── Pheatmap ──────────────────────────────────────────────────────
sig_metabs <- unique(sig_branch$metabolite)

r2_wide <- reshape(
  branch_regression_df[branch_regression_df$metabolite %in% sig_metabs,
                       c("metabolite", "branch", "R2")],
  idvar = "metabolite", timevar = "branch", direction = "wide"
)
rownames(r2_wide)  <- r2_wide$metabolite
r2_wide$metabolite <- NULL
colnames(r2_wide)  <- gsub("R2\\.", "", colnames(r2_wide))
r2_mat             <- as.matrix(r2_wide)
r2_mat[is.na(r2_mat)] <- 0

p_heat <- pheatmap(
  r2_mat,
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,
  color           = colorRampPalette(c("white", "#fee8c8", "#e34a33"))(50),
  display_numbers = round(r2_mat, 2),
  fontsize_number = 7,
  fontsize_row    = 8,
  angle_col       = 45,
  main            = expression(paste("GAM  ", R^2, "  by branch")),
  silent          = TRUE
)
p_heat

# Save for Python plotting — rename to match assoc_df convention
assoc_df_export <- branch_regression_df
colnames(assoc_df_export)[colnames(assoc_df_export) == "metabolite"] <- "variable"
colnames(assoc_df_export)[colnames(assoc_df_export) == "branch"]     <- "trajectory"

write.csv(assoc_df_export, file.path(out_dir, "MONOCLE_assoc_df.csv"), row.names = FALSE)
# 

write.csv(r2_mat, file.path(out_dir, "MONOCLE_assoc_df_r2_0.5.csv"), row.names = FALSE)
# 
