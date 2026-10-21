"""Confounder-adjusted ADNI clinical trajectory for Reviewer 1 comments 3-4.

The cognitive variables are residualized against age, sex, and education using
models fitted only in cognitively normal (CN) participants.  CN means and
standard deviations are then used to convert residuals to normative z-scores.
Subtype is deliberately excluded from the ClinTrajan feature matrix and is
used only in downstream validation models.

Imported by 02_clinical_trajectory.py. Install the ClinTrajan helper modules
in clintra/ beside this file.
"""

from pathlib import Path
import random
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy.stats
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parent
SOURCE_ROOT = ROOT / "inputs"
sys.path.insert(0, str(ROOT / "clintra"))

from clintraj_eltree import (  # noqa: E402
    ExtendLeaves_modified,
    extract_trajectories,
    project_on_tree,
    prune_the_tree,
    quantify_pseudotime,
)
from clintraj_qi import (  # noqa: E402
    SVDcomplete_imputation_method,
    detect_variable_type,
    quantify_dataframe_univariate,
)
import elpigraph  # noqa: E402


INPUT = SOURCE_ROOT / "cn_emci_lmci_ad_clin_pseudotime.csv"
ADJUSTED_SCORES_OUTPUT = (
    ROOT / "extdata" / "cn_emci_lmci_ad_clin_pseudotime_covariate_adjusted.csv"
)
CN_MODELS_OUTPUT = ROOT / "extdata" / "cn_normative_cognitive_models.csv"
PSEUDOTIME_OUTPUT = ROOT / "extdata" / "pseudotime_covariate_adjusted.csv"
MODEL_OUTPUT = ROOT / "extdata" / "pseudotime_subtype_adjusted_models.csv"
TRAJECTORY_FIGURE = ROOT / "figs" / "Figure_5C_covariate_adjusted_trajectory.png"
SUBTYPE_FIGURE = ROOT / "figs" / "Figure_5D_pseudotime_by_subtype_adjusted.png"
HEATMAP_FIGURE = ROOT / "figs" / "Figure_5B_education_adjusted_zscore_heatmap.png"
GROUP_CURVE_FIGURE = ROOT / "figs" / "Figure_5C_education_adjusted_group_pseudotime.png"
GROUP_LOCATION_FIGURE = (
    ROOT / "figs" / "Figure_5D_education_adjusted_group_locations.png"
)
COMBINED_BCD_FIGURE = ROOT / "figs" / "Figure_5BCD_education_adjusted.png"

COGNITIVE_VARIABLES = [
    "CDRSB",
    "ADAS11",
    "ADAS13",
    "ADASQ4",
    "MMSE",
    "RAVLT_immediate",
    "RAVLT_learning",
    "RAVLT_forgetting",
    "RAVLT_perc_forgetting",
    "LDELTOTAL",
    "TRABSCOR",
    "FAQ",
]

GROUP_ORDER = [
    "CN",
    "emcisubtype1",
    "emcisubtype2",
    "lmcisubtype1",
    "lmcisubtype2",
    "AD",
]

GROUP_LABELS = {
    "CN": "CN",
    "emcisubtype1": "EMCI-1",
    "emcisubtype2": "EMCI-2",
    "lmcisubtype1": "LMCI-1",
    "lmcisubtype2": "LMCI-2",
    "AD": "AD",
}

GROUP_COLORS = {
    "CN": "#4E79A7",
    "emcisubtype1": "#F28E2B",
    "emcisubtype2": "#E15759",
    "lmcisubtype1": "#59A14F",
    "lmcisubtype2": "#B07AA1",
    "AD": "#9C755F",
}

# Final configuration selected by the user; inference is conditional on selection.
PRIMARY_TREE_PARAMETERS = {
    "number_of_nodes": 35, "alpha": 0.01, "mu": 0.20,
    "elastic_lambda": 0.010, "extension_control": 0.25,
}


def fit_cn_normative_scores(clinical, include_education=True):
    """Return CN-normed residual z-scores from the requested CN model."""
    required_covariates = ["AGE", "PTGENDER"]
    if include_education:
        required_covariates.append("PTEDUCAT")
    keep = clinical[required_covariates].notna().all(axis=1)
    if not keep.all():
        dropped = clinical.loc[~keep, "PID"].astype(str).tolist()
        print(
            "Dropping participants with missing adjustment covariates: "
            + ", ".join(dropped)
        )
    adjusted = clinical.loc[keep].reset_index(drop=True).copy()
    adjusted["PTGENDER"] = adjusted["PTGENDER"].astype(int)

    audit_rows = []
    for variable in COGNITIVE_VARIABLES:
        cn = adjusted.loc[
            adjusted["SNF"].eq("CN") & adjusted[variable].notna(),
            [variable] + required_covariates,
        ].copy()
        design_columns = [
            np.ones(len(cn)),
            cn["AGE"].to_numpy(float),
            cn["PTGENDER"].to_numpy(float),
        ]
        full_design_columns = [
            np.ones(len(adjusted)),
            adjusted["AGE"].to_numpy(float),
            adjusted["PTGENDER"].to_numpy(float),
        ]
        if include_education:
            design_columns.append(cn["PTEDUCAT"].to_numpy(float))
            full_design_columns.append(adjusted["PTEDUCAT"].to_numpy(float))
        design = np.column_stack(design_columns)
        outcome = cn[variable].to_numpy(float)
        coefficients = np.linalg.lstsq(design, outcome, rcond=None)[0]

        full_design = np.column_stack(full_design_columns)
        predicted = full_design @ coefficients
        residual = adjusted[variable].astype(float) - predicted
        cn_residual = residual.loc[adjusted["SNF"].eq("CN")]
        cn_mean = cn_residual.mean(skipna=True)
        cn_sd = cn_residual.std(skipna=True, ddof=1)
        if not np.isfinite(cn_sd) or cn_sd <= 0:
            raise ValueError(f"Invalid CN residual SD for {variable}: {cn_sd}")
        adjusted[f"Zadj_{variable}"] = (residual - cn_mean) / cn_sd

        audit_rows.append(
            {
                "variable": variable,
                "cn_n": len(cn),
                "intercept": coefficients[0],
                "beta_age": coefficients[1],
                "beta_male": coefficients[2],
                "beta_education": (
                    coefficients[3] if include_education else np.nan
                ),
                "cn_residual_mean": cn_mean,
                "cn_residual_sd": cn_sd,
            }
        )

    return adjusted, pd.DataFrame(audit_rows)


def fit_clintrajan(adjusted):
    """Fit the final fixed ClinTrajan configuration with the CN-side root."""
    parameters = PRIMARY_TREE_PARAMETERS.copy()
    feature_names = [f"Zadj_{variable}" for variable in COGNITIVE_VARIABLES]
    trajectory_input = adjusted[["PID"] + feature_names].copy()

    variable_types, _, _, _ = detect_variable_type(
        trajectory_input, 10, verbose=False
    )
    if set(variable_types) != {"CONTINUOUS"}:
        raise ValueError(
            "All adjusted cognitive variables should be continuous; found "
            f"{variable_types}"
        )

    quantified, _ = quantify_dataframe_univariate(
        trajectory_input, variable_types
    )
    imputed = SVDcomplete_imputation_method(
        quantified, variable_types, verbose=True, num_components=-1
    )
    x_raw = imputed[feature_names].to_numpy(dtype=float)
    x_scaled = scipy.stats.zscore(x_raw, axis=0, ddof=0)

    # Preserve the original 12-component setting. With 12 cognitive variables,
    # this is a rotation of the full adjusted feature space rather than a
    # dimension-reducing or label-informed transformation.
    reduced_dimension = 12
    pca = PCA(n_components=reduced_dimension, svd_solver="full")
    x = pca.fit_transform(x_scaled)

    print(f"ClinTrajan tree parameters: {parameters}")
    tree = elpigraph.computeElasticPrincipalTree(
        x,
        parameters["number_of_nodes"],
        alpha=parameters["alpha"],
        Mu=parameters["mu"],
        Lambda=parameters["elastic_lambda"],
        FinalEnergy="Penalized",
        Do_PCA=False,
    )[0]
    prune_the_tree(tree)
    tree_extended = ExtendLeaves_modified(
        x,
        tree,
        Mode="QuantDists",
        ControlPar=parameters["extension_control"],
        DoSA=False,
    )

    # Node numbers are reassigned whenever the tree is refitted, so a literal
    # node number is not a stable definition of the original root. Preserve
    # the original biological root (the CN/least-impaired end) by selecting
    # the terminal node closest to the CN centroid in the fitted space.
    edges = np.asarray(tree_extended["Edges"][0], dtype=int)
    degrees = np.bincount(edges.ravel(), minlength=len(tree_extended["NodePositions"]))
    leaf_nodes = np.flatnonzero(degrees == 1)
    cn_centroid = x[adjusted["SNF"].eq("CN").to_numpy()].mean(axis=0)
    leaf_positions = np.asarray(tree_extended["NodePositions"])[leaf_nodes]
    root_node = int(
        leaf_nodes[np.argmin(np.linalg.norm(leaf_positions - cn_centroid, axis=1))]
    )
    print(f"CN-side root node: {root_node}")
    trajectories, trajectory_edges = extract_trajectories(
        tree_extended, root_node
    )
    projection = project_on_tree(x, tree_extended)
    pseudotime_trajectories = quantify_pseudotime(
        trajectories, trajectory_edges, projection
    )
    return x, pca, tree_extended, root_node, pseudotime_trajectories


def collapse_pseudotime(adjusted, pseudotime_trajectories):
    """Collapse root-to-leaf records to one root-distance per participant."""
    records = []
    for trajectory_id, trajectory in enumerate(pseudotime_trajectories):
        terminal_node = trajectory["Trajectory"][-1]
        for point, pseudotime in zip(
            trajectory["Points"], trajectory["Pseudotime"]
        ):
            records.append(
                {
                    "row_index": int(point),
                    "trajectory_id": int(trajectory_id),
                    "terminal_node": int(terminal_node),
                    "pseudotime": float(pseudotime),
                }
            )
    long = pd.DataFrame(records)
    spread = long.groupby("row_index")["pseudotime"].agg(
        lambda value: value.max() - value.min()
    )
    if spread.max() > 1e-7:
        raise ValueError(
            "Shared root-to-leaf paths produced inconsistent pseudotime values"
        )

    collapsed = (
        long.groupby("row_index", as_index=False)
        .agg(
            pseudotime=("pseudotime", "mean"),
            trajectory_ids=(
                "trajectory_id",
                lambda value: ";".join(map(str, sorted(set(value)))),
            ),
            terminal_nodes=(
                "terminal_node",
                lambda value: ";".join(map(str, sorted(set(value)))),
            ),
        )
        .sort_values("row_index")
    )
    metadata = adjusted[
        ["PID", "SNF", "AGE", "PTGENDER", "PTEDUCAT", "SNFSUB"]
    ].copy()
    result = metadata.merge(collapsed, on="row_index", how="left") \
        if "row_index" in metadata.columns else metadata.join(
            collapsed.set_index("row_index"), how="left"
        )
    if result["pseudotime"].isna().any():
        missing = result.loc[result["pseudotime"].isna(), "PID"].tolist()
        raise ValueError(f"Participants without pseudotime: {missing}")
    return result


def fit_subtype_models(pseudotime, include_education=True):
    """Fit stage-specific subtype models with the requested covariates."""
    specifications = [
        ("EMCI", "emcisubtype1", "emcisubtype2"),
        ("LMCI", "lmcisubtype1", "lmcisubtype2"),
    ]
    model_rows = []
    for stage, subtype1, subtype2 in specifications:
        data = pseudotime.loc[
            pseudotime["SNF"].isin([subtype1, subtype2])
        ].copy()
        data["subtype2"] = data["SNF"].eq(subtype2).astype(int)
        data["PTGENDER"] = data["PTGENDER"].astype("category")

        covariate_terms = "AGE + C(PTGENDER)"
        if include_education:
            covariate_terms += " + PTEDUCAT"
        full_formula = f"pseudotime ~ subtype2 + {covariate_terms}"
        reduced_formula = f"pseudotime ~ {covariate_terms}"
        full = smf.ols(full_formula, data=data).fit(cov_type="HC3")
        reduced = smf.ols(reduced_formula, data=data).fit()
        ordinary_full = smf.ols(full_formula, data=data).fit()
        partial_r2 = max(
            0.0,
            (reduced.ssr - ordinary_full.ssr) / reduced.ssr,
        )
        ci_low, ci_high = full.conf_int().loc["subtype2"]
        model_rows.append(
            {
                "stage": stage,
                "contrast": f"{GROUP_LABELS[subtype2]} vs {GROUP_LABELS[subtype1]}",
                "n": int(full.nobs),
                "beta_subtype2": full.params["subtype2"],
                "hc3_se": full.bse["subtype2"],
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "p_value": full.pvalues["subtype2"],
                "partial_r2": partial_r2,
                "model_r2": ordinary_full.rsquared,
                "covariates": (
                    "age + sex + education"
                    if include_education
                    else "age + sex"
                ),
            }
        )
    models = pd.DataFrame(model_rows)
    # Standard monotone Benjamini-Hochberg adjustment across the two planned
    # stage-specific tests.
    order = np.argsort(models["p_value"].to_numpy())
    ranked = models["p_value"].to_numpy()[order] * len(models) / (
        np.arange(len(models)) + 1
    )
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted_p = np.empty(len(models))
    adjusted_p[order] = np.minimum(ranked, 1.0)
    models["p_value_bh"] = adjusted_p
    return models






def adjusted_heatmap_clustering(adjusted):
    """Return group means, their clustered order, and both linkages."""
    feature_names = [f"Zadj_{variable}" for variable in COGNITIVE_VARIABLES]
    # Covariate coefficients were estimated in CN. For a readable heatmap, the
    # resulting residuals are re-standardized across the analysis cohort. This
    # display scaling does not alter the trajectory fit or inferential models.
    display_z = adjusted[feature_names].apply(
        lambda value: (value - value.mean(skipna=True)) / value.std(skipna=True)
    )
    display_z["SNF"] = adjusted["SNF"].to_numpy()
    means = (
        display_z.groupby("SNF", observed=True)[feature_names]
        .mean()
        .reindex(GROUP_ORDER)
    )
    means.index = [GROUP_LABELS[group] for group in means.index]
    means.columns = [name.removeprefix("Zadj_") for name in means.columns]

    row_distance = pdist(means.to_numpy(), metric="euclidean")
    column_distance = pdist(means.to_numpy().T, metric="euclidean")
    row_linkage = linkage(row_distance, method="average")
    row_linkage = rotate_linkage_to_order(
        row_linkage, desired_leaf_order=list(range(len(GROUP_ORDER)))
    )
    column_linkage = linkage(column_distance, method="average")
    row_order = leaves_list(row_linkage)
    column_order = leaves_list(column_linkage)
    ordered = means.iloc[row_order, column_order]
    return means, ordered, row_linkage, column_linkage


def rotate_linkage_to_order(linkage_matrix, desired_leaf_order):
    """Rotate dendrogram branches without changing clustering topology.

    This is the SciPy equivalent of dendextend::rotate used in
    5_ADNI_clinical_scores.ipynb. It only swaps the left and right children of
    internal nodes, preserving merge heights, memberships, and distances.
    """
    rotated = np.asarray(linkage_matrix, dtype=float).copy()
    number_of_leaves = rotated.shape[0] + 1
    desired_rank = {
        int(leaf): rank for rank, leaf in enumerate(desired_leaf_order)
    }

    def orient(node):
        node = int(node)
        if node < number_of_leaves:
            return [node]
        row = node - number_of_leaves
        left = int(rotated[row, 0])
        right = int(rotated[row, 1])
        left_leaves = orient(left)
        right_leaves = orient(right)
        if min(desired_rank[leaf] for leaf in left_leaves) > min(
            desired_rank[leaf] for leaf in right_leaves
        ):
            rotated[row, 0], rotated[row, 1] = rotated[row, 1], rotated[row, 0]
            left_leaves, right_leaves = right_leaves, left_leaves
        return left_leaves + right_leaves

    orient(2 * number_of_leaves - 2)
    observed_order = leaves_list(rotated).tolist()
    if observed_order != list(desired_leaf_order):
        raise ValueError(
            "The requested row order is incompatible with the fitted "
            f"dendrogram topology: observed {observed_order}"
        )
    return rotated




def draw_clustered_heatmap_panel(
    figure,
    subplot_spec,
    heatmap_data,
    row_linkage,
    column_linkage,
    adjustment_label,
):
    """Draw a heatmap with visible row and column dendrograms."""
    grid = subplot_spec.subgridspec(
        2,
        3,
        height_ratios=[0.20, 0.80],
        width_ratios=[0.17, 0.78, 0.05],
        hspace=0.02,
        wspace=0.04,
    )
    corner_ax = figure.add_subplot(grid[0, 0])
    column_ax = figure.add_subplot(grid[0, 1])
    row_ax = figure.add_subplot(grid[1, 0])
    heatmap_ax = figure.add_subplot(grid[1, 1])
    colorbar_ax = figure.add_subplot(grid[1, 2])
    corner_ax.axis("off")

    dendrogram(
        column_linkage,
        ax=column_ax,
        orientation="top",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#555555",
        link_color_func=lambda _: "#555555",
    )
    dendrogram(
        row_linkage,
        ax=row_ax,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#555555",
        link_color_func=lambda _: "#555555",
    )
    row_ax.invert_yaxis()
    for dendrogram_ax in [column_ax, row_ax]:
        dendrogram_ax.set_xticks([])
        dendrogram_ax.set_yticks([])
        for spine in dendrogram_ax.spines.values():
            spine.set_visible(False)

    limit = max(1.0, float(np.nanmax(np.abs(heatmap_data.to_numpy()))))
    sns.heatmap(
        heatmap_data,
        cmap="plasma",
        center=0,
        vmin=-limit,
        vmax=limit,
        linewidths=0,
        cbar=True,
        cbar_ax=colorbar_ax,
        cbar_kws={"label": "Mean adjusted z-score"},
        ax=heatmap_ax,
    )
    heatmap_ax.set_xlabel("")
    heatmap_ax.set_ylabel("Group")
    heatmap_ax.tick_params(axis="x", rotation=45, labelsize=8)
    heatmap_ax.tick_params(axis="y", rotation=0)
    colorbar_ax.tick_params(labelsize=7)
    colorbar_ax.yaxis.label.set_size(8)
    column_ax.set_title(
        f"{adjustment_label}-adjusted cognitive z-scores",
        fontsize=12,
        pad=5,
    )
    return heatmap_ax


def select_cn_to_ad_trajectory(x, adjusted, tree, pseudotime_trajectories):
    """Select the root-to-leaf trajectory whose endpoint is nearest the AD centroid."""
    ad_centroid = x[adjusted["SNF"].eq("AD").to_numpy()].mean(axis=0)
    node_positions = np.asarray(tree["NodePositions"])
    terminal_nodes = np.asarray(
        [trajectory["Trajectory"][-1] for trajectory in pseudotime_trajectories],
        dtype=int,
    )
    selected_index = int(
        np.argmin(
            np.linalg.norm(node_positions[terminal_nodes] - ad_centroid, axis=1)
        )
    )
    return selected_index, pseudotime_trajectories[selected_index]


def group_curve_data(x, adjusted, tree, pseudotime_trajectories):
    selected_index, trajectory = select_cn_to_ad_trajectory(
        x, adjusted, tree, pseudotime_trajectories
    )
    points = np.asarray(trajectory["Points"], dtype=int)
    pseudotime = np.asarray(trajectory["Pseudotime"], dtype=float)
    group_code = adjusted.loc[points, "SNF"].map(
        {group: index for index, group in enumerate(GROUP_ORDER)}
    ).to_numpy(float)

    order = np.argsort(pseudotime)
    x_sorted = pseudotime[order]
    y_sorted = group_code[order]
    smooth = lowess(y_sorted, x_sorted, frac=0.30, it=1, return_sorted=True)
    fitted = np.interp(x_sorted, smooth[:, 0], smooth[:, 1])
    total_sum_squares = np.sum((y_sorted - y_sorted.mean()) ** 2)
    descriptive_r2 = 1 - np.sum((y_sorted - fitted) ** 2) / total_sum_squares

    number_of_bins = 14
    breaks = np.linspace(x_sorted.min(), x_sorted.max(), number_of_bins + 1)
    bin_id = np.clip(np.digitize(x_sorted, breaks) - 1, 0, number_of_bins - 1)
    sliding_x = []
    sliding_y = []
    for index in range(number_of_bins):
        in_bin = bin_id == index
        if in_bin.any():
            sliding_x.append(x_sorted[in_bin].mean())
            sliding_y.append(y_sorted[in_bin].mean())

    return {
        "selected_index": selected_index,
        "root_node": int(trajectory["Trajectory"][0]),
        "terminal_node": int(trajectory["Trajectory"][-1]),
        "points": points,
        "pseudotime": pseudotime,
        "group_code": group_code,
        "smooth": smooth,
        "sliding_x": np.asarray(sliding_x),
        "sliding_y": np.asarray(sliding_y),
        "descriptive_r2": descriptive_r2,
    }


def draw_group_curve(ax, curve, models, test_label="age, sex, and education"):
    """Draw panel C using the visual style of the original figure."""
    ax.set_facecolor("#E5E5E5")
    ax.set_axisbelow(True)
    ax.grid(True, color="white", linewidth=1.0)
    ax.scatter(
        curve["pseudotime"],
        curve["group_code"],
        color="red",
        s=24,
        alpha=0.90,
        edgecolor="none",
        label="data",
    )
    ax.plot(
        curve["smooth"][:, 0],
        curve["smooth"][:, 1],
        color="#1B7837",
        linewidth=2.5,
        label="predicted",
    )
    ax.plot(
        curve["sliding_x"],
        curve["sliding_y"],
        color="#2166AC",
        linewidth=2.2,
        label="sliding av",
    )
    ax.set_yticks(range(len(GROUP_ORDER)))
    ax.set_yticklabels(["CN", "EMCI 1", "EMCI 2", "LMCI 1", "LMCI 2", "AD"])
    ax.set_ylim(-0.5, len(GROUP_ORDER) - 0.5)
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Conditions")
    ax.set_title(
        f"Trajectory: {curve['root_node']}->{curve['terminal_node']}, "
        f"r2={curve['descriptive_r2']:.2f}"
    )
    ax.legend(
        loc="lower right",
        fontsize=8,
        frameon=True,
        facecolor="#E5E5E5",
        edgecolor="#D0D0D0",
    )
    for spine in ax.spines.values():
        spine.set_color("white")


def clintrajan_2d_layout(x, tree, scatter_parameter=0.03):
    """Reproduce ClinTrajan's 2D Kamada-Kawai display projection.

    The elastic tree is fitted in the full cognitive space. ClinTrajan's
    visualization first lays out the graph in two dimensions, projects each
    participant to an edge, and adds a small orthogonal displacement based on
    distance from the fitted tree. This is display-only and does not alter
    pseudotime.
    """
    edges = np.asarray(tree["Edges"][0], dtype=int)
    projection = project_on_tree(x, tree)
    projection_values = np.asarray(projection["ProjectionValues"], dtype=float)
    edge_ids = np.asarray(projection["EdgeID"], dtype=int)
    x_projected = np.asarray(projection["X_projected"], dtype=float)

    distance_to_projection = np.sum((x - x_projected) ** 2, axis=1)
    distance_to_projection -= np.percentile(distance_to_projection, 20)

    graph = nx.Graph()
    graph.add_nodes_from(range(len(tree["NodePositions"])))
    graph.add_edges_from(edges)
    positions = nx.kamada_kawai_layout(graph, scale=2)
    node_xy = np.asarray([positions[index] for index in range(len(positions))])

    rng = np.random.default_rng(0)
    point_xy = np.zeros((len(x), 2), dtype=float)
    for index in range(len(x)):
        start, end = edges[edge_ids[index]]
        fraction = projection_values[index]
        clipped_fraction = np.clip(fraction, 0.0, 1.0)
        base = node_xy[start] + (node_xy[end] - node_xy[start]) * clipped_fraction

        edge_vector = node_xy[end] - node_xy[start]
        if fraction < 0 or fraction > 1:
            edge_vector = rng.random(2) - 0.5
        edge_norm = np.linalg.norm(edge_vector)
        if edge_norm == 0:
            edge_vector = np.array([1.0, 0.0])
            edge_norm = 1.0
        edge_unit = edge_vector / edge_norm
        orthogonal = np.array([edge_unit[1], -edge_unit[0]])
        radius = (
            np.sqrt(max(distance_to_projection[index], 0.0))
            * scatter_parameter
        )
        side = 1 if rng.random() < 0.5 else -1
        point_xy[index] = base + orthogonal * radius * side
    return point_xy, node_xy, edges


def draw_tree(
    ax,
    point_xy,
    node_xy,
    edges,
    point_colors,
    root_node=None,
    title=None,
    highlight_mask=None,
    show_node_numbers=False,
):
    ax.set_facecolor("#E5E5E5")
    ax.set_axisbelow(True)
    ax.grid(True, color="white", linewidth=0.7)
    if highlight_mask is None:
        ax.scatter(
            point_xy[:, 0],
            point_xy[:, 1],
            c=point_colors,
            s=6,
            alpha=0.72,
            linewidth=0,
            rasterized=True,
            zorder=2,
        )
    else:
        ax.scatter(
            point_xy[~highlight_mask, 0],
            point_xy[~highlight_mask, 1],
            color="#5A5A5A",
            s=4,
            alpha=0.48,
            linewidth=0,
            rasterized=True,
            zorder=1,
        )
        ax.scatter(
            point_xy[highlight_mask, 0],
            point_xy[highlight_mask, 1],
            color="#D73027",
            s=8,
            alpha=0.62,
            linewidth=0,
            rasterized=True,
            zorder=2,
        )
    for start, end in edges:
        ax.plot(
            node_xy[[start, end], 0],
            node_xy[[start, end], 1],
            color="#404040",
            linewidth=1.0,
            zorder=3,
        )
    ax.scatter(
        node_xy[:, 0],
        node_xy[:, 1],
        s=10,
        color="#202020",
        zorder=4,
    )
    if root_node is not None:
        ax.scatter(
            node_xy[root_node, 0],
            node_xy[root_node, 1],
            s=80,
            marker="*",
            color="#FFD700",
            edgecolor="black",
            linewidth=0.6,
            zorder=5,
        )
        ax.annotate(
            "Root Node",
            xy=node_xy[root_node],
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=5,
            color="#303030",
            zorder=7,
        )
    if show_node_numbers:
        for node, (x_coordinate, y_coordinate) in enumerate(node_xy):
            ax.text(
                x_coordinate,
                y_coordinate,
                str(node),
                fontsize=4.5,
                ha="center",
                va="center",
                color="white",
                zorder=6,
            )
    if title is not None:
        ax.set_xlabel(title, fontsize=11, labelpad=5)
    ax.tick_params(axis="both", labelsize=4, length=0, colors="#666666")
    ax.set_box_aspect(1)
    for spine in ax.spines.values():
        spine.set_color("white")


def draw_group_location_panels(
    axes, point_xy, node_xy, edges, root_node, adjusted
):
    overall_colors = np.repeat("#D73027", len(adjusted))
    draw_tree(
        axes[0],
        point_xy,
        node_xy,
        edges,
        overall_colors,
        root_node=root_node,
        title="Overall",
        show_node_numbers=True,
    )
    for ax, group in zip(axes[1:], GROUP_ORDER):
        highlight_mask = adjusted["SNF"].eq(group).to_numpy()
        draw_tree(
            ax,
            point_xy,
            node_xy,
            edges,
            "#303030",
            title=GROUP_LABELS[group],
            highlight_mask=highlight_mask,
        )

    all_x = np.concatenate([node_xy[:, 0], point_xy[:, 0]])
    all_y = np.concatenate([node_xy[:, 1], point_xy[:, 1]])
    x_padding = max(0.08, 0.04 * np.ptp(all_x))
    y_padding = max(0.08, 0.04 * np.ptp(all_y))
    for ax in axes:
        ax.set_xlim(all_x.min() - x_padding, all_x.max() + x_padding)
        ax.set_ylim(all_y.min() - y_padding, all_y.max() + y_padding)


def plot_requested_panels(
    x,
    tree,
    root_node,
    adjusted,
    pseudotime_trajectories,
    models,
    adjustment_label="Age-, sex-, and education",
    file_tag="education_adjusted",
    figure_title=None,
):
    (
        heatmap_means,
        heatmap_data,
        row_linkage,
        column_linkage,
    ) = adjusted_heatmap_clustering(adjusted)
    curve = group_curve_data(x, adjusted, tree, pseudotime_trajectories)
    point_xy, node_xy, edges = clintrajan_2d_layout(x, tree)
    test_label = (
        "age, sex, and education"
        if "education" in adjustment_label.lower()
        else "age and sex"
    )
    heatmap_figure = ROOT / "figs" / f"Figure_5B_{file_tag}_zscore_heatmap.png"
    group_curve_figure = (
        ROOT / "figs" / f"Figure_5C_{file_tag}_group_pseudotime.png"
    )
    group_location_figure = (
        ROOT / "figs" / f"Figure_5D_{file_tag}_group_locations.png"
    )
    combined_bcd_figure = ROOT / "figs" / f"Figure_5BCD_{file_tag}.png"

    limit = max(1.0, float(np.nanmax(np.abs(heatmap_data.to_numpy()))))
    cluster = sns.clustermap(
        heatmap_means,
        row_linkage=row_linkage,
        col_linkage=column_linkage,
        cmap="plasma",
        center=0,
        vmin=-limit,
        vmax=limit,
        linewidths=0,
        figsize=(9.5, 6.3),
        cbar_kws={"label": "Mean adjusted z-score"},
        dendrogram_ratio=(0.13, 0.18),
        cbar_pos=(0.95, 0.36, 0.018, 0.30),
    )
    cluster.ax_heatmap.set_xlabel("")
    cluster.ax_heatmap.set_ylabel("Group")
    cluster.ax_heatmap.tick_params(axis="x", rotation=45, labelsize=8)
    cluster.ax_heatmap.tick_params(
        axis="y", rotation=0, labelleft=True, labelright=False
    )
    cluster.ax_heatmap.yaxis.tick_left()
    cluster.ax_heatmap.yaxis.set_label_position("left")
    cluster.fig.suptitle(
        f"{adjustment_label}-adjusted cognitive z-scores",
        fontsize=15,
        y=1.02,
    )
    cluster.fig.savefig(heatmap_figure, dpi=600, bbox_inches="tight")
    cluster.fig.savefig(heatmap_figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(cluster.fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    draw_group_curve(ax, curve, models, test_label=test_label)
    fig.tight_layout()
    fig.savefig(group_curve_figure, dpi=600, bbox_inches="tight")
    fig.savefig(group_curve_figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 7, figsize=(18, 4.0))
    draw_group_location_panels(
        axes, point_xy, node_xy, edges, root_node, adjusted
    )
    fig.suptitle(
        f"Location of clinical groups on the {adjustment_label.lower()}-adjusted trajectory",
        fontsize=15,
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(group_location_figure, dpi=600, bbox_inches="tight")
    fig.savefig(group_location_figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    combined = plt.figure(figsize=(18, 10.5))
    outer = combined.add_gridspec(
        2,
        2,
        height_ratios=[1.15, 0.85],
        width_ratios=[1.05, 1.0],
        hspace=0.34,
        wspace=0.28,
    )
    curve_ax = combined.add_subplot(outer[0, 1])
    lower = outer[1, :].subgridspec(1, 7, wspace=0.13)
    tree_axes = [combined.add_subplot(lower[0, index]) for index in range(7)]

    heatmap_ax = draw_clustered_heatmap_panel(
        combined,
        outer[0, 0],
        heatmap_data,
        row_linkage,
        column_linkage,
        adjustment_label,
    )
    draw_group_curve(curve_ax, curve, models, test_label=test_label)
    draw_group_location_panels(
        tree_axes, point_xy, node_xy, edges, root_node, adjusted
    )
    combined.text(0.012, 0.865, "B", fontsize=25, fontweight="bold")
    combined.text(0.505, 0.865, "C", fontsize=25, fontweight="bold")
    combined.text(0.012, 0.405, "D", fontsize=25, fontweight="bold")
    if figure_title is None:
        figure_title = f"ADNI cognitive trajectory adjusted for {test_label}"
    combined.suptitle(figure_title, fontsize=18, y=0.995)
    combined.savefig(combined_bcd_figure, dpi=600, bbox_inches="tight")
    combined.savefig(combined_bcd_figure.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(combined)
    return curve, {
        "heatmap": heatmap_figure,
        "group_curve": group_curve_figure,
        "group_locations": group_location_figure,
        "combined": combined_bcd_figure,
    }


