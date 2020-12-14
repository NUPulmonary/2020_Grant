import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def expr_colormap():
    cdict = {
        'red':   [
            (0.0, 220/256, 220/256),
            (0.5, 42/256, 42/256),
            (1.0, 6/256, 6/256)
        ],

        'green': [
            (0.0, 220/256, 220/256),
            (0.5, 145/256, 145/256),
            (1.0, 37/256, 27/256)
        ],

        'blue':  [
            (0.0, 220/256, 220/256),
            (0.5, 174/256, 174/256),
            (1.0, 170/256, 170/256)
        ]
    }
    return mpl.colors.LinearSegmentedColormap('testCmap', segmentdata=cdict, N=256)


def get_markers(
    anndata,
    groupby,
    key="rank_genes_groups",
    p_val_cutoff=0.05,
    logfc_cutoff=0.5
):
    markers = pd.concat([
        pd.DataFrame(anndata.uns[key]["names"]).melt(),
        pd.DataFrame(anndata.uns[key]["pvals_adj"]).melt(),
        pd.DataFrame(anndata.uns[key]["logfoldchanges"]).melt()
    ], axis=1)
    markers.columns = ("cluster", "gene", "cluster2", "p_val_adj", "cluster3", "avg_logFC")
    markers = markers.loc[:, ["cluster", "gene", "avg_logFC", "p_val_adj"]]
    markers = markers.loc[markers.avg_logFC > logfc_cutoff, ]
    markers = markers.loc[markers.p_val_adj < p_val_cutoff, ]

    for cluster in markers.cluster.unique():
        cells = anndata.obs[groupby] == cluster
        in_cluster_selector = markers.cluster == cluster
        genes = markers.gene[in_cluster_selector]
        in_cluster = np.sum(anndata.raw[cells, genes].X > 0, axis=0) / cells.sum()
        markers.loc[in_cluster_selector, "pct.1"] = in_cluster.T
        other_cells = anndata.obs[groupby] != cluster
        other_clusters = np.sum(anndata.raw[other_cells, genes].X > 0, axis=0) / other_cells.sum()
        markers.loc[in_cluster_selector, "pct.2"] = other_clusters.T

    # all_clusters = markers.cluster.unique()
    # if all_clusters.size == 1:
    #     all_clusters = list(anndata.obs[groupby].astype(str).unique())
    #     if "nan" in all_clusters:
    #         all_clusters.remove("nan")
    #     all_clusters = pd.Series(all_clusters)

    # markers["pct.2"] = markers.apply(calc_pct_2, axis=1, args=(all_clusters,))
    markers["p_val"] = markers.p_val_adj
    markers = markers.loc[:, ["p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]]
    return markers


def find_markers(obj, groupby, groups):
    if len(groups) != 2:
        raise ValueError("Expecting exactly 2 elements in `groups`")
    key_pos = f"degs_{groups[0]}_{groups[1]}_pos"
    key_neg = f"degs_{groups[0]}_{groups[1]}_neg"
    sc.tl.rank_genes_groups(
        obj,
        groupby,
        groups=[groups[0]],
        reference=groups[1],
        key_added=key_pos,
        n_genes=0
    )
    sc.tl.rank_genes_groups(
        obj,
        groupby,
        groups=[groups[1]],
        reference=groups[0],
        key_added=key_neg,
        n_genes=0
    )
    pos = sc_utils.get_markers(dsc, "CD4_dbl", key="cd4_dbl_pos")
    neg = sc_utils.get_markers(dsc, "CD4_dbl", key="cd4_dbl_neg")


def feature_plot(ds, feature, gridsize=(180, 70), linewidths=0.15, figsize=None):
    if feature in ds.obs.columns:
        values = ds.obs_vector(feature)
    else:
        values = ds.raw.obs_vector(feature)

    kwargs = {}
    if figsize is not None:
        kwargs["figsize"] = figsize
    fig, ax = plt.subplots(**kwargs)
    hb = ax.hexbin(
        ds.obsm["X_umap"][:, 0],
        ds.obsm["X_umap"][:, 1],
        C=values,
        # cmap="YlOrRd",
        cmap=expr_colormap(),
        gridsize=gridsize,
        linewidths=linewidths
    )
    cax = fig.add_axes((0.92, 0.8, 0.02, 0.15))
    cb = fig.colorbar(hb, cax=cax, fraction=0.05, pad=0.02, aspect=40)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{feature}")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_composition(ds, group_by, color, relative=False, palette=None):
    left = np.zeros(len(ds.obs[group_by].unique()))
    total = None
    if relative:
        total = ds.obs[group_by].value_counts().sort_index()
    fig, ax = plt.subplots()
    num_colors = ds.obs[color].unique().size
    # TODO: adjust
    if palette is not None:
        colors = palette
    elif num_colors <= 10:
        colors = mpl.cm.tab10
    elif num_colors <= 20:
        colors = mpl.cm.tab20
    else:
        colors = sc.pl.palettes.default_28
    for i, s in enumerate(ds.obs[color].cat.categories):
        cnt = ds.obs[group_by][ds.obs[color] == s].value_counts().sort_index()
        if relative:
            cnt = cnt / total * 100
        c = isinstance(colors, list) and colors[i] or colors(i)
        ax.barh(cnt.index, cnt, left=left, label=s, color=c)
        left += cnt
    ax.legend(title=color.capitalize())
    ax.set_title(f"{group_by.capitalize()} by {color}")
    return ax
