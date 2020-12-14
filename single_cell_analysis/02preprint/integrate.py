import os
import sys
sys.path.insert(0, "../lib/")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import bbknn

import sankey
import sc_utils


SAMPLES = pd.read_csv("samples.csv")


def rename_genes(names):
    names = names.str.replace("^GRCh38_+", "")
    names = names.str.replace("^SARS-CoV-2i_", "SARS-CoV-2-")
    names = names.str.replace("SARS-CoV-2-antisense", "Antisense")
    return names


def load_ds(path):
    """
    H5 files are named like this GSM4698176_Sample_1_filtered_feature_bc_matrix.h5
    """
    fname = os.path.basename(path)
    sample = "_".join(fname.split("_")[1:3])
    ds = sc.read_10x_h5(path)
    ds.var_names = rename_genes(ds.var_names)
    ds.var_names_make_unique(join=".")
    ds.obs["orig.ident"] = sample
    ds.obs_names = sample + "_" + ds.obs_names.str.replace("-\d$", "")
    sc.pp.filter_cells(ds, min_genes=200)
    sc.pp.filter_genes(ds, min_cells=3)

    meta = SAMPLES.loc[SAMPLES.Sample == sample, :]
    ds.obs["Patient"] = meta.Patient.values[0]
    ds.obs["Day of intubation"] = meta["Day of intubation"].values[0]
    ds.obs["COVID-19"] = meta["COVID-19"].values[0]
    ds.obs["Sample"] = sample
    return ds


def integrate(
    ds_paths,
    h5ad_path,
    out_path,
    noribo,
    batch_hvg=False,
    hvg_genes=3000,
    hvg_remove=None,
    neighbors_total=50,
    cells_remove=None
):
    datasets = list(map(load_ds, ds_paths))
    ds = datasets[0].concatenate(datasets[1:], join="outer")
    if cells_remove:
        cells = pd.read_csv(cells_remove, header=None, index_col=0).iloc[:, 0]
        ds = ds[~ds.obs_names.isin(cells), :].copy()

    ds.var["mito"] = ds.var_names.str.startswith("MT-")
    ds.var["ribo"] = ds.var_names.str.match("^RP(L|S)")
    sc.pp.calculate_qc_metrics(
        ds,
        qc_vars=["mito", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    if noribo:
        ds = ds[:, ~ds.var["ribo"]]

    sc.pp.normalize_total(ds, target_sum=1e4)
    sc.pp.log1p(ds)
    if batch_hvg:
        sc.pp.highly_variable_genes(ds, n_top_genes=hvg_genes, batch_key="orig.ident")
    else:
        sc.pp.highly_variable_genes(ds, n_top_genes=2000)
    if hvg_remove is not None:
        ds.var.highly_variable[ds.var_names.str.match(hvg_remove)] = False
    ds.raw = ds

    sc.pp.scale(ds)
    sc.tl.pca(ds, svd_solver="arpack")
    n_neighbors = int(neighbors_total / len(datasets))
    bbknn.bbknn(ds, neighbors_within_batch=n_neighbors, n_pcs=30)
    sc.tl.leiden(ds, resolution=0.75)
    sc.tl.umap(ds)
    sc.tl.rank_genes_groups(ds, "leiden", method="wilcoxon", n_genes=0)

    ds.write(h5ad_path)

    os.makedirs(out_path)
    sc_utils.get_markers(ds, "leiden").to_csv(out_path + "/00markers.csv")

    mpl.rcParams["figure.figsize"] = (12, 8)
    ax = sc.pl.umap(ds, color="leiden", size=10, legend_loc="on data", show=False)
    ax.get_figure().savefig(out_path + "/01clusters.pdf")

    ax = sc.pl.umap(ds, color="orig.ident", size=10, show=False)
    ax.get_figure().savefig(out_path + "/02samples.pdf")

    bottom = np.zeros(len(ds.obs["leiden"].unique()))
    fig, ax = plt.subplots()
    for s in ds.obs["orig.ident"].unique():
        cnt = ds.obs["leiden"][ds.obs["orig.ident"] == s].value_counts().sort_index()
        ax.bar(cnt.index, cnt, bottom=bottom, label=s)
        bottom += cnt
    ax.legend()
    fig.suptitle("BBKNN Clusters by sample")
    fig.savefig(out_path + "/03composition.pdf")

    ax = sc.pl.stacked_violin(
        ds,
        ["MRC1", "FABP4", "CCL2", "SPP1", "CCL18", "CXCL10", "G0S2", "MKI67", "CD3E", "CD4", "CD8A", "JCHAIN", "FOXJ1"],
        groupby="leiden",
        rotation=90,
        figsize=(10, 10),
        show=False
    )
    ax[0].get_figure().savefig(out_path + "/04markers.pdf")

    ax = sc.pl.umap(ds, color="total_counts", size=10, show=False)
    ax.get_figure().savefig(out_path + "/05nUMI.pdf")

    ax = sc.pl.umap(ds, color="pct_counts_mito", size=10, show=False)
    ax.get_figure().savefig(out_path + "/06mito.pdf")

    ax = sc.pl.umap(ds, color="pct_counts_ribo", size=10, show=False)
    ax.get_figure().savefig(out_path + "/06ribo.pdf")
