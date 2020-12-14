import scrublet as scr
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import gzip
import pandas as pd
import numpy as np


def run_scrublet(tenx_h5, doublet_rate=0.06, npca=40, save_to=None):
    if not save_to:
        raise ValueError("Please, specify prefix path where to save results to")
    if tenx_h5.endswith(".h5"):
        ds = sc.read_10x_h5(tenx_h5)
        counts_matrix = ds.X.tocsc().astype(np.longlong)
        obs = ds.obs.reset_index()
        obs.columns = ["0"]
    else:
        counts_matrix = scipy.io.mmread(gzip.open(tenx_h5 + '/matrix.mtx.gz')).T.tocsc()
        obs = pd.read_table(gzip.open(tenx_h5 + '/barcodes.tsv.gz'), header=None)
    #features = pd.read_table(gzip.open(input_dir + '/features.tsv.gz'), header=None)
    #genes = scr.make_genes_unique(features[1])
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=doublet_rate)
    doublet_scores, doublets = scrub.scrub_doublets(min_counts=2,
                                                    min_cells=3,
                                                    min_gene_variability_pctl=85,
                                                    n_prin_comps=npca)
    save_dir = os.path.dirname(save_to)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    obs['doublet'] = doublet_scores
    obs.to_csv(save_to + 'doublets.csv')
    scrub.plot_histogram()
    plt.savefig(save_to + 'doublet_hist.pdf')
    if not os.path.exists(save_to + 'threshold.txt'):
        with open(save_to + 'threshold.txt', 'w') as f:
            f.write(str(scrub.threshold_))
