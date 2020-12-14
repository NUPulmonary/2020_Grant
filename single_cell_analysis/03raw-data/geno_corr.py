import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os.path


def plot_corr(inputs, out):
    assert len(inputs) == 2

    sample1 = os.path.splitext(os.path.basename(inputs[0]))[0]
    sample2 = os.path.splitext(os.path.basename(inputs[1]))[0]

    t1 = pd.read_table(inputs[0], header=None, sep=r"\s+")
    t1.columns = ("Chr", "Pos", "Ref", "Alt", "C0", "C1")

    t2 = pd.read_table(inputs[1], header=None, sep=r"\s+")
    t2.columns = ("Chr", "Pos", "Ref", "Alt", "C0", "C1")

    t1 = t1[(t1.C0 != "./.") & (t1.C1 != "./.")]
    t2 = t2[(t2.C0 != "./.") & (t2.C1 != "./.")]

    t1['id'] = t1.Chr.astype(str) + "/" + t1.Pos.astype(str) + "/" + t1.Ref + "/" + t1.Alt
    t2['id'] = t2.Chr.astype(str) + "/" + t2.Pos.astype(str) + "/" + t2.Ref + "/" + t2.Alt

    common_ids = t1.id[t1.id.isin(t2.id)]

    common = pd.DataFrame(index=common_ids)
    t1.set_index('id', inplace=True)
    t2.set_index('id', inplace=True)

    common["{}-0".format(sample1)] = t1.C0[common.index]
    common["{}-1".format(sample1)] = t1.C1[common.index]
    common["{}-0".format(sample2)] = t2.C0[common.index]
    common["{}-1".format(sample2)] = t2.C1[common.index]
    common.replace({'0/0': 0, '0/1': 0.5, '1/1': 1}, inplace=True)
    corr = common.corr()

    plt.figure(figsize=(10, 8))
    plt.title("Correlation between genotype variants in {}/{}\n{} common variants".format(sample1, sample2, len(common_ids)))
    ax = sns.heatmap(
        corr,
        vmin=-1, vmax=1, center=0,
        cmap=sns.diverging_palette(20, 220, n=200),
        square=True,
        annot=True
    )
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=16, verticalalignment = 'center')
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=16)
    ax.title.set_size(18)
    plt.savefig(out)
