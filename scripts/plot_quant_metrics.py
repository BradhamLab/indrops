import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def plot_distributions(series, unit, norm=True):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=False,
                                   figsize=(13, 8))
    # plot pdf 
    sns.distplot(series, kde=False, ax=ax1)
    if norm:
        locs = ax1.get_yticks()
        new_locs = np.array(locs) / metrics.shape[0] * 100
        new_labels = ["{:0.0f}".format(x) for x in new_locs]
        ax1.set_yticklabels(new_labels)
        ax1.set_ylabel("Percent {}".format(unit))
    # plot cdf
    sns.distplot(series, hist_kws={'cumulative': True}, kde=False,
                 norm_hist=True, ax=ax2)
    ax2.set_ylabel("Cumulative Percent {}".format(unit))
    plt.tight_layout()

if __name__ == '__main__':
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        metrics = pd.read_csv(snakemake.input['metrics'],
                              index_col=0,
                              delimiter='\t')
        metrics['Percent Aligned'] = metrics.apply(lambda x:\
                    x["Reads with unique alignment"] / x['Reads'] * 100, axis=1)
        plot_distributions(metrics['Percent Aligned'], 'Barcodes')
        plt.savefig(snakemake.output['aligned'])
        plt.close()
        plt.clf()
        plt.cla()
        plot_distributions(metrics['UMIFM'], 'Barcodes')
        plt.savefig(snakemake.output['umi'])
        plt.close()
        plt.clf()
        plt.cla()