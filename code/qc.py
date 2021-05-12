import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

def plot_qc(data, sample_meta, gene_meta, to_plot=True):
    sample_meta['Library_size'] = data.sum(0)
    sample_meta['Genes_expressed'] = (data != 0).sum()
    sample_meta['Percent MT'] = data[pd.Series(data.index).apply(lambda x: x[:3] == 'MT-').values].sum() / data.sum()
    gene_meta['Total_count'] = data.sum(1)
    gene_meta['Expressed_samples'] = (data!= 0).sum(1)
    if to_plot:
        f, axs = plt.subplots(2, 2, figsize=(6, 6), dpi=300)

        axs[0, 0].hist(gene_meta['Expressed_samples'], bins=min([len(sample_meta), 100]), log=True)
        axs[0, 0].set_xlabel('Expressed in samples', size=10)
        axs[0, 0].set_ylabel('Number of genes', size=10)

        axs[0, 1].plot(sample_meta['Library_size'], sample_meta['Percent MT'], '.')
        axs[0, 1].set_xlabel('Library size', size=10)
        axs[0, 1].set_ylabel('Pecent MT', size=10)

        axs[1, 0].plot(sample_meta['Library_size'], sample_meta['Genes_expressed'], '.')
        axs[1, 0].set_xlabel('Library size', size=10)
        axs[1, 0].set_ylabel('Expressed genes', size=10)

        axs[1, 1].plot(np.log1p(data).mean(1), np.log1p(data).std(1), '.')
        axs[1, 1].set_xlabel('Mean expression', size=10)
        axs[1, 1].set_ylabel('Std expression', size=10)

        plt.tight_layout()
    return sample_meta, gene_meta
