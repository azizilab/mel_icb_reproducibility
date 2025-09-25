import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import scvi
import seaborn as sns
import scipy.stats
import pickle
import seaborn as sns

# custom decipher plotting function
def plot_decipher_v(adata, color, title="", show_axis="arrow", figsize=(3.5, 3.5), palette=None, subsample_frac=1.0, basis="decipher_v_corrected",
        x_label="Decipher 1", y_label="Decipher 2", ax_label_only_on_edges=False, ncols=2, **kwargs ):
    with plt.rc_context({"figure.figsize": figsize}):
        fig = sc.pl.embedding(
            sc.pp.subsample(adata, subsample_frac, copy=True), basis=basis, color=color, palette=palette, return_fig=True, frameon=(show_axis != "no"),
            ncols=ncols, **kwargs )
    ax = fig.axes[0]
    if type(color) == str or len(color) == 1:
        ax.set_title(title)

    for i, ax in enumerate(fig.axes[::2]):
        if show_axis == "arrow":
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.plot(1, 0, ">k", transform=ax.transAxes, clip_on=False)
            ax.plot(0, 1, "^k", transform=ax.transAxes, clip_on=False)

            if i % ncols == 0 or not ax_label_only_on_edges:
                ax.set_ylabel(y_label)
            else:
                ax.set_ylabel(None)
            if i // ncols == (len(color) - 1) // ncols or not ax_label_only_on_edges:
                ax.set_xlabel(x_label)
            else:
                ax.set_xlabel(None)
    return fig

def plot_trajectory(
    ax,
    trajectory,
    color="blue",
    o_size=7,
    star_size=20,
    linewidth=3
):
    ax.plot(
        trajectory['checkpoints'][:, 0],
        trajectory['checkpoints'][:, 1],
        marker="o",
        c="black",
        markerfacecolor=color,
        markersize=o_size,
        linewidth=linewidth,
    )
    ax.plot(
        trajectory['checkpoints'][:1, 0],
        trajectory['checkpoints'][:1, 1],
        marker="*",
        markersize=star_size,
        c="black",
        markerfacecolor=color,
    )

def calculate_mean_expression(ad, marker_genes, gene_set, layer = None):
    """
    Creates obs column with gene_set as name, containing mean expression and z scores
    for a given gene set for each cell
    """
    genes = marker_genes[gene_set] # set list of genes of interest
    print('genes not in adata from', gene_set, ':', list(set(genes).difference(set(ad.var_names))))
    genes = list(set(genes) & set(ad.var_names))# added 041924 due to genes not in HVG list
    
    # create dataframe to index
    if layer == None: # in case unintegrated logarithmized counts aren't saved as layer
        X = ad.X
        if not isinstance(X, np.ndarray): # added 081924 to check if already sparsified and made into np
            X = ad.X.todense()
        ad_df = pd.DataFrame(X)
    else:
        # ad_df = pd.DataFrame(ad.layers[layer]) # commented out 041924 due to errors with pd
        ad_df = ad.to_df(layer=layer)
    ad_df.columns = ad.var_names
    
    array = ad_df.loc[:,genes].to_numpy() # recover array of expression values for each gene of interest
    means = np.mean(array, axis = 1) # take row average
    z_scores = scipy.stats.zscore(means) # calculate z scores across each row
    
    if layer == None:
        key = gene_set + "_mean_expression"
    else:
        key = gene_set + "_mean_expression_" + layer
    ad.obs[key] = means # save mean expression value for each cell in obs of anndata
    
    if layer == None:
        key = gene_set + "_z_score"
    else:
        key = gene_set + "_z_score_" + layer
    ad.obs[key] = z_scores # save z scores for each cell in obs of anndata
    
def calculate_mean_expression_signed(ad, marker_genes, gene_set, layer = None):
    """
    Similar to above function, but calculate mean expression and z scores for a given
    gene set for each cell, assuming the gene_set has signed genes ("+" or "-" at end)
    """
    genes = marker_genes[gene_set] # set list of genes of interest
    
    # create dataframe to index
    if layer == None: # in case unintegrated logarithmized counts aren't saved as layer
        ad_df = pd.DataFrame(ad.X.todense())
    else:
        ad_df = pd.DataFrame(ad.layers[layer])
    ad_df.columns = ad.var_names
    
    # create list to track positive and negative genes, and strip +/- from gene name
    pos = []
    neg = []
    for gene in genes:
        if gene[-1] == '+':
            unsigned = gene.split('+')[0]
            pos.append(unsigned)
        else:
            unsigned = gene.split('-')[0]
            neg.append(unsigned)

    # recover array of expression values for each positive and negative gene of interest
    pos_array = ad_df.loc[:,pos].to_numpy() # keep the sign
    neg_array = -ad_df.loc[:,neg].to_numpy() # flip the sign
    array = np.concatenate((pos_array,neg_array), axis=1) # join columns (genes)
    
    means = np.mean(array, axis = 1) # take row average
    z_scores = scipy.stats.zscore(means) # calculate z scores across each row
    
    if layer == None:
        key = gene_set + "_mean_expression"
    else:
        key = gene_set + "_mean_expression_" + layer
    ad.obs[key] = means # save mean expression value for each cell in obs of anndata
    
    if layer == None:
        key = gene_set + "_z_score"
    else:
        key = gene_set + "_z_score_" + layer
    ad.obs[key] = z_scores # save z scores for each cell in obs of anndata