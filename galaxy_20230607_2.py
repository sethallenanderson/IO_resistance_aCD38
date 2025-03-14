#!/usr/bin/env python3.6

'''
Galaxy plot
v1.0

Purpose: creates Kevin's galaxy plot
Input: h5ad scanpy file with experiment labels in adata.obs and embedding already performed
Output: galaxy plots

Version history:
v1.0: implemented

Last modified by Peter Du on Mar 2019

Modified by Celeste Nobrega April 2024
'''


import sys, os, argparse
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def getParser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    required = parser.add_argument_group('required')
    required.add_argument('--input',
                          type=str,
                          help='H5AD scanpy file with experiment annotations in .obs',
                          required=True)
    parser.add_argument('--colname',
                        type=str,
                        help='Name of experiment label column in .obs; default=batch',
                        default='batch')
    parser.add_argument('--embed',
                        type=str,
                        help='embedding method; default=umap',
                        choices=['umap', 'tsne', 'fa'],
                        default='umap')
    return parser


def galaxy_plot(adata, outdir, colname='batch', embed='X_umap'):
    keys = adata.obs[colname].unique()
    fig, ax = plt.subplots(1, len(keys) + 1, figsize=(8 * (len(keys) + 1), 8), sharex=True, sharey=True)
    cmap = plt.cm.get_cmap('magma')
    n_levels = 50

    global_bounds={}

    df = pd.DataFrame({'x': adata.obsm[embed][:, 0],
                       'y': adata.obsm[embed][:, 1]})
    print(df.x)
    print(df.y)
    g = sns.kdeplot(x=df.x, y=df.y, n_levels=n_levels, cmap=cmap, ax=ax[0], shade=True)
    g.scatter(x=df.x, y=df.y, s=8, c='w', marker='.', linewidths=0, edgecolors='none')
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].axis('on')
    ax[0].set_title('all')
    ax[0].set_facecolor(cmap(1.0 / (n_levels * 2)))
    global_bounds['x'] = (ax[0].get_xlim())
    global_bounds['y'] = (ax[0].get_ylim())
    for i, key in enumerate(keys):
        df = pd.DataFrame({'x': adata[adata.obs[colname] == key].obsm[embed][:, 0],
                           'y': adata[adata.obs[colname] == key].obsm[embed][:, 1]})
        g = sns.kdeplot(x=df.x, y=df.y, n_levels=n_levels, cmap=cmap, ax=ax[i + 1], shade=True)
        g.scatter(x=df.x, y=df.y, s=8, c='w', marker='.', linewidths=0, edgecolors='none')
        ax[i + 1].set_xticks([])
        ax[i + 1].set_yticks([])
        ax[i + 1].axis('on')
        ax[i + 1].set_title(key)
        ax[i + 1].set_facecolor(cmap(1.0 / (n_levels * 2)))

    plt.tight_layout()
    outpath = os.path.join(outdir, colname + '_galaxy.pdf')
    plt.savefig(outpath)

    # Plot each key's data in separate files
    for key in keys:
        df_key = pd.DataFrame({'x': adata[adata.obs[colname] == key].obsm[embed][:, 0],
                               'y': adata[adata.obs[colname] == key].obsm[embed][:, 1]})
        fig, ax = plt.subplots(figsize=(8, 8))
        g = sns.kdeplot(x=df_key.x, y=df_key.y, n_levels=n_levels, cmap=cmap, ax=ax, shade=True)
        g.scatter(x=df_key.x, y=df_key.y, s=10, c='w', marker='.', linewidths=0, edgecolors='none')
        ax.set_xlim(global_bounds['x'])
        ax.set_ylim(global_bounds['y'])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('on')
        ax.set_title(key)
        ax.set_facecolor(cmap(1.0 / (n_levels * 2)))

        plt.tight_layout()
        outpath = os.path.join(outdir, key + '_galaxy.svg')  # Save as SVG
        plt.savefig(outpath)
        plt.close(fig)
    plt.close('all')


def main():
    args = getParser().parse_args()
    infile = args.input
    colname = args.colname
    embed = 'X_' + args.embed

    adata = sc.read_h5ad(infile)
    outdir = os.path.dirname(infile)
    galaxy_plot(adata, outdir, colname=colname, embed=embed)


if __name__ == '__main__':
    main()
