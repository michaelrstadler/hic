#!/usr/bin/env python

"""Foobar.py: Description of what foobar does."""

__author__      = "Michael Stadler"
__copyright__   = "Copyright 2022, California, USA, Planet Earth"

import numpy as np
from sklearn.decomposition import PCA
from flymovie import mesh_like
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
def norm_hic_matrix(mat):
    """Apply vanilla normalization to a matrix.
    
    Args:
        mat: ndarray, Hi-C matrix

    Returns:
        mat2: ndarray, normalized HiC-matrix    
    """
    mat2 = mat.copy()
    mean_ = np.mean(mat)
    mat2 = mat2 / (mat.mean(axis=0) + 1) * mean_
    mat2 = mat2.T
    mat2 = mat2 / (mat.mean(axis=1) + 1) * mean_
    mat2 = mat2.T
    return mat2

#-----------------------------------------------------------------------
def distnorm(mat):
    """Normalize a Hi-C matrix for distance to diagonal.
    
    Mean signal at each distance from diagonal is calculated, each
    position is divided by that mean.

    Args:
        mat: ndarray, Hi-C matrix

    Returns:
        ndarray, distance-normalized HiC-matrix   
    """
    x,y = mesh_like(mat, 2)
    distmesh = y - x
    distmeans = np.zeros_like(mat)
    for n in range(np.max(distmesh)):
        distmean = np.mean(mat[distmesh == n])
        distmeans[distmesh == n] = distmean
        distmeans[distmesh == (-1 * n)] = distmean
    distmeans[distmeans == 0] = np.mean(distmeans) / 1000
    distnormed = (mat / distmeans) * 1000
    distnormed = matrix_trim(distnormed, 1, 99)
    return  distnormed

#-----------------------------------------------------------------------
def pca_matrix(mat, n_components=2):
    """Perform PCA on a Hi-C matrix.
    
    Matrix data is first centered and trimmed.

    Args:
        mat: ndarray, Hi-C matrix, likely distance-normalized

    Return:
        centered: ndarray, centered version of input matrix
        pca: pca object, fitted PCA
        X: ndarray, PCA weight matrix applied to centered matrix. Rows
            are rows(/columns) of matrix, columns are PCs
    
    """
    # Make and center covariance matrix.
    cov = np.cov(mat)
    centered = (cov - np.mean(cov)) / np.std(cov)
    # Apply PCA.
    pca = PCA(n_components=n_components)
    pca.fit(centered)
    print(pca.explained_variance_ratio_)
    X = pca.transform(centered)
    return centered, pca, X

def matrix_trim(mat, pctlow=1, pcthigh=99):
    upper_cutoff = np.percentile(mat, pcthigh)
    lower_cutoff = np.percentile(mat, pctlow)
    mat_trimmed = mat.copy()
    mat_trimmed[mat > upper_cutoff] = upper_cutoff
    mat_trimmed[mat < lower_cutoff] = lower_cutoff
    return mat_trimmed
#-----------------------------------------------------------------------
def plot_compartments(cov, pca_X, vmax=None, filename=None, invert=False):
    if vmax == None:
        vmax = np.percentile(cov, 99)
    comp = np.linalg.eig(cov)[1][:,0]
    #comp = cov[:,0]
    if invert:
        comp = -1 * comp
    fig = plt.figure(constrained_layout=True, figsize=(10,10))
    subfigs = fig.subfigures(2, 1, wspace=0.07, height_ratios=[1, 6.])
    axs0 = subfigs[0].subplots()
    axs0.plot(comp, linewidth=2)
    #axs0.set_ylim(-1,1)
    axs0.margins(x=0)
    axs0.set_yticklabels([])
    axs0.set_xticklabels([])
    axs0.set_xticks([])
    axs1 = subfigs[1].subplots()
    axs1.imshow(matrix_trim(cov), vmax=vmax, aspect="auto", cmap='cividis')
    axs1.set_yticklabels([])
    axs1.margins(x=0)

    #plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, dpi=300)