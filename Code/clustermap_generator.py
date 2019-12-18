# Perturbseq library for generating cluster maps for single-cell experiments
# Copyright (C) 2019  Max Land

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from scipy.spatial import distance_matrix
import pandas as pd
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
import seaborn as sns
import matplotlib

pd.options.display.float_format = '{:.4f}'.format
matplotlib.rcParams['pdf.fonttype'] = 42


def clusterAll(normalized_matrix, method='euclidean'):
    """ Performs hierchical clustering on an normalized matrix. Default distance metric is euclidean.

            Args: expressionMatrix = (pandas dataframe) Normalized expression matrix
                    method = (String) distance metric used to create distance matrix"""
    matrix = normalized_matrix

    # check if matrix contains any nans
    if matrix.isnull().values.any():
        print('Warning: Expression Matrix contains nan values')

    # custom color of heatmap
    color = sns.diverging_palette(240, 10, n=7)

    if method == 'pearson':
        # convert matrix to pandas dataframe
        df = pd.DataFrame(matrix.transpose().values)
        # calculate Pearson distance matrix
        dist_matrix = df.corr('pearson')
        dist_matrix = pd.DataFrame(dist_matrix.values, index=matrix.index, columns=matrix.index)
    else:
        # calculate Euclidean distance matrix
        dist_matrix = pd.DataFrame(distance_matrix(matrix.values, matrix.values), index=matrix.index, columns=matrix.index)
        # plot cluster heatmap
        linkage = hc.linkage(sp.distance.squareform(dist_matrix.values), method='average')
        sns.clustermap(dist_matrix, row_linkage=linkage, col_linkage=linkage, cmap=color)

    # custom color for heatmap to make it more aesthetic
    sns.clustermap(dist_matrix, cmap=color)


def clusterRandomSubset(normalized_matrix, num_row, num_col, method='euclidean'):
    """ Performs hierchical clustering on a random subset of normalized matrix. Default distance metric is euclidean.

        Args: expressionMatrix = (pandas dataframe) Normalized expression matrix to be subsetted
                num_row = (int) number of rows to be included in subset (in this case cell samples)
                num_col = (int) number of columns to be included in subset (in this case genes)
                method = (String) distance metric used to create distance matrix"""

    matrix = normalized_matrix

    # check if matrix contains any nans
    if matrix.isnull().values.any():
        print('Warning: Expression Matrix contains nan values')

    # subset number rows
    matrix = matrix.sample(n=num_row)

    # subset number of columns
    matrix = matrix.transpose().sample(n=num_col)

    # custom color of heatmap
    color = sns.diverging_palette(240, 10, n=7)

    if method == 'pearson':
        # convert matrix to pandas dataframe
        df = pd.DataFrame(matrix.transpose().values)
        # calculate Pearson distance matrix
        dist_matrix = df.corr('pearson')
        dist_matrix = pd.DataFrame(dist_matrix.values, index=matrix.index, columns=matrix.index)
    else:
        # calculate Euclidean distance matrix
        dist_matrix = pd.DataFrame(distance_matrix(matrix.values, matrix.values), index=matrix.index, columns=matrix.index)
        # plot cluster heatmap
        linkage = hc.linkage(sp.distance.squareform(dist_matrix.values), method='average')
        sns.clustermap(dist_matrix, row_linkage=linkage, col_linkage=linkage, cmap=color)

    # custom color for heatmap to make it more aesthetic
    sns.clustermap(dist_matrix, cmap=color)
