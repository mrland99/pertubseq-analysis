# Driver code for generating cluster maps for uprPerturb single-cell experiments
# Copyright (C) 2019  Max Land

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import os
import pandas as pd
import matplotlib
from cell_population import CellPopulation
from expression_normalization import z_normalize_expression, strip_low_expression
from clustermap_generator import clusterRandomSubset, clusterAll
from differential_expression import find_noisy_genes

pd.options.display.float_format = '{:.4f}'.format
matplotlib.rcParams['pdf.fonttype'] = 42

""" ONLY RUN ONCE. If you have already generated HDF file, you can just read in HDF. Much faster """

# read in rawMatrix data, filter it, and store it in HDF5 format. It is much
# faster then matrix market format
path = os.path.join(os.getcwd(), 'uprPerturbFiles')
pop = CellPopulation.from_file(path, filtered=False, raw_umi_threshold=2000)
# remove multiplets
pop = pop.subpopulation(cells='single_cell')
# For memory efficiency, remove genes with no counts from expression matrix
strip_low_expression(pop, threshold=0)
# store it as HDF
hdfPath = os.path.join(os.getcwd(), 'upr_sub_stripped_population.hdf')
pop.to_hdf(hdfPath, store_normalized_matrix=True)

""" OTHERWISE, load from hdf file """
pop = CellPopulation.from_hdf(hdfPath)

""" Generate expression matrix of differentially expressed genes"""
# combine cell samples by guide_target
mean_pop = pop.average('guide_target', verbose=True)

# find and subset differentially expressed genes
noisy_genes = find_noisy_genes(mean_pop)
mean_pop.matrix = z_normalize_expression(mean_pop, scale_by_total=True)
diff_expressed_matrix = mean_pop.where(genes=noisy_genes, gene_names=True)

# normalize expression distribution by Z-scoring

""" Generate expression matrix of normal pop """
# Normalize expression distribution by Z-scoring
pop.matrix = z_normalize_expression(pop, scale_by_total=True)
# store pop normalized matrix for future analysis if desired
pop_normalized_matrix = pop.matrix
# convert Enseml ids into human readable gene names
pop_normalized_matrix.columns = pop.gene_names(pop_normalized_matrix.columns)

""" Perform cluster analysis """
clusterAll(diff_expressed_matrix, method='pearson')
clusterAll(diff_expressed_matrix, method='euclidean')
# Example: cluster random subset of matrix and display heatmap with dendrogram
clusterRandomSubset(pop_normalized_matrix, 20000, 100)