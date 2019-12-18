# Driver code for generating cluster maps for epistasis single-cell experiments
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
import matplotlib
import pandas as pd
from cell_population import CellPopulation
from clustermap_generator import clusterRandomSubset
from expression_normalization import strip_low_expression, normalize_matrix_to_control

pd.options.display.float_format = '{:.4f}'.format
matplotlib.rcParams['pdf.fonttype'] = 42

""" ONLY RUN ONCE. If you have already generated HDF file, you can just read in HDF. Much faster """

# read in rawMatrix data, filter it, and store it in HDF5 format. It is much
# faster then matrix market format
path = os.path.join(os.getcwd(), 'epistasisExpFiles')
pop = CellPopulation.from_file(path, filtered=False, raw_umi_threshold=2000)
# remove multiplets
pop = pop.subpopulation(cells='single_cell')
# For memory efficiency, remove genes with no counts from expression matrix
strip_low_expression(pop, threshold=0)
# store it as HDF
pop.to_hdf(os.getcwd() + '/epi_sub_stripped_population.hdf', store_normalized_matrix=True)

""" OTHERWISE, load from hdf file """
pop = CellPopulation.from_hdf(os.getcwd() + '/epi_sub_stripped_population.hdf')

"""Normalize expression distribution by control group"""

# select all control cell samples
controlMatrix = pop.where(cells='guide_target == "3x"')
popNormMatrix = normalize_matrix_to_control(pop.matrix, controlMatrix)

# convert Enseml ids into human readable gene names
popNormMatrix.columns = pop.gene_names(popNormMatrix.columns)

# remove columns (genes) that contain nans
popNormMatrix = popNormMatrix.dropna(axis=1)

# EXAMPLE: cluster random subset of matrix and display heatmap with dendrogram
# There are many things you can do, such as average by guide_targets before clustering (see uprExperimentAnalysis)
clusterRandomSubset(popNormMatrix, 13000, 100)
