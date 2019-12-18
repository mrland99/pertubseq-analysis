# Perturbseq
Cluster analysis of gene expression data from large scale mammalian unfolded protein response Perturb-seq data.

This repository contains code reproducing single-cell analysis from:

Adamson, B., Norman, T.M., Jost, M., Parnas, O., Regev, A., Weissman, J.S., "A Multiplexed Single-Cell CRISPR Screening Platform Enables Systematic Dissection of the Unfolded Protein Response", Cell, 2016. DOI: https://doi.org/10.1016/j.cell.2016.11.048.

This version contains Perturbseq library code provided by @thomasmaxwellnorman. These library codes are useful processing single-cell expression matrix data. The actual code used to read in the data and perform the clustering is written by me (@mrland99). I have code performing clustering on two experiments: "UPR Epistasis Experiment" and "UPR Perturb-Seq Experiment". A demo can be found in `uprExperimentAnalysisDemo.ipynb`.

To run the code, you need to download the sequencing data from:
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406677 (Epistasis Experiment)
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681 (UPR Perturb-seq Experiment)

Only the outputs from `cell ranger` are necessary. There should be four files per experiment: `barcodes.tsv`, `genes.tsv`, `matrix.mtx`, `raw_cell_identities.csv`.  I recommend storing them in folders named `uprPerturbFiles` and `epistasisExpFiles`. You might have different structure, but it should be easy to adjust the code in `cell_population.py` and maybe some others as necessary.
