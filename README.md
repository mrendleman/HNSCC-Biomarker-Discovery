# HNSCC-Biomarker-Discovery
A repository that contains some code related to my Dissertation, focused on feature engineering in large oncology data.


## DAE_preproc.R
Pre- and post-processing of DAE-constructed features.

## TCGA_data_split.R
Holdout data selection using a continuous stratification approach.

## data_split_simulations.R
Monte carlo simulations for empirical validation of the continuous stratification approach; including generation of simulated mRNAseq data, simulated outcomes using Weibull survival time generators, and evaluation of extra-sample performance estimation approaches with multiple survival model types and metrics.

This script is designed to run in parallel on SGE, so it is called with a command-line parameter indicating whether to: a) generate simulated data, b) run a specific simulation and write the results to a file, or c) aggregate output files of all completed simulations for analysis.

## install_bioc_rtcga.R
Script for downloading TCGA-HNSC data.
