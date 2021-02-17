# Reproducing results for "Comparing N-mixture models and GLMMs for relative abundance estimation in eBird"
Reproducible code for paper comparing GLMM and N-mixture outcomes from the eBird dataset.

### How to use this repository

This repository contains code for processing eBird code, retrieving data for environmental 
covariates and running the analysis in the associated paper.

If you encounter any problems or have any questions, feel free to contact the corresponding 
author Ben Goldstein.

Due to copyright, this repository does **not** include the eBird dataset, nor
the LANDFIRE landscape covariate dataset. To use these scripts, these files
must be downloaded and added to the `data/` directory. Instructions for filepathing can
be found in the files where these inputs are read.

Once that's been done, numbered files can be run in sequence. The following scripts are included:

00_process_raw_eBird.R - converts the .txt eBird file into a CA-level relational database
01_prepare_eBird_data.R - filters eBird data for quality, chooses species-subregion datasets (SSRs)
02_prepare_covariate_data.R - associates checklists with spatial covariates (climate, landcover)
03_fit_one_ssr_final.R - contains the functions called for fitting models. This is the bulk of the work
04_batch_many_ssr.R - batch-runs the functions in 03 over all SSRs in parallel
05_visualize_model_results.R - Creates figures for primary results
06_post_tests.R - contains code defining post-tests (goodness-of-fit, stability, autocorrelation)
07_batch_posttests.R - runs the posttests over all SSRs
08_visualize_model_results.R - visualize posttest results

The following auxilliary scripts are also included:
model_code.R - contains code defining NIMBLE models
nmixgof_manual.R - support code for computing RQ residuals. Code adapted from Knape et al. 2018 (DOI: 10.1111/2041-210X.13062)
nmixgof_manual.cpp - c++ code for above
read_results_helper_file.R - contains code for translating output .RDS files in a summary df
