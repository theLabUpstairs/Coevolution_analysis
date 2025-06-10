# Mantel Test Analysis

## Description

This directory contains an R-based workflow for performing Mantel tests between
breviate phylogenetic distance matrices and beta diversity matrices of
Arcobacteraceae species in each microbial community, in order to assess
potential coevolutionary patterns between breviates and Arcobacteraceae.

## Contents

- `01_do_all.R`: Main script to load data, run Mantel tests, and generate figures.
- `Mantel_test.Rproj`: RStudio project file.
- `figures/`: Output folder for PCoA plots.
- `_raw/`: Directory for raw or input data.

## Required Packages
- `tidyverse`
- `vegan`
- `ape`

## How to run

1. Open the project in RStudio via `Mantel_test.Rproj`.
2. Run `01_do_all.R` step by step or source the entire script.
